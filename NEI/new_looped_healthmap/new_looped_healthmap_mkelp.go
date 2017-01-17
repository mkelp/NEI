package main

import (
	"fmt"
	"image/color"
	"log"
	"os"
	"runtime"
	"sort"

	"net/http"
	_ "net/http/pprof" // Performance profiling

	"bitbucket.org/ctessum/sparse"
	"bitbucket.org/ctessum/sr/sr"
	"github.com/BurntSushi/toml"
	"github.com/ctessum/aep"
	"github.com/ctessum/geom"
	"github.com/ctessum/geom/carto"
	"github.com/ctessum/geom/encoding/shp"
	"github.com/ctessum/geom/proj"
	"github.com/gonum/floats"
	"github.com/gonum/plot/vg"
	"github.com/gonum/plot/vg/draw"
	"github.com/gonum/plot/vg/vgimg"
)

//script creates maps for each pollutant and its subsequent impact on health to specified populations

// Config holds configuration information.
type Config struct {

	// SRFile gives the location of the InMAP SR matrix data file.
	SRFile string

	// NEIFiles lists National Emissions Inventory emissions files to use
	// for making SCC-based spatial surrogates. The file names can include
	// environment variables. The format is map[sector name][list of files].
	NEIFiles map[string][]string

	// PolsToKeep lists pollutants from the NEI that should be kept.
	PolsToKeep map[string]*aep.PolHolder

	// NEIData holds information needed for processing the NEI.
	NEIData struct {
		// SrgSpec gives the location of the surrogate specification file.
		SrgSpec string

		// SrgShapefileDirectory gives the location of the directory holding
		// the shapefiles used for creating spatial surrogates.
		SrgShapefileDirectory string

		// SCCExactMatch specifies whether SCC codes must match exactly when processing
		// emissions.
		SCCExactMatch bool

		// GridRef specifies the locations of the spatial surrogate gridding
		// reference files used for processing the NEI.
		GridRef []string

		// OutputSR specifies the output spatial reference in Proj4 format.
		OutputSR string

		// InputSR specifies the input spatial reference in Proj4 format.
		InputSR string

		// SimplifyTolerance is the tolerance for simplifying spatial surrogate
		// geometry, in units of OutputSR.
		SimplifyTolerance float64
	}

	// SpatialCache specifies the location for storing spatial emissions
	// data for quick access. If this is left empty, no cache will be used.
	SpatialCache string

	// MaxCacheEntries specifies the maximum number of emissions and concentrations
	// surrogates to hold in a memory cache. Larger numbers can result in faster
	// processing but increased memory usage.
	MaxCacheEntries int

	sp *aep.SpatialProcessor

	sr        *sr.SR
	mr        []float64            // Baseline mortality rate.
	pop       map[string][]float64 // Population by demographic type.
	gridCells []geom.Polygonal     // InMAP grid cell geometry

}

//Currently not used in this code... yet
var popTypes = []string{"TotalPop", "WhiteNoLat", "Black",
	"Native", "Asian", "Other", "Latino", "Poverty", "TwoXPov"}

func main() {

	go func() {
		runtime.GOMAXPROCS(1)

		log.Println(http.ListenAndServe("localhost:6060", nil))
	}()

	// Create a new configuration variable
	c := new(Config)

	// Open the configuration file, filepath to config file
	r, err := os.Open("cstref.toml")
	if err != nil {
		panic(err) // Deal with any errors that may have come up
	}

	// Read the configuration file into the configuration variable.
	if _, err = toml.DecodeReader(r, c); err != nil {
		panic(err)
	}

	// Open sr file, update filepath to sr file
	f, err := os.Open("/home/mkelp/srMatrix/sr.ncf")
	if err != nil {
		panic(err)
	}
	// create sr object from file
	c.sr, err = sr.New(f)
	if err != nil {
		panic(err)
	}
	// load population
	c.pop = make(map[string][]float64)
	for _, p := range popTypes {
		c.pop[p], err = c.sr.Population(p)
		if err != nil {
			panic(err)
		}
	}

	c.mr, err = c.sr.MortalityBaseline()
	if err != nil {
		panic(err)
	}

	// load grid cells.
	gridCells, err := c.sr.GridCells()
	if err != nil {
		panic(err)
	}

	c.gridCells = gridCells.Cells

	if c.sp, err = c.setupSpatialProcessor(); err != nil {
		panic(err)
	}

	//DATA Manipulation
	concentrationholder := make(map[string]map[string][]float64) //creates outside data holder
	recordholder := make(map[string]aep.Record)                  //creates outside data holder

	for sector, fileTemplates := range c.NEIFiles {
		log.Printf("Started Sector %s", sector)
		r, err2 := aep.NewEmissionsReader(c.PolsToKeep, aep.Annually, aep.Ton)
		if err2 != nil {
			panic(err2)
		}
		r.Group = sector
		var files []*aep.InventoryFile
		for _, filetemplate := range fileTemplates {
			var tempFiles []*aep.InventoryFile
			tempFiles, err = r.OpenFilesFromTemplate(filetemplate)
			if err != nil {
				panic(err)
			}
			fmt.Println(tempFiles)
			files = append(files, tempFiles...)
		}
		fmt.Println(files)

		data, _, err4 := r.ReadFiles(files, nil)
		if err4 != nil {
			panic(err4)
		}
		for _, rec := range data {
			scc := rec.GetSCC()

			// if scc >= "2310000230" {
			// 	continue
			// }
			//
			// if scc >= "2801700003" {
			// 	continue
			// }

			pm := pmtotal(rec, c)
			if len(pm) == 0 {
				continue
			}

			if _, ok := concentrationholder[scc]; !ok {
				concentrationholder[scc] = pm
				recordholder[scc] = rec
			} else {
				recordholder[scc].CombineEmissions(rec) // adds concentrations to the data holder
				for pol, subpm := range pm {
					if _, ok := concentrationholder[scc][pol]; !ok {
						concentrationholder[scc][pol] = subpm //initializes empty pollutants to have some length
					} else {
						floats.Add(concentrationholder[scc][pol], subpm) // aggregates all emissions to the data holder
					}
				}
			}

		}
		for _, f := range files {
			// Close files.
			f.ReadSeeker.(*os.File).Close()
		}
	}
	pollutants := []string{}

	//Define Pollutants array to be analyzed
	for scc := range recordholder {

		concentrationx := concentrationholder[scc]
		concentrationkey := make(map[string]string)
		// //used in loop to allow only unique pollutants to show up in header
		// //assign a generic key so that redundant pollutants are overwritten
		//
		// // for pol := range concentrationx {
		// // 	pollutants = append(pollutants, pol)
		// // }
		for pol := range concentrationx {
			concentrationkey[pol] = ""
		}
		pollutants = make([]string, 0, len(concentrationx))
		for polx := range concentrationkey {
			pollutants = append(pollutants, polx)
		}
	}
	sort.Strings(pollutants)

	for _, x := range concentrationholder {
		for _, w := range x {
			fmt.Println("Sum of data at withdraw", floats.Sum(w))
		}
	}

	//FOR MAP
	//file path to scc descriptions, very nice to do this so that we can link the SCC code to its very descript config title
	q, err := os.Open("/home/mkelp/nei2011Dir/ge_dat/smkreport/sccdesc_pf31_05jan2015_v22.txt")
	sccDesc, err := aep.SCCDescription(q)

	//HEALTH IMPACT Calculation
	for scc := range recordholder {

		concentration := concentrationholder[scc]

		// healthimpacts := make(map[string]map[string][]float64)
		// for poptype, pop := range c.pop {
		// 	healthimpacts[poptype] = make(map[string][]float64)
		// 	for pol, conc := range concentration {
		// 		h := make([]float64, len(pop))
		// 		for i, w := range conc {
		// 			rr := aqhealth.RRpm25Linear(-w)
		// 			h[i] = -aqhealth.Deaths(rr, pop[i], c.mr[i])
		// 		}
		// 		healthimpacts[poptype][pol] = h
		// 	}
		// }

		//MKELP
		// create color map.

		var states []geom.Geom
		states = getShapes("cb_2015_us_state_500k")
		cmap := carto.NewColorMap(carto.LinCutoff)

		// These loops will give the health effects of each pollutant to each race Group
		//See healthmap_mkelp.go for pollutant concentration maps

		//NOTE: take out NEXT line for CONCENTRATIONS OF EACH POLLUTANT
		//for _, pop := range popTypes {
		for _, p := range pollutants {
			cmap.AddArray(concentration[p]) // map of pollutant concentrations
			//cmap.AddArray(healthimpacts[pop][p]) // map of healthimpacts
			cmap.Set()

			const (
				figWidth     = 4 * vg.Inch
				figHeight    = 3.1 * vg.Inch
				legendHeight = 0.3 * vg.Inch
			)

			// create image
			ic := vgimg.New(figWidth, figHeight)
			dc := draw.New(ic)
			legendc := draw.Crop(dc, 0, 0, 0, legendHeight-figHeight)
			plotc := draw.Crop(dc, 0, 0, legendHeight, 0)

			// creates legend with SCC description
			cmap.Legend(&legendc, sccDesc[scc])

			// create map of the USA
			v := carto.NewCanvas(1944000, -2088000, 2592000, -2736000, plotc)

			for i, cell := range gridCells.Cells {

				// NOTE: change for concentrations
				colorx := cmap.GetColor(concentration[p][i]) // map of concentrations
				//colorx := cmap.GetColor(healthimpacts[pop][p][i]) // map of healthimpacts
				ls := draw.LineStyle{
					Color: colorx,
					Width: vg.Points(0.5),
				}
				err = v.DrawVector(cell, colorx, ls, draw.GlyphStyle{})
				if err != nil {
					panic(err)
				}
			}
			mk := color.NRGBA{0, 0, 0, 0}
			sl := draw.LineStyle{
				Color: color.Black,
				Width: 0.1 * vg.Millimeter,
			}
			for _, g := range states {
				v.DrawVector(g, mk, sl, draw.GlyphStyle{})
			}

			// saves image to file name, outputs SCC code and whatever title you would like to give it
			// NOTE: Change for concentrations
			//hey, err := os.Create("20161118 " + pop + "by " + p + scc + ".png")
			hey, err := os.Create("Rail,Marine 20161208" + p + scc + ".png")

			if err != nil {
				panic(err)
			}
			_, err = vgimg.PngCanvas{Canvas: ic}.WriteTo(hey)
			if err != nil {
				panic(err)
			}

			fmt.Println(sccDesc[scc], " Map made")
		}
	}
}

//}

// setupSpatialProcessor reads in the necessary information to initialize
// a processor for spatializing emissions, and then does so.
func (c *Config) setupSpatialProcessor() (*aep.SpatialProcessor, error) {
	f, err := os.Open(os.ExpandEnv(c.NEIData.SrgSpec))
	if err != nil {
		return nil, err
	}
	srgSpecs, err := aep.ReadSrgSpec(f, os.ExpandEnv(c.NEIData.SrgShapefileDirectory), c.NEIData.SCCExactMatch)
	if err != nil {
		return nil, err
	}
	if err = f.Close(); err != nil {
		return nil, err
	}

	var gridRef *aep.GridRef
	for _, gf := range c.NEIData.GridRef {
		f, err = os.Open(os.ExpandEnv(gf))
		if err != nil {
			return nil, err
		}
		gridRefTemp, err2 := aep.ReadGridRef(f)
		if err2 != nil {
			return nil, err2
		}
		if err = f.Close(); err != nil {
			return nil, err
		}
		if gridRef == nil {
			gridRef = gridRefTemp
		} else {
			err = gridRef.Merge(*gridRefTemp)
			if err != nil {
				return nil, err
			}
		}
	}

	outSR, err := proj.Parse(os.ExpandEnv(c.NEIData.OutputSR))
	if err != nil {
		return nil, err
	}
	inSR, err := proj.Parse(os.ExpandEnv(c.NEIData.InputSR))
	if err != nil {
		return nil, err
	}
	gridCells, err := c.sr.GridCells()
	if err != nil {
		return nil, err
	}
	grid, err := aep.NewGridIrregular("InMAP", gridCells.Cells, outSR, outSR)
	if err != nil {
		return nil, err
	}
	matchFullSCC := false
	sp := aep.NewSpatialProcessor(srgSpecs, []*aep.GridDef{grid}, gridRef, inSR, matchFullSCC)
	sp.DiskCachePath = c.SpatialCache
	sp.SimplifyTolerance = c.NEIData.SimplifyTolerance
	return sp, nil
}

// inmapResultVars is a list of variables that we are keeping from the
// InMAP output.
var inmapResultVars = []string{"soa", "primarypm2_5", "pnh4", "pso4", "pno3"}

var pollutantCrosswalk = map[string]string{
	"VOC":      "soa",
	"PM25-PRI": "primarypm2_5",
	"PM25":     "primarypm2_5",
	"NOX":      "pno3",
	"NH3":      "pnh4",
	"SO2":      "pso4",
}

// runSR creates an InMAP surrogate using the SR matrix.
func (c *Config) runSR(srg *sparse.SparseArray) map[string][]float64 {

	const (
		// SR inputs are in units of ton/year, so convert kg/year to ton/year.
		tonPerKg = 0.00110231
	)

	o := make(map[string][]float64)
	for _, pol := range inmapResultVars {
		oo := make([]float64, len(c.gridCells))
		//log.Printf("Getting concentrations for pol %s (%d total)", pol, len(srg.Elements))
		for i, val := range srg.Elements {
			conc, err := c.sr.Source(pol, i, 0) // TODO: account for plume rise here.
			if err != nil {
				panic(err)
			}
			for i, cc := range conc {
				oo[i] += cc * val * tonPerKg
			}
		}
		o[pol] = oo
	}
	return o
}

// pm2.5 data withdraw

func pmtotal(rec aep.Record, c *Config) map[string][]float64 {
	// Get normalized spatial locations of emissions for this record.
	spatialSurrogate, _, inGrid, err := rec.Spatialize(c.sp, 0)
	if err != nil {
		panic(err)
	}
	if !inGrid {
		return nil
	}

	// Get normalized spatial locations of the changes in concentrations
	// caused by emissions in this record
	concentrationSurrogates := c.runSR(spatialSurrogate)

	totals := rec.Totals()

	concentration := make(map[string][]float64)

	for pollutant, value := range totals {
		inmapPol, ok := pollutantCrosswalk[pollutant.Name]
		if !ok {
			panic(fmt.Errorf("missing pollutant %s", pollutant))
		}
		concSrg, ok := concentrationSurrogates[inmapPol]
		if !ok {
			panic("concentrationSurrogates error")
		}
		conc := make([]float64, len(concSrg))
		for i, v := range concSrg {
			conc[i] = v * value.Value()
		}
		concentration[pollutant.Name] = conc
	}
	//fmt.Println(concentration)

	return concentration
}

func getShapes(filename string) []geom.Geom {
	dec, err := shp.NewDecoder(filename + ".shp")
	handle(err)

	// This is the spatial projection that the SR matrix is in.
	sr1, err := proj.Parse("+proj=lcc +lat_1=33.000000 +lat_2=45.000000 +lat_0=40.000000 +lon_0=-97.000000 +x_0=0 +y_0=0 +a=6370997.000000 +b=6370997.000000 +to_meter=1")
	handle(err)

	// This is the projection of the shapefile.
	sr2, err := dec.SR()
	handle(err)
	ct, err := sr2.NewTransform(sr1)
	handle(err)

	// Read all the shapes from the shapefile and store them.
	var counties []geom.Geom
	for {
		var d struct{ geom.Geom }
		if more := dec.DecodeRow(&d); !more {
			break
		}
		d.Geom, err = d.Geom.Transform(ct)
		handle(err)
		counties = append(counties, d.Geom)
	}
	handle(dec.Error())
	return counties
}

func handle(err error) {
	if err != nil {
		panic(err)
	}
}
