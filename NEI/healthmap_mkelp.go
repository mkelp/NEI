package main

import (
	"fmt"
	"image/color"
	"log"
	"os"
	"runtime"

	"net/http"
	_ "net/http/pprof" // Performance profiling

	"bitbucket.org/ctessum/aqhealth"
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

const totalPop = "TotalPop"

//Currently not used in this code... yet
var popTypes = []string{totalPop, "WhiteNoLat", "Black",
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
	concentrationholder := make(map[string][]float64) //creates outside data holder for concentrations
	recordholder := make(map[string]aep.Record)       //creates outside data holder for records
	//pols := make(map[string]int)										// pollutant names for csv
	//dims := make(map[string]unit.Dimensions)				// pollutant dimensions for ccsv

	//routine for processing NEI files
	for scc, fileTemplates := range c.NEIFiles {
		log.Printf("Started SCC %s", scc)
		r, err2 := aep.NewEmissionsReader(c.PolsToKeep, aep.Annually, aep.Ton)
		if err2 != nil {
			panic(err2)
		}
		r.Group = scc
		var files []*aep.InventoryFile
		for _, filetemplate := range fileTemplates {
			var tempFiles []*aep.InventoryFile
			tempFiles, err = r.OpenFilesFromTemplate(filetemplate)
			if err != nil {
				panic(err)
			}
			files = append(files, tempFiles...)
		}
		data, _, err4 := r.ReadFiles(files, nil)
		if err4 != nil {
			panic(err4)
		}
		for _, rec := range data {
			// runtime.GOMAXPROCS(8)

			scc := rec.GetSCC()

			//mkelp, filter by scc code to reduce run time!!!
			//Can uncomment below if you know the SCC codes you are intertested in
			// if scc >= "2310000100" {
			// 	continue
			// }

			//this routine calculates the total change in PM2.5 concentration caused by the
			// emissions in this record.
			pm := pmtotal(rec, c)
			if len(pm) == 0 {
				continue
			}
			if _, ok := concentrationholder[scc]; !ok {
				concentrationholder[scc] = pm
				recordholder[scc] = rec
			} else {
				floats.Add(concentrationholder[scc], pm) //adds concentrations to data holder
				recordholder[scc].CombineEmissions(rec)  // aggregates all emissions to data holder
			}
			for _, f := range files {
				// Close files.
				f.ReadSeeker.(*os.File).Close()
			}
		}
	}

	//FOR MAP
	//file path to scc descriptions, very nice to do this so that we can link the SCC code to its very descript config title
	q, err := os.Open("/home/mkelp/nei2011Dir/ge_dat/smkreport/sccdesc_pf31_05jan2015_v22.txt")
	sccDesc, err := aep.SCCDescription(q)

	//HEALTH IMPACT Calculation
	for scc := range recordholder {
		pm25Concentration := concentrationholder[scc]

		// Calculate health impacts of the emissions in this record.
		healthImpacts := make([]float64, len(pm25Concentration))
		pop := c.pop[totalPop] // map[string][]float64
		//popSum := floats.Sum(pop)
		//you can manipulate this however you would like, such as a populatation weighted health impact,
		//or, you can simply make a map of concentrations if you would like
		for i, conc := range pm25Concentration {
			rr := aqhealth.RRpm25Linear(conc)
			healthImpacts[i] = aqhealth.Deaths(rr, pop[i], c.mr[i])
			//pm25Concentration[i] = (conc * pop[i]) / (popSum)
		}

		//mk := floats.Sum(healthImpacts)

		// conc, err := c.sr.Source(pollutant, sourceIndex, layer)
		// if err != nil {
		// 	panic(err)
		// }
		// fmt.Println("Sum of all pollution concentrations: ", floats.Sum(conc))

		//MKELP
		// create color map.

		var states []geom.Geom
		states = getShapes("cb_2015_us_state_500k")
		cmap := carto.NewColorMap(carto.LinCutoff)

		//cmap.AddArray(healthImpacts)							//map of health function
		cmap.AddArray(pm25Concentration) // map of pollutant concentrations
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

			//color := cmap.GetColor(healthImpacts[i])					// map of health function
			colorx := cmap.GetColor(pm25Concentration[i]) // map of concentrations
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
		hey, err := os.Create(scc + "_sr_examplemk_mkelpmap.png")
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
		log.Printf("Getting concentrations for pol %s (%d total)", pol, len(srg.Elements))
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
func pmtotal(rec aep.Record, c *Config) []float64 {

	//2 d array? spatialize each pol

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

	// Calculate the total change in PM2.5 concentration caused by the
	// emissions in this record.
	pm25Concentration := make([]float64, len(c.gridCells))
	for pollutant, value := range totals {
		inmapPol, ok := pollutantCrosswalk[pollutant.Name]
		if !ok {
			panic(fmt.Errorf("missing pollutant %s", pollutant.Name))
		}
		floats.AddScaled(pm25Concentration, value.Value(),
			concentrationSurrogates[inmapPol])
	}
	return pm25Concentration
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
