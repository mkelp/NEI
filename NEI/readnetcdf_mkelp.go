package main

import (

	//"github.com/ctessum/unit"

	"encoding/csv"
	"fmt"
	"image/color"
	"log"
	"os"
	"sort"
	"strconv"

	"bitbucket.org/ctessum/aqhealth"
	"bitbucket.org/ctessum/cdf"
	"bitbucket.org/ctessum/sparse"
	"bitbucket.org/ctessum/sr/sr"
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

//this script takes a netcdf file created by another script and creates maps of the US for either
//emissions or health impacts per pollutant

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

const (
	totalPop   = "TotalPop"
	whitenolat = "WhiteNoLat"
	black      = "Black"
	native     = "Native"
	asian      = "Asian"
	other      = "Other"
	latino     = "Latino"
	poverty    = "Poverty"
	twoxpov    = "TwoXPov"
)

var popTypes = []string{totalPop, whitenolat, black,
	native, asian, other, latino, poverty, twoxpov}

const numCells int = 50105

func main() {

	// Create a new configuration variable
	c := new(Config)

	// Open sr file, update filepath to sr file
	e, err := os.Open("/home/mkelp/srMatrix/sr.ncf")
	if err != nil {
		panic(err)
	}
	// create sr object from file
	c.sr, err = sr.New(e)
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
	gridCells, err := c.sr.GridCells()
	if err != nil {
		panic(err)
	}
	c.gridCells = gridCells.Cells

	//open netcdf file
	ff, err := os.Open("/home/mkelp/work/src/NEI/spatialResults_mkelp_test_ptoil.ncf")
	if err != nil {
		panic(err)
	}

	f, err := cdf.Open(ff)
	if err != nil {
		panic(err)
	}
	pol := f.Header.Variables()
	fmt.Println("list of pollutants", pol)
	concentrationholder := make(map[string]map[string][]float64)
	var value []float64

	//open csv file for the scc index
	ff2, err := os.Open("/home/mkelp/work/src/NEI/index_numbers_test_ptoil.csv")
	if err != nil {
		panic(err)
	}
	r := csv.NewReader(ff2)
	lines, err := r.ReadAll()
	if err != nil {
		panic(err)
	}

	//fmt.Println(lines)
	sccIndices := make(map[string]int)
	for _, line := range lines {
		//sector := (line[0])
		scc := (line[1])
		//index2 := (line[2])

		fmt.Println("scc", scc)
		index2, err2 := strconv.ParseInt(line[2], 10, 32)
		if err2 != nil {
			panic(err2)
		}
		index := int(index2)

		sccIndices[scc] = index
	}
	ff2.Close()

	var data []float64

	//	totalholder := make(map[string][]float64)
	//smallholder := make(map[string][]float64)
	for scc := range sccIndices {
		// start := []int{i, 0}
		// end := []int{i, numCells - 1}
		concentrationholder[scc] = make(map[string][]float64)
		for _, v := range pol {
			r := f.Reader(v, nil, nil)

			//Make start and end values 'nil' to read in all the data at once

			buf := r.Zero(-1)
			//fmt.Println(buf)
			n, err2 := r.Read(buf)
			if err2 != nil {
				panic(err2)
			}

			// fmt.Println("scc and index", scc, i)
			//add i as data object in loop declaration to view index from csv
			data = buf.([]float64)
			fmt.Println("sum of data", floats.Sum(data))

			value = append(value, float64(n))
			// fmt.Println("sum of inmap rows", floats.Sum(value))

			concentrationholder[scc][v] = data
		}
	}

	// Re-Write data to a CSV

	//new csv produced
	w, err := os.Create("mkelp_readin_netcdf_20161118.csv")
	if err != nil {
		panic(err)
	}

	ww := csv.NewWriter(w)

	q, err := os.Open("/home/mkelp/nei2011Dir/ge_dat/smkreport/sccdesc_pf31_05jan2015_v22.txt")
	if err != nil {
		panic(err)
	}
	sccDesc, err := aep.SCCDescription(q)
	if err != nil {
		panic(err)
	}

	// Write header for CSV
	row0 := []string{}

	header := []string{}

	sort.Strings(pol)
	fmt.Println(pol)

	for _, populationtype := range popTypes {
		for _, poltitle := range pol { //ranging over poltitle might be more reliable than pollutant
			header = append(header, fmt.Sprintf("Deaths from %s to %s", poltitle, populationtype))
		}
	}

	// row0 = append([]string{"Sector: "}, "SCC: ", "SCC Descriptions: ")
	row0 = append([]string{"SCC: "}, "SCC Descriptions: ")

	//row0 = append(row0, pol...)
	row0 = append(row0, header...)

	//Get emissions labels
	// for _, mk := range data {
	// 	for i, pol := range row0 {
	// 		if mk.dims[pol] != nil {
	// 			row0[i] += fmt.Sprintf(" (%s)", mk.dims[pol].String())
	// 			row0[i] += " Emissions"
	// 		}
	// 	}
	// }

	// //writes the header for the csv file
	err5 := ww.Write(row0)
	if err5 != nil {
		panic(err5)
	}
	fmt.Println("writing the header")

	//data writing coded

	for scc, concentration := range concentrationholder {

		//	fmt.Println("storedata hello channel3")

		healthimpacts := make(map[string]map[string][]float64)
		for poptype, pop := range c.pop {
			healthimpacts[poptype] = make(map[string][]float64)
			for pol, conc := range concentration {
				h := make([]float64, len(pop))
				for i, w := range conc {
					rr := aqhealth.RRpm25Linear(-w)
					h[i] = -aqhealth.Deaths(rr, pop[i], c.mr[i])
				}
				healthimpacts[poptype][pol] = h
			}
		}

		//Data to CSV

		var linedeaths []string
		for _, pop := range popTypes {
			for _, p := range pol {
				linedeaths = append(linedeaths, fmt.Sprintf("%g", floats.Sum(healthimpacts[pop][p])))
			}
		}

		// totals := rec.Totals()
		//fmt.Println("totals 2", mk.sector, scc, mk.totals)
		// fmt.Println("scc codes", mk.scc)

		line := make([]string, len(row0))
		//line[0] = sector //sort?
		line[0] = scc //sort?
		line[1] = sccDesc[scc]

		// for i, pol := range pol {
		// 	for pollutant, value := range mk.totals[scc] {
		// 		if pollutant.Name == pol {
		// 			line[i+3] = fmt.Sprintf("%g", value.Value())
		// 		}
		// 	}
		// }
		line = line[:len(line)-len(header)]

		line = append(line, linedeaths...)

		err5 := ww.Write(line)
		if err5 != nil {
			panic(err5)
		}
	}

	ww.Flush()
	w.Close()
	fmt.Println("storeData done")

	//FOR MAP
	//file path to scc descriptions, very nice to do this so that we can link the SCC code to its very descript config title
	// q, err := os.Open("/home/mkelp/nei2011Dir/ge_dat/smkreport/sccdesc_pf31_05jan2015_v22.txt")
	// sccDesc, err := aep.SCCDescription(q)

	//MKELP
	// create color map.

	// These loops will give the health effects of each pollutant to each race Group
	//See healthmap_mkelp.go for pollutant concentration maps

	//HEALTH IMPACT Calculation
	for scc, polConc := range concentrationholder {

		var states []geom.Geom
		states = getShapes("cb_2015_us_state_500k")
		cmap := carto.NewColorMap(carto.LinCutoff)

		//fmt.Println(concentration)
		// healthimpacts := make(map[string]map[string][]float64)
		// for poptype, pop := range c.pop {
		// 	healthimpacts[poptype] = make(map[string][]float64)
		// 	for pol, conc := range concentration {
		// 		h := make([]float64, len(pop))
		// 		//fmt.Println(conc)
		// 		for i, w := range conc {
		// 			rr := aqhealth.RRpm25Linear(-w)
		// 			h[i] = -aqhealth.Deaths(rr, pop[i], c.mr[i])
		// 		}
		// 		healthimpacts[poptype][pol] = h
		// 	}
		// }

		//NOTE: take out NEXT line for CONCENTRATIONS OF EACH POLLUTANT
		//for _, pop := range popTypes {
		for p, conc := range polConc {

			// if floats.Sum(conc) == float64(0) {
			// 	continue
			// }

			cmap.AddArray(conc) // map of pollutant concentrations
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
				colorx := cmap.GetColor(conc[i]) // map of concentrations
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

			hey, err := os.Create("Emissions  for " + p + " in " + scc + ".png")
			//hey, err := os.Create("20161031mkelp to " + pop + "by " + p + scc + ".png")
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

//setupSpatialProcessor reads in the necessary information to initialize
//a processor for spatializing emissions, and then does so.
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

func checkPol(pol string) error {
	found := false
	for _, p := range []string{"soa", "primarypm2_5", "pnh4", "pso4", "pno3"} {
		if pol == p {
			found = true
			break
		}
	}
	if !found {
		return fmt.Errorf("unknown pollutant %s", pol)
	}
	return nil
}
