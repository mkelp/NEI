package main

import (
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"strconv"

	//"github.com/ctessum/unit"

	"bitbucket.org/ctessum/cdf"
	"bitbucket.org/ctessum/sparse"
	"bitbucket.org/ctessum/sr/sr"
	"github.com/BurntSushi/toml"
	"github.com/ctessum/aep"
	"github.com/ctessum/geom"
	"github.com/ctessum/geom/proj"
	"github.com/gonum/floats"
)

//this script produces a CSV of a specified sector's SCC code, total emissions per pollutant (kg), and health impact caused by each pollutant to a
//specifed population type, without using channels

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

const inmapNrows = 50105

func main() {

	// go func() {
	// 	// runtime.GOMAXPROCS(8)
	//
	// 	log.Println(http.ListenAndServe("localhost:6060", nil))
	// }()

	// Create a new configuration variable
	c := new(Config)

	// Open the configuration file
	r, err := os.Open("cstref.toml")
	if err != nil {
		panic(err) // Deal with any errors that may have come up
	}

	// Read the configuration file into the configuration variable.
	if _, err = toml.DecodeReader(r, c); err != nil {
		panic(err)
	}

	f, err := os.Open(c.SRFile)
	if err != nil {
		panic(err)
	}
	c.sr, err = sr.New(f)
	if err != nil {
		print(err)
	}
	c.mr, err = c.sr.MortalityBaseline()
	if err != nil {
		panic(err)
	}
	c.pop = make(map[string][]float64)
	for _, p := range popTypes {
		c.pop[p], err = c.sr.Population(p)
		if err != nil {
			panic(err)
		}
	}
	gridCells, err := c.sr.GridCells()
	if err != nil {
		panic(err)
	}
	c.gridCells = gridCells.Cells

	// Initialize a spatial processor. We're not actually using this now because
	// we're just counting the total number of records,
	// but it will be important in the future.
	if c.sp, err = c.setupSpatialProcessor(); err != nil {
		panic(err)
	}

	// Now we will read through all of the NEI files.

	//totalReport := new(aep.InventoryReport) // initiliazes the files

	concentrationholder := make(map[string]map[string][]float64) //creates outside data holder
	recordholder := make(map[string]aep.Record)                  //creates outside data holder
	pols := []string{}                                           // pollutant names for csv
	totalRows2 := 0
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
					pols = append(pols, pol)
					if _, ok := concentrationholder[scc][pol]; !ok {
						concentrationholder[scc][pol] = subpm //initializes empty pollutants to have some length
					} else {
						floats.Add(concentrationholder[scc][pol], subpm) // aggregates all emissions to the data holder
					}
				}
			}

			// 	concentrationholder[scc] = make(map[string][]float64)
			// 	concentrationholder[scc] = pm
			//
			// 	for pol, conc := range pm {
			// 		h := make([]float64, len(conc))
			// 		for i, w := range conc {
			// 			h[i] = (h[i] + w)
			// 		}
			// 		concentrationholder[scc][pol] = h
			// 	}

		}

		for _, f := range files {
			// Close files.
			f.ReadSeeker.(*os.File).Close()
		}
	}

	// for _, rec := range recordholder {
	// 	//fmt.Println("scc: ", scc)
	// 	recordTotals := rec.Totals() // map[Pollutant]Value
	// 	for pollutant, amount := range recordTotals {
	// 		pols[pollutant.Name] = 0
	// 		//fmt.Println(pollutant.Name, amount.Value(), amount.Dimensions())
	//
	// 		//Dimensions check
	// 		if _, ok := dims[pollutant.Name]; !ok {
	// 			dims[pollutant.Name] = amount.Dimensions()
	// 		}
	// 		//else { // SOMETHING WRONG, dim not the same
	// 		// 	if !amount.Dimensions().Matches(d) {
	// 		// 		panic(fmt.Errorf("dimensions mismatch: '%v' != '%v'", amount.Dimensions(), d))
	// 		// 	}
	// 		// }
	// 	}
	// 	//fmt.Println("totals: ", rec.Totals())
	// }

	//NEW CSV writer

	// ww := csv.NewWriter(w)
	//
	// q, err := os.Open("/home/mkelp/nei2011Dir/ge_dat/smkreport/sccdesc_pf31_05jan2015_v22.txt")
	// if err != nil {
	// 	panic(err)
	// }
	// sccDesc, err := aep.SCCDescription(q)
	// if err != nil {
	// 	panic(err)
	// }
	//
	// // Write header for CSV
	//
	// row0 := []string{}
	// pollutants := []string{}
	//
	// header := []string{}
	//
	// for scc := range recordholder {
	//
	// 	concentrationx := concentrationholder[scc]
	// 	concentrationkey := make(map[string]string)
	// 	// //used in loop to allow only unique pollutants to show up in header
	// 	// //assign a generic key so that redundant pollutants are overwritten
	// 	//
	// 	// // for pol := range concentrationx {
	// 	// // 	pollutants = append(pollutants, pol)
	// 	// // }
	// 	for pol := range concentrationx {
	// 		concentrationkey[pol] = ""
	// 	}
	// 	pollutants = make([]string, 0, len(concentrationx))
	// 	for polx := range concentrationkey {
	// 		pollutants = append(pollutants, polx)
	// 	}

	//
	// fmt.Println(pollutants)

	//fmt.Println(header)

	tempPols := make(map[string]string)
	for _, pol := range pols {
		tempPols[pol] = ""
	}
	poltitle := []string{}
	//fmt.Println("tempPols", tempPols)

	for p := range tempPols {
		poltitle = append(poltitle, p)
	}

	fmt.Println("poltitle", poltitle)
	// fmt.Println("pollutantkey", polkey)
	//fmt.Println("poltitle", poltitle)
	// fmt.Println("pollutants", pollutants)
	// fmt.Println("concentrationkey", concentrationkey)

	// for scc, rec := range recordholder {
	//
	// 	concentration := concentrationholder[scc]
	//
	// 	healthimpacts := make(map[string]map[string][]float64)
	// 	for poptype, pop := range c.pop {
	// 		healthimpacts[poptype] = make(map[string][]float64)
	// 		for pol, conc := range concentration {
	// 			h := make([]float64, len(pop))
	// 			for i, w := range conc {
	// 				rr := aqhealth.RRpm25Linear(-w)
	// 				h[i] = -aqhealth.Deaths(rr, pop[i], c.mr[i])
	// 			}
	// 			healthimpacts[poptype][pol] = h
	// 		}
	// 	}

	totalRows2 = len(concentrationholder)

	p, err := os.Create("/home/mkelp/work/src/NEI/spatialResults_mkelp_loop.ncf")

	h := cdf.NewHeader([]string{"scc", "cell"}, []int{totalRows2, inmapNrows})
	fmt.Println("Netcdf header", []int{totalRows2, inmapNrows})
	//
	for _, w := range poltitle {
		h.AddVariable(w, []string{"scc", "cell"}, []float64{0})
	}

	h.Define()

	cf, err := cdf.Create(p, h)
	if err != nil {
		panic(err)
	}

	//NEW CSV writer

	m, err := os.Create("/home/mkelp/work/src/NEI/index_numbers_loop.csv")
	if err != nil {
		panic(err)
	}
	mm := csv.NewWriter(m)
	line := make([]string, totalRows2)

	// index := 0
	index2 := 0

	for scc, pollutants := range concentrationholder {

		line[0] = scc
		line[1] = strconv.Itoa(index2)

		err5 := mm.Write(line)
		if err5 != nil {
			panic(err5)
		}

		for pol, concArray := range pollutants {

			if len(concArray) != inmapNrows {
				panic(fmt.Errorf("incorrect number of rows: %d", len(concArray)))
			}

			fmt.Println("sum of data at netcdf writing", floats.Sum(concArray))

			start := []int{index2, 0}
			end := []int{index2, inmapNrows}

			// fmt.Println("len of concArray", len(concArray))
			// fmt.Println("pol", pol)
			w := cf.Writer(pol, start, end)
			_, err := w.Write(concArray)
			if err != nil {
				panic(err)
			}
		}

		index2++
	}
	// index++

	// fmt.Println("index- number of sectors", index)
	fmt.Println("index2- number of SCC codes", index2)

	mm.Flush()
	m.Close()
	if err := cdf.UpdateNumRecs(p); err != nil {
		panic(err)
	}
	f.Close()

}

//If no errors, should ouput a CSV file with emissions and deaths for each SCC code

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
	"PM2_5":    "primarypm2_5",
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
