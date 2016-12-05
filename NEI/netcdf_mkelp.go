package main

import (
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"runtime"
	"sort"
	"strconv"
	"sync"

	//"github.com/ctessum/unit"

	"bitbucket.org/ctessum/cdf"
	"bitbucket.org/ctessum/sparse"
	"bitbucket.org/ctessum/sr/sr"
	"github.com/BurntSushi/toml"
	"github.com/ctessum/aep"
	"github.com/ctessum/geom"
	"github.com/ctessum/geom/proj"
	"github.com/ctessum/unit"
	"github.com/gonum/floats"
)

//this script reads in NEI data, spatializes it and writes teh data to a netcdf file

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

	//define channel to send and recieve information for each sector
	dataChan := make(chan *dataHolder)

	var wgProcess sync.WaitGroup
	//var wgStore sync.WaitGroup
	var wgnetcdf sync.WaitGroup

	// Generates a list of all the sectors of interst!

	fmt.Println("length of sectors: ", len(c.NEIFiles))
	wgProcess.Add(len(c.NEIFiles)) // We need to wait for 4 sectors
	//wgStore.Add(1)
	wgnetcdf.Add(1)

	runtime.GOMAXPROCS(8)

	//go storeData(dataChan, c, &wgStore) // Process all of the sectors concurrently
	go spatialOutputter(dataChan, &wgnetcdf)

	for sector, fileTemplates := range c.NEIFiles {

		log.Printf("Started Sector %s", sector)

		go process(sector, c, fileTemplates, dataChan, &wgProcess)

	}
	fmt.Println("waiting for processing to finish")
	wgProcess.Wait() // wait for processing to finish

	fmt.Println("Sending signal that processing is finished")
	close(dataChan) // let the storage function know that processing is finished
	fmt.Println("Waiting for writing to finish")
	//	wgStore.Wait() // wait for the storage function to finish

	wgnetcdf.Wait()

	fmt.Println("Writing is finished")

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

// struct to store values from process function and then to pass it onto the store function without having to re-process each time

type dataHolder struct {
	sector              string
	concentrationholder map[string]map[string][]float64 //creates outside data holder
	recordholder        map[string]aep.Record           //creates outside data holder
	pols                map[string]int                  // pollutant names for csv
	dims                map[string]unit.Dimensions      // pollutant dimensions for csv
	poltitle            []string
	totals              map[string]map[aep.Pollutant]*unit.Unit
	scccodes            []string
	totalRows           int
}

func process(sector string, c *Config, fileTemplates []string, dataChan chan *dataHolder, wg *sync.WaitGroup) {
	var files []*aep.InventoryFile

	//try to take out all code above

	r, err2 := aep.NewEmissionsReader(c.PolsToKeep, aep.Annually, aep.Ton)
	if err2 != nil {
		panic(err2)
	}

	for _, filetemplate := range fileTemplates {
		var tempFiles []*aep.InventoryFile

		tempFiles, err := r.OpenFilesFromTemplate(filetemplate)
		if err != nil {
			panic(err)
		}
		files = append(files, tempFiles...)
	}
	data, _, err4 := r.ReadFiles(files, nil)
	if err4 != nil {
		panic(err4)
	}

	//initializes the variables in the struct dataHolder

	d := new(dataHolder)
	d.concentrationholder = make(map[string]map[string][]float64)
	d.pols = make(map[string]int)             // pollutant names for csv
	d.dims = make(map[string]unit.Dimensions) // pollutant dimensions for csv
	d.sector = sector
	d.totals = make(map[string]map[aep.Pollutant]*unit.Unit)
	d.poltitle = []string{}
	// d.totalRows = int
	sccx := []string{}
	fmt.Println("process break")

	for _, rec := range data {

		scc := rec.GetSCC()

		sccx = append(sccx, scc)

		pm := pmtotal(rec, c)
		if len(pm) == 0 {
			continue
		}

		recordTotals := rec.Totals()

		for pollutant, amount := range recordTotals {
			d.pols[pollutant.Name] = 0
			d.dims[pollutant.Name] = amount.Dimensions()
		}

		if _, ok := d.concentrationholder[scc]; !ok {
			d.concentrationholder[scc] = pm
			d.totals[scc] = recordTotals
		} else {
			//recordholder[scc].CombineEmissions(rec) // adds concentrations to the data holder
			for pol, subpm := range pm {
				if _, ok := d.concentrationholder[scc][pol]; !ok {
					d.concentrationholder[scc][pol] = subpm //initializes empty pollutants to have some length
				} else {
					floats.Add(d.concentrationholder[scc][pol], subpm) // aggregates all emissions to the data holder
				}
				d.poltitle = append(d.poltitle, pol)
			}
		}
		for pollutant, amount := range recordTotals {
			if _, ok := d.totals[scc][pollutant]; !ok {
				d.totals[scc][pollutant] = amount.Clone()
			} else {
				d.totals[scc][pollutant].Add(amount)
			}
		}

	}

	//generates a list of unique SCC codes

	tempscc := make(map[string]string)
	for _, pol := range sccx {
		tempscc[pol] = ""
	}

	sccnum := make([]string, 0, len(tempscc))
	for p := range tempscc {
		sccnum = append(sccnum, p)
	}

	d.scccodes = sccnum
	//fmt.Println("scc codes", d.scccodes)
	d.totalRows = len(d.concentrationholder)
	fmt.Println("total Rows", d.totalRows)
	fmt.Println("len of scc", len(d.scccodes))
	for _, f := range files {
		// Close files.
		f.ReadSeeker.(*os.File).Close()
	}
	fmt.Println("process number2", len(d.concentrationholder))
	dataChan <- d
	fmt.Println("process number3", len(d.concentrationholder))

	wg.Done()
}

//creates a netcdf file with concentration for each pollutant for each of the 50105 inmap grid Cells
//in addition, a csv file of SCC codes with an asscoiated index number are produced to be used in the readnetcdf script

func spatialOutputter(dataChan chan *dataHolder, wg *sync.WaitGroup) {

	//	initialize data holder from different sectors and recieve all the data to use in this function

	data := make(map[string]*dataHolder)
	for d := range dataChan {
		data[d.sector] = d
	}

	const inmapNrows = 50105

	//creates a string of unique pollutants
	var totalRows2 int
	tempPols := make(map[string]string)
	for _, d := range data {
		totalRows2 = d.totalRows
		for _, pol := range d.poltitle {
			tempPols[pol] = ""
		}
	}

	//fmt.Println("tempPols", tempPols)
	poltitle := make([]string, 0, len(tempPols))

	for p := range tempPols {
		poltitle = append(poltitle, p)
	}

	sort.Strings(poltitle)
	//fmt.Println("poltitle", poltitle)

	//creates netcdf file
	f, err := os.Create("/home/mkelp/work/src/NEI/spatialResults_mkelp_test_ptegu.ncf")

	h := cdf.NewHeader([]string{"scc", "cell"}, []int{totalRows2, inmapNrows})
	fmt.Println("Netcdf header", []int{totalRows2, inmapNrows})
	//
	for _, w := range poltitle {
		h.AddVariable(w, []string{"scc", "cell"}, []float64{0})
	}

	h.Define()

	cf, err := cdf.Create(f, h)
	if err != nil {
		panic(err)
	}

	//netcdf := make(map[string][]float64)
	//totalnetcdf := make(map[string][]float64)

	//assigns concentration values for all grid cells per pollutant
	//index := 0

	// m, err := os.Create("index_numbers_test.csv")
	// if err != nil {
	// 	panic(err)
	// }
	// line := make([]string, len(strconv.Itoa(totalRows)))

	//NEW CSV writer

	m, err := os.Create("index_numbers_test_ptegu.csv")
	if err != nil {
		panic(err)
	}
	mm := csv.NewWriter(m)
	line := make([]string, totalRows2)

	index := 0
	index2 := 0

	for sector, mk := range data {
		for scc, pollutants := range mk.concentrationholder {
			// fmt.Println(sector, "sector")
			// fmt.Println(scc, "scc")
			// fmt.Println(index, "index")
			// fmt.Println(index2, "index2")

			line[0] = sector
			line[1] = scc
			line[2] = strconv.Itoa(index2)

			err5 := mm.Write(line)
			if err5 != nil {
				panic(err5)
			}

			for pol, concArray := range pollutants {

				if len(concArray) != inmapNrows {
					panic(fmt.Errorf("incorrect number of rows: %d", len(concArray)))
				}
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
		index++
	}

	mm.Flush()
	m.Close()
	wg.Done()
	if err := cdf.UpdateNumRecs(f); err != nil {
		panic(err)
	}
	f.Close()

}

// writes pollutant data to netcdf file
// for _, code := range scc {
// 	start := []int{sccindex[code], 0}
// 	end := []int{sccindex[code], inmapNrows}
// 	for _, d := range poltitle {
// 		if len(netcdf[d]) != inmapNrows {
// 			panic(fmt.Errorf("incorrect number of rows: %d", len(netcdf[d])))
// 		}
// 		w := cf.Writer(d, start, end)
// 		_, err := w.Write(netcdf[d])
// 		if err != nil {
// 			panic(err)
// 		}
// 	}
// }
// wg.Done()
// f.Close()
