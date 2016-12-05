package main

import (
	"fmt"
	"log"
	"os"
	"runtime"
	"sync"

	"bitbucket.org/ctessum/sr/sr"
	"github.com/BurntSushi/toml"
	"github.com/ctessum/aep"
	"github.com/ctessum/geom"
	"github.com/ctessum/geom/proj"
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

func main() {
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

	// Initialize a spatial processor. We're not actually using this now because
	// we're just counting the total number of records,
	// but it will be important in the future.
	//if c.sp, err = c.setupSpatialProcessor(); err != nil {
	//	panic(err)
	//}
	// Now we will read through all of the NEI files.

	// totalRecords := 0
	var wgProcess sync.WaitGroup
	var wgStore sync.WaitGroup

	fmt.Println("length of sectors: ", len(c.NEIFiles))
	wgProcess.Add(len(c.NEIFiles)) // We need to wait for 4 sectors
	dataChan := make(chan string)
	// totalReport := new(aep.InventoryReport)
	//
	wgStore.Add(1)

	runtime.GOMAXPROCS(1)
	go storeData(dataChan, &wgStore) // Process all of the sectors concurrently

	for sector, fileTemplates := range c.NEIFiles {
		log.Printf("Started Sector %s", sector)

		go process(sector, c, fileTemplates, dataChan, &wgProcess)
	}
	// 	go process(sector, fileTemplates, dataChan, wg)
	//
	// 	log.Printf("Started sector %s", sector)
	// 	r, err2 := aep.NewEmissionsReader(c.PolsToKeep, aep.Annually, aep.Ton)
	// 	if err2 != nil {
	// 		panic(err2)
	// 	}
	// r.Group = sector
	// var files []*aep.InventoryFile
	// for _, filetemplate := range fileTemplates {
	// 	var tempFiles []*aep.InventoryFile
	// 	tempFiles, err = r.OpenFilesFromTemplate(filetemplate)
	// 	if err != nil {
	// 		panic(err)
	// 	}
	// 	files = append(files, tempFiles...)
	// }
	// data, report, err3 := r.ReadFiles(files, nil)
	// if err != nil {
	// 	panic(err3)
	fmt.Println("waiting for processing to finish")
	wgProcess.Wait() // wait for processing to finish

	fmt.Println("Sending signal that processing is finsihed")
	close(dataChan) // let the storage function know that processing is finished
	fmt.Println("Waiting for writing to finish")
	wgStore.Wait() // wait for the storage function to finish
	fmt.Println("Writing is fininsehd")
}

// log.Printf("Finished sector %s", sector)
// totalReport.AddData(report.Files...)
// totalRecords += len(data)

// for _, f := range files {
// 	// Close files.
// 	f.ReadSeeker.(*os.File).Close()
// }
// }

// Write out the totals to a file.
// w, err := os.Create("totals.csv")
// if err != nil {
// 	panic(err)
// }
// ww := csv.NewWriter(w)
// ww.WriteAll(totalReport.TotalsTable())
// ww.Flush()
// w.Close()
//
// log.Printf("############## The total number of records is %d.", totalRecords)

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
	sp := aep.NewSpatialProcessor(srgSpecs, []*aep.GridDef{grid}, gridRef, inSR, true)
	sp.DiskCachePath = c.SpatialCache
	sp.SimplifyTolerance = c.NEIData.SimplifyTolerance
	return sp, nil
}

func process(sector string, c *Config, fileTemplates []string, dataChan chan string, wg *sync.WaitGroup) {
	var files []*aep.InventoryFile

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
	println("troubleshoot", len(data))

	d := sector
	dataChan <- d
}

// This function receives all of the sector and stores it.
func storeData(dataChan chan string, wg *sync.WaitGroup) {
	// c := new(Config)
	fmt.Println("storedata started")

	for mk := range dataChan {

		fmt.Println("storedata hello channel2")

		fmt.Println("storeData", mk)
	}
	wg.Done()
}
