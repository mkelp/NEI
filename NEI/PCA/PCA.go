package main

import (
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"strconv"

	"github.com/gonum/matrix/mat64"
	"github.com/gonum/stat"

	"bitbucket.org/ctessum/cdf"
	"bitbucket.org/ctessum/sr/sr"
	"github.com/BurntSushi/toml"
	"github.com/ctessum/aep"
	"github.com/ctessum/geom"
	"github.com/gonum/floats"
)

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
const numCells int = 50105

var popTypes = []string{totalPop, whitenolat, black,
	native, asian, other, latino, poverty, twoxpov}

func main() {
	c := new(Config)

	// Open the configuration file
	j, err := os.Open("cstref.toml")
	if err != nil {
		panic(err) // Deal with any errors that may have come up
	}

	// Read the configuration file into the configuration variable.
	if _, err = toml.DecodeReader(j, c); err != nil {
		panic(err)
	}

	g, err := os.Open(c.SRFile)
	if err != nil {
		panic(err)
	}
	c.sr, err = sr.New(g)
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

	NCDFPOP := make([]float64, numCells)

	for sector, fileTemplates := range c.NEIFiles {

		log.Printf("Started Sector %s", sector)
		j, err2 := aep.NewEmissionsReader(c.PolsToKeep, aep.Annually, aep.Ton)
		if err2 != nil {
			panic(err2)
		}
		//r.Group = sector
		var files []*aep.InventoryFile
		for _, filetemplate := range fileTemplates {
			var tempFiles []*aep.InventoryFile
			tempFiles, err = j.OpenFilesFromTemplate(filetemplate)
			if err != nil {
				panic(err)
			}
			fmt.Println(tempFiles)
			files = append(files, tempFiles...)
		}
		fmt.Println(files)

		for popt, w := range c.pop {
			if popt == totalPop {
				NCDFPOP = w
			} else {
				continue
			}
		}

	}

	fmt.Println("sum of pop", len(NCDFPOP), floats.Sum(NCDFPOP))

	ff, err := os.Open("/home/mkelp/work/src/NEI/spatialResults_mkelp_20161229FINAL.ncf")
	f, err := cdf.Open(ff)
	if err != nil {
		panic(err)
	}
	pol := f.Header.Variables()
	fmt.Println("list of pollutants", pol)
	var value []float64

	//open csv file for the scc index
	ff2, err := os.Open("/home/mkelp/work/src/NEI/index_numbers_20161229FINAL.csv")
	if err != nil {
		panic(err)
	}
	r := csv.NewReader(ff2)
	lines, err := r.ReadAll()
	if err != nil {
		panic(err)
	}
	index0 := 0
	sccIndices := make(map[string]int)
	for _, line := range lines {

		//sector := (line[0])
		scc := (line[1])
		//index2 := (line[2])

		// fmt.Println("scc", scc)
		index2, err2 := strconv.ParseInt(line[2], 10, 32)
		if err2 != nil {
			panic(err2)
		}
		index := int(index2)
		sccIndices[scc] = index
		index0++
	}
	ff2.Close()

	var data2 []float64
	// start := []int{i, 0}
	// end := []int{i, numCells - 1}

	for _, v := range pol {
		if v == "VOC" {
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
			data2 = buf.([]float64)
			// fmt.Println("sum of data", floats.Sum(data))

			value = append(value, float64(n))
			// fmt.Println("sum of inmap rows", floats.Sum(value))

		} else {
			continue
		}
	}
	fmt.Println("index0", index0)
	// NewDense(row, column)
	hi := mat64.NewDense(50105, index0, data2)

	rows, cols := hi.Dims()

	totalConc := make([]float64, rows)
	for r := 0; r < rows; r++ {
		var sum float64
		for c := 0; c < cols; c++ {
			sum += hi.At(r, c)
		}
		totalConc[r] = sum
	}
	fmt.Println("pop len", len(NCDFPOP))
	fmt.Println("conc len", len(totalConc))
	fmt.Println("rows", rows, "columns", cols)
	//make weights NCDFPOP

	floats.Mul(NCDFPOP, totalConc)
	vecs, vars, ok := stat.PrincipalComponents(hi, NCDFPOP)
	if !ok {
		panic("idnt' work")
	}
	fmt.Printf("variances = %.4f\n\n", vars)

	// Project the data onto the first k principal components.
	k := 7118
	var proj mat64.Dense
	proj.Mul(hi, vecs.View(0, 0, index0, k))

	fmt.Printf("proj = %.4f", mat64.Formatted(&proj, mat64.Prefix("       ")))

	//NEW UPDATE
	// pc := new(stat.PC)
	// ok := pc.PrincipalComponents(hi, totalConc)
	// if !ok {
	// 	panic("idnt' work")
	// }
	// vars := pc.Vars(nil)
}
