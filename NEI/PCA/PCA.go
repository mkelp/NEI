package main

import (
	"encoding/csv"
	"fmt"
	"os"
	"strconv"

	"github.com/gonum/matrix/mat64"

	"bitbucket.org/ctessum/cdf"
	"github.com/gonum/stat"
)

func main() {

	const numCells int = 50105

	ff, err := os.Open("/home/mkelp/work/src/NEI/spatialResults_mkelp_test_afdust_to_npoilgas.ncf")
	f, err := cdf.Open(ff)
	if err != nil {
		panic(err)
	}
	pol := f.Header.Variables()
	fmt.Println("list of pollutants", pol)
	var value []float64

	//open csv file for the scc index
	ff2, err := os.Open("/home/mkelp/work/src/NEI/index_numbers_test_afdust_to_npoilgas.csv")
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
		if v == "PM25-PRI" {
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
	hi := mat64.NewDense(50105, index0, data2)
	vecs, vars, ok := stat.PrincipalComponents(hi, nil)
	if !ok {
		return
	}
	fmt.Printf("variances = %.4f\n\n", vars)

	// Project the data onto the first 2 principal components.
	k := 5
	var proj mat64.Dense
	proj.Mul(hi, vecs.View(0, 0, index0, k))

	//fmt.Printf("proj = %.4f", mat64.Formatted(&proj, mat64.Prefix("       ")))
}
