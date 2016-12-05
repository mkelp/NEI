// You can edit this code!
// Click here and start typing.
package main

import (
	"fmt"
	"sort"
)

func main() {

	populations := []string{"total", "white"}

	population := map[string][]float64{
		"total": []float64{0, 1, 0, 5},
		"white": []float64{0, 0, 0, 3},
	}

	totals := map[string]float64{
		"pollutant 1": 1,
		"pollutant 2": 6,
		"pollutant 3": 10,
	}

	crosswalk := map[string]string{
		"pollutant 1": "a",
		"pollutant 2": "b",
		"pollutant 3": "c",
	}

	concentrationSurrogates := map[string][]float64{
		"a": []float64{1, 2, 3, 8},
		"b": []float64{34, 12, 33, 12},
		"c": []float64{11, 21, 16, 39},
	}

	concentration := make(map[string][]float64)

	for pol, totalVal := range totals {
		inmapPol, ok := crosswalk[pol]
		if !ok {
			panic("xxxxxxx")
		}
		concSrg, ok := concentrationSurrogates[inmapPol]
		if !ok {
			panic("xxxxxxx")
		}
		conc := make([]float64, len(concSrg))
		for i, v := range concSrg {
			conc[i] = v * totalVal
		}
		concentration[pol] = conc
	}
	fmt.Println(concentration)

	health := make(map[string]map[string][]float64)

	for popType, pop := range population {
		health[popType] = make(map[string][]float64)
		for pol, conc := range concentration {
			h := make([]float64, len(pop))
			for i, c := range conc {
				h[i] = c * pop[i]
			}
			health[popType][pol] = h
		}
	}
	fmt.Println(health)

	var pollutants []string
	for pol := range concentration {
		pollutants = append(pollutants, pol)
	}
	sort.Strings(pollutants)

	var header []string
	for _, pol := range pollutants {
		for _, popType := range populations {
			header = append(header, fmt.Sprintf("deaths to %s from %s", popType, pol))
		}
	}
	fmt.Println(header)
	var line []string
	for _, p := range pollutants {
		for _, popType := range populations {
			line = append(line, fmt.Sprintf("%g", health[popType][p]))
		}
	}
	fmt.Println(line)
}
