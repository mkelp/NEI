package main

import (
	"fmt"
	"os"

	"bitbucket.org/ctessum/aqhealth"
	"bitbucket.org/ctessum/sr/sr"
	"github.com/ctessum/aep"
	"github.com/ctessum/geom"
	"github.com/ctessum/geom/carto"
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

var popTypes = []string{totalPop, "WhiteNoLat", "Black",
	"Native", "Asian", "Other", "Latino", "Poverty", "TwoXPov"}

func main() {

	// Open sr file
	f, err := os.Open("/home/mkelp/srMatrix/sr.ncf")
	if err != nil {
		panic(err)
	}
	// create sr object from file
	data, err := sr.New(f)
	if err != nil {
		panic(err)
	}
	// load population

	pop, err := data.Population("TotalPop")
	if err != nil {
		panic(err)
	}

	mr, err := data.MortalityBaseline()
	if err != nil {
		panic(err)
	}

	// load grid cells.
	cells, err := data.GridCells()
	if err != nil {
		panic(err)
	}

	// c.mr, err = c.sr.MortalityBaseline()
	// if err != nil {
	// 	panic(err)
	// }
	// c.pop = make(map[string][]float64)
	// for _, p := range popTypes {
	// 	c.pop[p], err = c.sr.Population(p)
	// 	if err != nil {
	// 		panic(err)
	// 	}

	const (
		x         = 0              // meters; center of domain
		y         = 0              // meters; center of domain
		r         = 500000000      // meters
		layer     = 2              // elevated release
		pollutant = "primarypm2_5" // emissions of SOx causing pso4 concentrations
	)
	sourceIndex := cells.Index(x, y)
	cellsOfInterest := cells.IndicesWithinRadius(x, y, r)

	conc, err := data.Source(pollutant, sourceIndex, layer)
	if err != nil {
		panic(err)
	}
	fmt.Println("Sum of all pollution concentrations: ", floats.Sum(conc))

	healthImpacts := make([]float64, len(conc))
	for i, c := range conc {
		rr := aqhealth.RRpm25Linear(c)
		healthImpacts[i] = aqhealth.Deaths(rr, pop[i], mr[i])
	}
	fmt.Println("Sum of all Deaths: ", floats.Sum(healthImpacts))

	//MKELP
	// create color map.
	cmap := carto.NewColorMap(carto.Linear)
	cmap.AddArray(healthImpacts)
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

	// create legend
	cmap.Legend(&legendc, "Population-weighted concentration (μg/m3/ton)")

	// create map
	c := carto.NewCanvas(600000, -600000, 600000, -600000, plotc)
	for _, i := range cellsOfInterest {
		// only plot cells within 500,000 m of source.
		color := cmap.GetColor(healthImpacts[i])
		ls := draw.LineStyle{
			Color: color,
			Width: vg.Points(0.5),
		}
		err = c.DrawVector(cells.Cells[i], color, ls, draw.GlyphStyle{})
		if err != nil {
			panic(err)
		}
	}

	// save image to file
	hey, err := os.Create("sr_examplemk.png")
	if err != nil {
		panic(err)
	}
	_, err = vgimg.PngCanvas{Canvas: ic}.WriteTo(hey)
	if err != nil {
		panic(err)
	}

	fmt.Println("Output image can be viewed at https://bytebucket.org/ctessum/sr/raw/master/sr/sr_example.png")
}
