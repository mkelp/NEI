package main

import (
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"sort"

	"bitbucket.org/ctessum/aqhealth"
	"bitbucket.org/ctessum/sparse"
	"bitbucket.org/ctessum/sr/sr"
	"github.com/BurntSushi/toml"
	"github.com/ctessum/aep"
	"github.com/ctessum/geom"
	"github.com/ctessum/geom/proj"
	"github.com/ctessum/unit"
	"github.com/gonum/floats"
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

	concentrationholder := make(map[string][][]float64) //creates outside data holder
	recordholder := make(map[string]aep.Record)         //creates outside data holder
	pols := make(map[string]int)                        // pollutant names for csv
	dims := make(map[string]unit.Dimensions)            // pollutant dimensions for csv

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

			// if scc >= "2310000100" {
			// 	continue
			// }
			pm := pmtotal(rec, c) //outputs a 2D array

			if len(pm) == 0 {
				continue
			}
			if _, ok := concentrationholder[scc]; !ok {
				concentrationholder[scc] = pm
				recordholder[scc] = rec
			} else {
				recordholder[scc].CombineEmissions(rec) // adds concentrations to the data holder

				for i, subpm := range pm {
					floats.Add(concentrationholder[scc][i], subpm) // aggregates all emissions to the data holder
				}
			}
		}
		for _, f := range files {
			// Close files.
			f.ReadSeeker.(*os.File).Close()
		}
	}

	for _, rec := range recordholder {
		//fmt.Println("scc: ", scc)
		recordTotals := rec.Totals() // map[Pollutant]Value
		for pollutant, amount := range recordTotals {
			pols[pollutant.Name] = 0
			//fmt.Println(pollutant.Name, amount.Value(), amount.Dimensions())

			//Dimensions check
			if _, ok := dims[pollutant.Name]; !ok {
				dims[pollutant.Name] = amount.Dimensions()
				// } else { // SOMETHING WRONG, dim not the same
				// 	if !amount.Dimensions().Matches(d) {
				// 		panic(fmt.Errorf("dimensions mismatch: '%v' != '%v'", amount.Dimensions(), d))
				// 	}
				// }
			}
			//fmt.Println("totals: ", rec.Totals())
		}

		// Write out the totals to a csv file.
		w, err := os.Create("totalsmkelp_race_health.csv")
		if err != nil {
			panic(err)
		}

		//NEW CSV writer

		ww := csv.NewWriter(w)

		poltitle := []string{}
		for pol := range pols {
			poltitle = append(poltitle, pol)
		}
		sort.Strings(poltitle)

		row0 := append([]string{"SCC: "}, poltitle...)
		row0 = append(row0, "SCC_descriptions", "total[Emissions]", "total[NOx]", "total[PRI]", "total[SO2]", "total[VOC]",
			"tot_deathsE", "deaths_whiteE", "death_blackE", "death_nativeE", "death_asianE", "death_otherE", "death_latinoE", "death_povertyE", "death_twoxpovE",
			"tot_deathsnox", "deaths_whitenox", "death_blacknox", "death_nativenox", "death_asiannox", "death_othernox", "death_latinonox", "death_povertynox", "death_twoxpovnox",
			"tot_deathspri", "deaths_whitepri", "death_blackpri", "death_nativepri", "death_asianpri", "death_otherpri", "death_latinopri", "death_povertypri", "death_twoxpovpri",
			"tot_deathsvoc", "deaths_whitevoc", "death_blackvoc", "death_nativevoc", "death_asianvoc", "death_othervoc", "death_latinovoc", "death_povertyvoc", "death_twoxpovvoc",
			"tot_deathsso2", "deaths_whiteso2", "death_blackso2", "death_nativeso2", "death_asianso2", "death_otherso2", "death_latinoso2", "death_povertyso2", "death_twoxpovso2",
		)

		//links units to pollutants
		for i, pol := range row0 {
			if dims[pol] != nil {
				row0[i] += fmt.Sprintf(" (%s)", dims[pol].String())
			}
		}

		err5 := ww.Write(row0)
		if err5 != nil {
			panic(err5)
		}

		q, err := os.Open("/home/mkelp/nei2011Dir/ge_dat/smkreport/sccdesc_pf31_05jan2015_v22.txt")
		if err != nil {
			panic(err)
		}
		sccDesc, err := aep.SCCDescription(q)
		if err != nil {
			panic(err)
		}

		for scc, rec := range recordholder {
			pm25Concentration := concentrationholder[scc][0]
			noxConcentration := concentrationholder[scc][1]
			PRIConcentration := concentrationholder[scc][2]
			so2Concentration := concentrationholder[scc][3]
			vocConcentration := concentrationholder[scc][4]

			// Calculate health impacts of the emissions in this record.
			//Total PM2.5

			healthImpactspm := make([]float64, len(pm25Concentration))
			whitehealthImpactspm := make([]float64, len(pm25Concentration))
			blackhealthImpactspm := make([]float64, len(pm25Concentration))
			nativehealthImpactspm := make([]float64, len(pm25Concentration))
			asianealthImpactspm := make([]float64, len(pm25Concentration))
			otherhealthImpactspm := make([]float64, len(pm25Concentration))
			latinohealthImpactspm := make([]float64, len(pm25Concentration))
			povertyealthImpactspm := make([]float64, len(pm25Concentration))
			twoxpovalthImpactspm := make([]float64, len(pm25Concentration))

			pop := c.pop[totalPop]
			whitepop := c.pop[whitenolat]
			blackpop := c.pop[black]
			nativepop := c.pop[native]
			asianpop := c.pop[asian]
			otherpop := c.pop[other]
			latinopop := c.pop[latino]
			povertyop := c.pop[poverty]
			twoxpovpop := c.pop[twoxpov]

			for i, conc := range pm25Concentration {
				rr := aqhealth.RRpm25Linear(conc)
				healthImpactspm[i] = aqhealth.Deaths(rr, pop[i], c.mr[i])
				whitehealthImpactspm[i] = aqhealth.Deaths(rr, whitepop[i], c.mr[i])
				blackhealthImpactspm[i] = aqhealth.Deaths(rr, blackpop[i], c.mr[i])
				nativehealthImpactspm[i] = aqhealth.Deaths(rr, nativepop[i], c.mr[i])
				asianealthImpactspm[i] = aqhealth.Deaths(rr, asianpop[i], c.mr[i])
				otherhealthImpactspm[i] = aqhealth.Deaths(rr, otherpop[i], c.mr[i])
				latinohealthImpactspm[i] = aqhealth.Deaths(rr, latinopop[i], c.mr[i])
				povertyealthImpactspm[i] = aqhealth.Deaths(rr, povertyop[i], c.mr[i])
				twoxpovalthImpactspm[i] = aqhealth.Deaths(rr, twoxpovpop[i], c.mr[i])
			}

			//NOX health effects

			healthImpactsnox := make([]float64, len(noxConcentration))
			whitehealthImpactsnox := make([]float64, len(noxConcentration))
			blackhealthImpactsnox := make([]float64, len(noxConcentration))
			nativehealthImpactsnox := make([]float64, len(noxConcentration))
			asianealthImpactsnox := make([]float64, len(noxConcentration))
			otherhealthImpactsnox := make([]float64, len(noxConcentration))
			latinohealthImpactsnox := make([]float64, len(noxConcentration))
			povertyealthImpactsnox := make([]float64, len(noxConcentration))
			twoxpovalthImpactsnox := make([]float64, len(noxConcentration))

			for i, conc := range noxConcentration {
				rr := aqhealth.RRpm25Linear(conc)
				healthImpactsnox[i] = aqhealth.Deaths(rr, pop[i], c.mr[i])
				whitehealthImpactsnox[i] = aqhealth.Deaths(rr, whitepop[i], c.mr[i])
				blackhealthImpactsnox[i] = aqhealth.Deaths(rr, blackpop[i], c.mr[i])
				nativehealthImpactsnox[i] = aqhealth.Deaths(rr, nativepop[i], c.mr[i])
				asianealthImpactsnox[i] = aqhealth.Deaths(rr, asianpop[i], c.mr[i])
				otherhealthImpactsnox[i] = aqhealth.Deaths(rr, otherpop[i], c.mr[i])
				latinohealthImpactsnox[i] = aqhealth.Deaths(rr, latinopop[i], c.mr[i])
				povertyealthImpactsnox[i] = aqhealth.Deaths(rr, povertyop[i], c.mr[i])
				twoxpovalthImpactsnox[i] = aqhealth.Deaths(rr, twoxpovpop[i], c.mr[i])
			}

			//PRI health impacts
			healthImpactspri := make([]float64, len(PRIConcentration))
			whitehealthImpactspri := make([]float64, len(PRIConcentration))
			blackhealthImpactspri := make([]float64, len(PRIConcentration))
			nativehealthImpactspri := make([]float64, len(PRIConcentration))
			asianealthImpactspri := make([]float64, len(PRIConcentration))
			otherhealthImpactspri := make([]float64, len(PRIConcentration))
			latinohealthImpactspri := make([]float64, len(PRIConcentration))
			povertyealthImpactspri := make([]float64, len(PRIConcentration))
			twoxpovalthImpactspri := make([]float64, len(PRIConcentration))

			for i, conc := range noxConcentration {
				rr := aqhealth.RRpm25Linear(conc)
				healthImpactspri[i] = aqhealth.Deaths(rr, pop[i], c.mr[i])
				whitehealthImpactspri[i] = aqhealth.Deaths(rr, whitepop[i], c.mr[i])
				blackhealthImpactspri[i] = aqhealth.Deaths(rr, blackpop[i], c.mr[i])
				nativehealthImpactspri[i] = aqhealth.Deaths(rr, nativepop[i], c.mr[i])
				asianealthImpactspri[i] = aqhealth.Deaths(rr, asianpop[i], c.mr[i])
				otherhealthImpactspri[i] = aqhealth.Deaths(rr, otherpop[i], c.mr[i])
				latinohealthImpactspri[i] = aqhealth.Deaths(rr, latinopop[i], c.mr[i])
				povertyealthImpactspri[i] = aqhealth.Deaths(rr, povertyop[i], c.mr[i])
				twoxpovalthImpactspri[i] = aqhealth.Deaths(rr, twoxpovpop[i], c.mr[i])
			}

			//SO2 health impacts

			healthImpactsso2 := make([]float64, len(so2Concentration))
			whitehealthImpactsso2 := make([]float64, len(so2Concentration))
			blackhealthImpactsso2 := make([]float64, len(so2Concentration))
			nativehealthImpactsso2 := make([]float64, len(so2Concentration))
			asianealthImpactsso2 := make([]float64, len(so2Concentration))
			otherhealthImpactsso2 := make([]float64, len(so2Concentration))
			latinohealthImpactsso2 := make([]float64, len(so2Concentration))
			povertyealthImpactsso2 := make([]float64, len(so2Concentration))
			twoxpovalthImpactsso2 := make([]float64, len(so2Concentration))

			for i, conc := range noxConcentration {
				rr := aqhealth.RRpm25Linear(conc)
				healthImpactsso2[i] = aqhealth.Deaths(rr, pop[i], c.mr[i])
				whitehealthImpactsso2[i] = aqhealth.Deaths(rr, whitepop[i], c.mr[i])
				blackhealthImpactsso2[i] = aqhealth.Deaths(rr, blackpop[i], c.mr[i])
				nativehealthImpactsso2[i] = aqhealth.Deaths(rr, nativepop[i], c.mr[i])
				asianealthImpactsso2[i] = aqhealth.Deaths(rr, asianpop[i], c.mr[i])
				otherhealthImpactsso2[i] = aqhealth.Deaths(rr, otherpop[i], c.mr[i])
				latinohealthImpactsso2[i] = aqhealth.Deaths(rr, latinopop[i], c.mr[i])
				povertyealthImpactsso2[i] = aqhealth.Deaths(rr, povertyop[i], c.mr[i])
				twoxpovalthImpactsso2[i] = aqhealth.Deaths(rr, twoxpovpop[i], c.mr[i])
			}

			//VOC health impacts
			healthImpactsvoc := make([]float64, len(vocConcentration))
			whitehealthImpactsvoc := make([]float64, len(vocConcentration))
			blackhealthImpactsvoc := make([]float64, len(vocConcentration))
			nativehealthImpactsvoc := make([]float64, len(vocConcentration))
			asianealthImpactsvoc := make([]float64, len(vocConcentration))
			otherhealthImpactsvoc := make([]float64, len(vocConcentration))
			latinohealthImpactsvoc := make([]float64, len(vocConcentration))
			povertyealthImpactsvoc := make([]float64, len(vocConcentration))
			twoxpovalthImpactsvoc := make([]float64, len(vocConcentration))

			for i, conc := range noxConcentration {
				rr := aqhealth.RRpm25Linear(conc)
				healthImpactsvoc[i] = aqhealth.Deaths(rr, pop[i], c.mr[i])
				whitehealthImpactsvoc[i] = aqhealth.Deaths(rr, whitepop[i], c.mr[i])
				blackhealthImpactsvoc[i] = aqhealth.Deaths(rr, blackpop[i], c.mr[i])
				nativehealthImpactsvoc[i] = aqhealth.Deaths(rr, nativepop[i], c.mr[i])
				asianealthImpactsvoc[i] = aqhealth.Deaths(rr, asianpop[i], c.mr[i])
				otherhealthImpactsvoc[i] = aqhealth.Deaths(rr, otherpop[i], c.mr[i])
				latinohealthImpactsvoc[i] = aqhealth.Deaths(rr, latinopop[i], c.mr[i])
				povertyealthImpactsvoc[i] = aqhealth.Deaths(rr, povertyop[i], c.mr[i])
				twoxpovalthImpactsvoc[i] = aqhealth.Deaths(rr, twoxpovpop[i], c.mr[i])
			}

			totals := rec.Totals()
			line := make([]string, len(row0))
			line[0] = scc
			for i, pol := range poltitle {
				for pollutant, value := range totals {
					if pollutant.Name == pol {
						line[i+1] = fmt.Sprintf("%g", value.Value())
					}
				}
			}
			line = line[:len(line)-51]

			line = append(line, sccDesc[scc])

			line = append(line, fmt.Sprintf("%g", floats.Sum(pm25Concentration)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(noxConcentration)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(PRIConcentration)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(so2Concentration)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(vocConcentration)))

			line = append(line, fmt.Sprintf("%g", floats.Sum(healthImpactspm)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(whitehealthImpactspm)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(blackhealthImpactspm)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(nativehealthImpactspm)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(asianealthImpactspm)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(otherhealthImpactspm)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(latinohealthImpactspm)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(povertyealthImpactspm)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(twoxpovalthImpactspm)))

			line = append(line, fmt.Sprintf("%g", floats.Sum(healthImpactsnox)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(whitehealthImpactsnox)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(blackhealthImpactsnox)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(nativehealthImpactsnox)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(asianealthImpactsnox)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(otherhealthImpactsnox)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(latinohealthImpactsnox)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(povertyealthImpactsnox)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(twoxpovalthImpactsnox)))

			line = append(line, fmt.Sprintf("%g", floats.Sum(healthImpactspri)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(whitehealthImpactspri)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(blackhealthImpactspri)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(nativehealthImpactspri)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(asianealthImpactspri)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(otherhealthImpactspri)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(latinohealthImpactspri)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(povertyealthImpactspri)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(twoxpovalthImpactspri)))

			line = append(line, fmt.Sprintf("%g", floats.Sum(healthImpactsso2)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(whitehealthImpactsso2)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(blackhealthImpactsso2)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(nativehealthImpactsso2)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(asianealthImpactsso2)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(otherhealthImpactsso2)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(latinohealthImpactsso2)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(povertyealthImpactsso2)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(twoxpovalthImpactsso2)))

			line = append(line, fmt.Sprintf("%g", floats.Sum(healthImpactsvoc)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(whitehealthImpactsvoc)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(blackhealthImpactsvoc)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(nativehealthImpactsvoc)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(asianealthImpactsvoc)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(otherhealthImpactsvoc)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(latinohealthImpactsvoc)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(povertyealthImpactsvoc)))
			line = append(line, fmt.Sprintf("%g", floats.Sum(twoxpovalthImpactsvoc)))

			err5 := ww.Write(line)
			if err5 != nil {
				panic(err5)
			}
		}
		ww.Flush()
		w.Close()

	}
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

func pmtotal(rec aep.Record, c *Config) [][]float64 {
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
	var pm25Concentration [][]float64
	row1 := make([]float64, len(c.gridCells))
	row2 := make([]float64, len(c.gridCells))
	row3 := make([]float64, len(c.gridCells))
	row4 := make([]float64, len(c.gridCells))
	row5 := make([]float64, len(c.gridCells))

	for pollutant, value := range totals {

		inmapPol, ok := pollutantCrosswalk[pollutant.Name]
		if !ok {
			panic(fmt.Errorf("missing pollutant %s", pollutant.Name))
		}
		floats.AddScaled(row1, value.Value(), concentrationSurrogates[inmapPol])
	}

	for _, value := range totals {

		inmapPol, ok := pollutantCrosswalk["NOX"]
		if !ok {
			panic(fmt.Errorf("missing pollutant NOX"))
		}
		floats.AddScaled(row2, value.Value(), concentrationSurrogates[inmapPol])

	}

	for _, value := range totals {

		inmapPol, ok := pollutantCrosswalk["PM25-PRI"]
		if !ok {
			panic(fmt.Errorf("missing pollutant PM25-PRI"))
		}
		floats.AddScaled(row3, value.Value(), concentrationSurrogates[inmapPol])

	}

	for _, value := range totals {

		inmapPol, ok := pollutantCrosswalk["SO2"]
		if !ok {
			panic(fmt.Errorf("missing pollutant SO2"))
		}
		floats.AddScaled(row4, value.Value(), concentrationSurrogates[inmapPol])

	}

	for _, value := range totals {

		inmapPol, ok := pollutantCrosswalk["VOC"]
		if !ok {
			panic(fmt.Errorf("missing pollutant VOC"))
		}
		floats.AddScaled(row5, value.Value(), concentrationSurrogates[inmapPol])

	}
	//	pm25Concentration = append(pm25Concentration, row1, row2)
	//pm25Concentration = append(pm25Concentration, row2)
	pm25Concentration = append(pm25Concentration, row1, row2, row3, row4, row5)
	// fmt.Println("row 1", floats.Sum(row1))
	// fmt.Println("row2", floats.Sum(row2))
	//
	// for i, pm := range pm25Concentration {
	// 	fmt.Println("pmdata", i, floats.Sum(pm))
	// }

	return pm25Concentration
}
