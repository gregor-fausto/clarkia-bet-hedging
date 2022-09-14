# clarkia-bet-hedging

### Bet hedging is not sufficient to explain intraspecific variation in germination patterns of a winter annual plant

### Authors

  - Gregor-Fausto Siegmund, Cornell University, <gs589@cornell.edu>
  - Vincent M. Eckhart, Grinnell College
  - David A. Moeller, University of Minnesota
  - Monica A. Geber, Cornell University

  This repository contains the scripts for the project on bet hedging in <i>Clarkia xantiana</i> ssp. <i>xantiana</i>.

-----

### Abstract

Bet hedging consists of life history strategies that buffer against environmental variability by trading off immediate and long-term fitness. Delayed germination in annual plants is a classic example of bet hedging, and is often invoked to explain low germination fractions. We examined whether bet hedging explains low and variable germination fractions among 20 populations of the winter annual plant <i>Clarkia xantiana</i> ssp. <i>xantiana</i> that experience substantial variation in reproductive success among years. Leveraging 15 years of demographic monitoring and 3 years of field germination experiments, we assessed the fitness consequences of seed banks and compared optimal germination fractions from a density-independent bet-hedging model to observed germination fractions. We did not find consistent evidence of bet hedging or the expected trade-off between arithmetic and geometric mean fitness, though delayed germination increased long-term fitness in 7 of 20 populations. Optimal germination fractions were 2 to 5 times higher than observed germination fractions, and among-population variation in germination fractions were not correlated with risks across the life cycle. Our comprehensive test suggests that bet hedging is insufficient to explain the observed germination patterns. Understanding variation in germination strategies will likely require integrating bet hedging with complementary forces shaping the evolution of delayed germination.

-----

### Contributions

GS and MAG conceived of the ideas and analysis, using data collected by MAG, VME, and DAM. GS wrote the scripts, analyzed the data, and wrote the manuscript with input from MAG. All authors contributed critically to drafts of the manuscript.

-----

### Repository Directory

The repository is organized so that the analyses in the paper can be replicated.

- `data`: Contains data used in the study. The data files listed below are found in this Dryad repository: [link]. Contents of data files are documented further down in the README, as well as a README in the data repository.
    + `metadata`
      * `attributes.csv`: Describes the variables for all data files.
      * `creators.csv`: Documents the creators of the data files.
    + `countFruitsPerPlantAllPlants.csv`
    + `countFruitsPerPlantFromPermanentPlots.csv`
    + `countSeedPerFruit.csv`
    + `countUndamagedDamagedFruitsPerPlantAllPlants.csv`
    + `countUndamagedDamagedFruitsPerPlantFromPermanentPlots.csv`
    + `seedBagsData.csv`
    + `seedlingFruitingPlantCountsPermanentPlots.csv`
    + `siteAbioticData.csv`
    + `viabilityData.csv`
    + `README.md`: README file for the data folder.
- `models`: Contains the statistical models written in JAGS language.
- `outputs`: Folder to save output from R scripts.
- `products`: Folder to save figures and diagrams for manuscript.
- `scripts`: Contains R scripts to process data, fit models, and analyze output.
    + `001_prepareDataForModels`: scripts to prepare data for model fitting
    + `002_fitStatisticalModels`: scripts to fit statistical models
    + `003_runStatisticalModelDiagnostics`: scripts to run model diagnostics
    + `004_checkStatisticalModels`: scripts to perform model checks
    + `005_calculatePopulationModelParameters`: scripts to calculate parameters for population model
    + `006_testHypotheses`: scripts to test hypotheses in the manuscript
    + `007_createFiguresDiagramsTables`: scripts to create figures, diagrams, and some tables for paper

Running `primaryScript.R` in the appropriate directory will create the folders `outputs` and `products` with the following file structure. Note that replicating the simulation and model fitting may be slow. We recommend testing the code in `003_statisticalModelFitting` on a smaller number of replicates than the default.

- `outputs`: Folder for output of simulations and model fitting
    + `001_prepareDataForModels`: holds data in format ready for fitting models with JAGS
    + `002_fitStatisticalModels`: holds (1) data used to fit models and (2) posterior samples obtained via MCMC
    + `003_runStatisticalModelDiagnostics`: holds trace plots and histograms of R-hat
    + `004_checkStatisticalModels`: holds PDFs of model checks
    + `005_calculatePopulationModelParameters`: holds (1) posterior samples of the statistical model parameters and derived quantities including (2) parameters of the population model in lists, (3) parameters of the population model in matrices, and (4) computed annual values for per-capita reproductive success
    + `006_hypothesisTesting`: holds files with optimal germination fractions and results of demographic bet hedging test
    + `007_createFiguresDiagrams`: holds shapefiles for creating the map in Figure 1
- `products`: Folder for figures
    + `figures`: directory to hold figures produced by scripts   

---

### Data

The data associated with the scripts in this Github folder are archived in the following Dryad repository: [link]. Once the data files from the Dryad repository are downloaded, they can be added to the `data` folder. The contents of the data files are documented here.

Briefly, we used field surveys and experiments to observe components of above- and below-ground demography for 20 populations across the range of <i>Clarkia xantiana</i> ssp. <i>xantiana</i>. To collect data on seedling survival, fruit production, and seed set, we used field surveys. In each population, these surveys included observations in both permanent plots as well as additional, haphazardly sampled plots arrayed across the population. To observe emergence of seedlings and seeds remaining intact in the soil seed bank, we conducted field experiments, which were complemented with lab experiments in order to assay viability of seeds. Brief details for each dataset are provided here, and the manuscript associated with these datasets describes the survey and experimental methods in further detail.

#### countFruitsPerPlantAllPlants.csv

Observations of total fruit equivalents per plant from field surveys in 2006-2012. Data come from up to 15 plants per permanent plot, plus plants found across the site. The counts are of “undamaged fruit equivalents” per plant, and were made by counting the number of undamaged fruits, and then counting the damaged fruits and estimating how many undamaged fruits these damaged fruits amounted to.

| variableName                      | description   |
| --------------------------------- | ------------- |
| site                              | Site acronym              |
| countFruitNumberPerPlant          | Count of fruits on a plant; the count combined undamaged and damaged fruits into a single value when fruits were surveyed in the field  |
| permanentPlot                     | Binary variable indicating whether the plant on which fruits were counted was growing in one of 30 permanent plots at the site (1=true) or located in either an additional, haphazardly located plot or found by surveying the site (0=false) |
| damage                            | Binary variable indicating whether the fruit that was used to count seeds was undamaged (0=false) or damaged (1=true); all observations in this dataset have NA here because the count is of total fruit equivalents                               |
| year                              | Year in which observations were made

#### countFruitsPerPlantFromPermanentPlots.csv

Observations of total fruit equivalents per plant from field surveys in 2007-2012. Data come from up to 15 plants per permanent plot. The counts are of “undamaged fruit equivalents” per plant, and were made by counting the number of undamaged fruits, and then counting the damaged fruits and estimating how many undamaged fruits these damaged fruits amounted to.

| variableName          | description |
| ------------          | ----------- |
| site                  | Site acronym |
| transect              | Transect number |
| position              | Plot number |
| plantNumber           | Variable indexing the plants on which fruits per plant was recorded; the index starts at 1 in each permanent plot in each year at each site         |
| countFruitsPerPlant   | Count of fruits on a plant; the count combined undamaged and damaged fruits into a single value when fruits were surveyed in the field   |
| year                  | Year in which observations were made     |

#### countSeedPerFruit.csv

Observations of seeds per fruit from seeds collected in the field in 2006-2020. From 2006-2020, data are for counts of seeds per undamaged fruit for up to roughly 30 fruits per site. From 2013-2020, the data also include the counts of seeds per damaged fruit.

| variableName                      | description   |
| --------------------------------- | ------------- |
| site                              | Site acronym  |
| year                              | Year in which observations were made  |
| damaged                           | Binary variable indicating whether seeds were counted in a fruit that was undamaged (0=false) or damaged (1=true)                            |
| sdno                              | Count of seeds    |

#### countUndamagedDamagedFruitsPerPlantAllPlants.csv

Observations of undamaged and damaged fruits per plant from field surveys in 2013-2020. Data come from up to 15 plants per permanent plot, plus plants found across the site. The counts are of “undamaged fruits” and "damaged fruits" on each plant.

| variableName                      | description  |
| --------------------------------- | ------------ |
| site                              | Site acronym  |
| countUndamagedFruitNumberPerPlant | Count of undamaged fruits on a plant  |
| countDamagedFruitNumberPerPlant   | Count of damaged fruits on a plant  |
| permanentPlot                     | Binary variable indicating whether the plant on which fruits were counted was growing in one of 30 permanent plots at the site (1=true, 0=false) |
| year                              | Year in which observations were made  |

#### countUndamagedDamagedFruitsPerPlantFromPermanentPlots.csv

Observations of undamaged and damaged fruits per plant from field surveys in 2013-2020. Data come from up to 15 plants per permanent plot. The counts are of “undamaged fruits” and "damaged fruits" on each plant.

| variableName                 | description  |
| ---------------------------- | ------------ |
| site                         | Site acronym  |
| transect                     | Transect number    |
| position                     | Plot number     |
| plantNumber                  | Variable indexing the plants on which fruits per plant was recorded; the index starts at 1 in each permanent plot in each year at each site  |
| countUndamagedFruitsPerPlant | Count of undamaged fruits on a plant      |
| year                         | Year in which observations were made |
| countDamagedFruitsPerPlant   | Count of damaged fruits on a plant  |

#### seedBagsData.csv

Observations of intact seeds and germinants from the seed bag burial experiment, conducted from October 2005-October 2008. The experiment consisted of 3 experimental rounds, starting in October 2005, 2006, and 2007. The data in this file includes counts of intact seeds and germinants. The data for the lab component of the seed bag burial experiment is found in viabilityData.csv.

| variableName   | description                                                                                               |
| -------------- | --------------------------------------------------------------------------------------------------------- |
| site           | Site acronym                                                                                              |
| transect       | Transect number                                                                                           |
| position       | Plot number                                                                                               |
| bagNo          | Bag number                                                                                                |
| round          | Experimental round in which the bag was buried                                                            |
| yearStart      | Year in which the experimental round started; bags were buried In October of this year                    |
| age            | Number of years the bag had been buried when it was collected from the field in October                   |
| yearData       | Year in which observations were made for this bag; the bag was dug up in January and October of this year |
| seedlingJan    | Count of seedlings in the bag when it was dug up in January                                               |
| intactJan      | Count of intact, ungerminated seeds in the bag when it was dug up in January                              |
| totalJan       | Sum of the count of seedlings and intact, ungerminated seeds in the bag when it was dug up in January     |
| intactOct      | Count of intact seeds in the bag when it was dug up in October                                            |

#### seedlingFruitingPlantCountsPermanentPlots.csv

Observations of seedlings and fruiting plants in permanent plots from 2006-2020. Data are counts of seedlings in permanent plots in January/February censuses, and counts of fruiting plants in permanent plots in June/Jly censuses.

| variableName   | description                                                                             |
| -------------- | --------------------------------------------------------------------------------------- |
| site           | Site acronym                                                                            |
| transect       | Transect number                                                                         |
| position       | Plot number                                                                             |
| year           | Year in which observations were made                                                    |
| seedlingNumber | Count of seedlings in 0.5 m x 1 m permanent plot; NA if not recorded                    |
| fruitplNumber  | Count of fruiting plants in 0.5 m x 1 m permanent plot; NA if not recorded              |

#### siteAbioticData.csv

Summary of abiotic variables associated with each study site in the long-term study of Clarkia xantiana ssp. xantiana demography.  

| variableName | description                                                                             |
| ------------ | --------------------------------------------------------------------------------------- |
| site         | Site acronym                                                                            |
| siteName     | Full name for site                                                                              |
| easting      | Geographic position reported as eastward measured distance, reported for UTM zone 11, NAD 1927 (units of 100 meters)                                         |
| northing          | Geographic position reported as northward measured distance, reported for UTM zone 11, NAD 1927 (units of 100 meters) |
| elevation        | Height above sea level (meters)                           |
| area    | Size of the Clarkia xantiana ssp. xantiana population at the study site (hectares)                                          |
| surfaceRock    | Dominant soil parent material                    |


#### viabilityData.csv

Observations from lab germination assays and viability assays that are associated with the seed bag burial experiment. The lab experiments were conducted with seeds from seed bags that were recovered from the field in October 2006, 2007, and 2008. The data for the field component of the seed bag burial experiment is found in seedBagsData.csv.

| variableName | description                                                                             |
| ------------ | --------------------------------------------------------------------------------------- |
| site         | Site acronym                                                                            |
| bagNo        | Bag number                                                                              |
| round        | Experimental round in which the bag was buried                                          |
| age          | Number of years the bag had been buried when it was collected from the field in October |
| block        | Block number for the Petri dishes used in germination trials                            |
| germStart    | Count of seeds starting the germination trials                                          |
| germCount    | Count of seeds germinating over the course of the germination trials                    |
| viabStart    | Count of seeds tested in the viability trials                                           |
| viabStain    | Count of seeds staining red in the viability trials                                     |
