# clarkia-bet-hedging

### Bet hedging is not sufficient to explain intraspecific variation in germination patterns of a winter annual plant

### Authors

  - Gregor-Fausto Siegmund, Cornell University, <gs589@cornell.edu>
  - Vincent M. Eckhart, Grinnell College
  - David A. Moeller, University of Minnesota
  - Monica A. Geber, Cornell University

  This repository contains the scripts for the project on bet hedging in <i>Clarkia xantiana</i>.

### Contributions

GS and MAG conceived of the ideas and analysis, using data collected by MAG, VME, and DAM. GS wrote the scripts, analyzed the data, and wrote the manuscript with input from MAG. All authors contributed critically to drafts of the manuscript.

-----

### Abstract

[include abstract here]

-----

### Repository Directory

The repository is organized so that the analysis in the paper can be replicated.

- `data`: Contains data used in the study. Contents of data files are documented further down in the README.
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
- `outputs`: Folders to save output from R scripts.
- `scripts`: Contains R scripts to process data, fit models, and analyze output.
    + `primaryScript.R`: script
    + `001_processData`: scripts to prepare data for model fitting
    + `002_fitStatisticalModels`: scripts to fit statistical models
    + `003_runStatisticalModelDiagnostics`: scripts to run model diagnostics
    + `004_checkStatisticalModels`: scripts to perform model checks
    + `005_calculatePopulationModelParameters`: scripts to calculate parameters for structured population model from statistical model output
    + `006_testHypotheses`: scripts to test hypotheses in the manuscript
    + `007_createFiguresDiagrams`: scripts to create figures for paper
    + `008_createTexDiagrams`: scripts to create diagrams for paper with LaTeX

Running `primaryScript.R` in the appropriate directory will create the folders `outputs` and `products` with the following file structure. Note that replicating the simulation and model fitting may be slow. We recommend testing the code in `003_statisticalModelFitting` on a smaller number of replicates than the default.

- `outputs`: Folder for output of simulations and model fitting
    + `001_simulateObservations`: simulated observations from seed bag burial  and seed addition experiments
        * `01_identifiability`: simulated observations for showing effect of identifiability on parameter estimation.         
        * `02_estimability`: simulated observations to compare statistical properties of estimates from seed bag burial vs. seed addition experiments
        * `03_misspecification`: simulated observations to evaluate consequences of incorrect assumptions about seed mortality
    + `002_statisticalModelFitting`: scripts to fit JAGS models to simulated observations
        * `01_identifiability`: samples from posterior distribution of model fits         
        * `02_estimability`: samples from posterior distribution of model fits
        * `03_misspecification`: samples from posterior distribution of model fits
- `products`: Folder for figures
    + `figures`: directory to hold figures produced by scripts     
    + `tables`: directory to hold tables produced by scripts    

---

### Data

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

Observations of total fruit equivalents per plant from field surveys in 2006-2012. Data come from up to 15 plants per permanent plot. The counts are of “undamaged fruit equivalents” per plant, and were made by counting the number of undamaged fruits, and then counting the damaged fruits and estimating how many undamaged fruits these damaged fruits amounted to.

| variableName          | description |
| ------------          | ----------- |
| site                  | Site acronym |
| transect              | Transect number |
| position              | Plot number |
| plantNumber           | Variable indexing the plants on which fruits per plant was recorded; the index starts at 1 in each permanent plot in each year at each site         |
| countFruitsPerPlant   | Count of fruits on a plant; the count combined undamaged and damaged fruits into a single value when fruits were surveyed in the field   |
| year                  | Year in which observations were made     |

#### countSeedPerFruit.csv

Observations of seeds per fruit from seeds collected in the field in 2006-2021. From 2006-2021, data are for counts of seeds per undamaged fruit for up to roughly 30 fruits per site. From 2013-2021, the data also include the counts of seeds per damaged fruit.

| variableName                      | description   |
| --------------------------------- | ------------- |
| site                              | Site acronym  |
| year                              | Year in which observations were made  |
| damaged                           | Binary variable indicating whether seeds were counted in a fruit that was undamaged (0=false) or damaged (1=true)                            |
| sdno                              | Count of seeds    |

#### countUndamagedDamagedFruitsPerPlantAllPlants.csv

Observations of undamaged and damaged fruits per plant from field surveys in 2013-2021. Data come from up to 15 plants per permanent plot, plus plants found across the site. The counts are of “undamaged fruits” and "damaged fruits" on each plant.

| variableName                      | description  |
| --------------------------------- | ------------ |
| site                              | Site acronym  |
| countUndamagedFruitNumberPerPlant | Count of undamaged fruits on a plant  |
| countDamagedFruitNumberPerPlant   | Count of damaged fruits on a plant  |
| permanentPlot                     | Binary variable indicating whether the plant on which fruits were counted was growing in one of 30 permanent plots at the site (1=true, 0=false) |
| year                              | Year in which observations were made  |

#### countUndamagedDamagedFruitsPerPlantFromPermanentPlots.csv

Observations of undamaged and damaged fruits per plant from field surveys in 2013-2021. Data come from up to 15 plants per permanent plot. The counts are of “undamaged fruits” and "damaged fruits" on each plant.

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

Observations of seedlings and fruiting plants in permanent plots from 2006-2021. Data are counts of seedlings in permanent plots in January/February censuses, and counts of fruiting plants in permanent plots in June/Jly censuses.

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
