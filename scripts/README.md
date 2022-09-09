# clarkia-bet-hedging
## README for scripts folder

### Repository Directory

- `scripts`: Contains R scripts to process data, fit models, and analyze output.
    + `001_prepareDataForModels`: scripts to prepare data for model fitting
    + `002_fitStatisticalModels`: scripts to fit statistical models
    + `003_runStatisticalModelDiagnostics`: scripts to run model diagnostics
    + `004_checkStatisticalModels`: scripts to perform model checks
    + `005_calculatePopulationModelParameters`: scripts to calculate parameters for population model
    + `006_testHypotheses`: scripts to test hypotheses in the manuscript
    + `007_createFiguresDiagrams`: scripts to create diagrams for paper

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
    + `figuresForManuscript`: directory to hold figures composed with LaTeX   
    + `tables`: directory to hold tables produced by scripts    
    + `texDiagrams`: directory to TeX diagrams in the manuscript

---

### Scripts

#### `001_prepareDataForModels`: scripts to prepare data for model fitting

Each file is an R script to prepare observations from field surveys or experiments to be fit by the statistical models. Preparing observations for model fitting involves organizing raw data so that it can be read and used by the model scripts, which are written in the JAGS language.

#### `002_fitStatisticalModels`: scripts to fit statistical models

Each file is an R script to fit the statistical models (written in JAGS; found in the folder `models/`) to the data (output from the files in `001_prepareDataForModels`). The table below describes the relationship between the raw data, data prep scripts, JAGS models, and model fitting scripts.

| Raw data            | Script to prepare raw data for model fitting | JAGS model  | Script to fit JAGS model to data  |
| --------------------------------- | ------------- | ------------- | ------------- |
| `seedBagsData.csv`; `seedlingFruitingPlantCountsPermanentPlots.csv`; `countFruitsPerPlantFromPermanentPlots.csv`; `countSeedPerFruit.csv`; `countFruitsPerPlantAllPlants.csv`                        | `01_prepDataForSeedModel.R`              | `jags-seedBagExperiment.R`              | `01_modelScriptsSeedBagExperiment.R`
