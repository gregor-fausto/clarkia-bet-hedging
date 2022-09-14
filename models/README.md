# clarkia-bet-hedging
## README for models folder

### Repository Directory

- `models`: Contains the statistical models written in JAGS language.

---

### Models

Each file is an R script containing a statistical models written in JAGS to data (output from the files in `../scripts/001_prepareDataForModels`). The table below describes the relationship between the raw data, data prep scripts, JAGS models, and model fitting scripts. The `models` folder contains the JAGS model scripts.

| Raw data            | Script to prepare raw data for model fitting | JAGS model  | Script to fit JAGS model to data  |
| --------------------------------- | ------------- | ------------- | ------------- |
| `seedBagsData.csv`; `seedlingFruitingPlantCountsPermanentPlots.csv`; `countFruitsPerPlantFromPermanentPlots.csv`; `countSeedPerFruit.csv`; `countFruitsPerPlantAllPlants.csv`                        | `01_prepDataForSeedModel.R`              | `jags-seedBagExperiment.R`              | `01_modelScriptsSeedBagExperiment.R` |
| `viabilityData.csv`    | `02_prepDataForViabilityModel.R`              | `jags-viabilityTrials.R`              | `02_modelScriptsViabilityTrials.R` |
| `seedlingFruitingPlantCountsPermanentPlots.csv`    | `03_prepDataForSeedlingSurvivalModel.R`              | `jags-seedlingSurvival.R`              | `03_modelScriptsSeedlingSurvival.R` |
| `countFruitsPerPlantAllPlants.csv`; `countUndamagedDamagedFruitsPerPlantAllPlants.csv`    | `04_prepDataForFruitsPerPlantModel.R`              | `jags-fruitsPerPlant.R`              | `04_modelScriptsFruitsPerPlant.R` |
| `countSeedPerFruit.csv`    | `05_prepDataForSeedsPerFruitModel.R`              | `jags-seedsPerFruit.R`              | `05_modelScriptsSeedsPerFruit.R` |
