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
    + `007_createFiguresDiagramsTables`: scripts to create figures, diagrams, and some tables for paper
    + `008_exploratoryAnalysis`: scripts to run exploratory analyses conducted during peer review
    + `primaryScript.R`: run this script FIRST so that the scripts in this directory populate the appropriate directories.

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
    + `008_sampleSizeSummaries`: holds text files summarizing sample sizes of datasets
    + `009_exploratoryAnalysis`: holds text files related to the exploratory analyses conducted
    + `figures`: directory to hold figures produced by scripts   

---

### Scripts

#### `001_prepareDataForModels`: scripts to prepare data for model fitting

Each file is an R script to prepare observations from field surveys or experiments to be fit by the statistical models. Preparing observations for model fitting involves organizing raw data so that it can be read and used by the model scripts, which are written in the JAGS language.

#### `002_fitStatisticalModels`: scripts to fit statistical models

Each file is an R script to fit the statistical models (written in JAGS; found in the folder `models/`) to the data (output from the files in `001_prepareDataForModels`). The table below describes the relationship between the raw data, data prep scripts, JAGS models, and model fitting scripts.

| Raw data            | Script to prepare raw data for model fitting | JAGS model  | Script to fit JAGS model to data  |
| --------------------------------- | ------------- | ------------- | ------------- |
| `seedBagsData.csv`; `seedlingFruitingPlantCountsPermanentPlots.csv`; `countFruitsPerPlantFromPermanentPlots.csv`; `countSeedPerFruit.csv`; `countFruitsPerPlantAllPlants.csv`                        | `01_prepDataForSeedModel.R`              | `jags-seedBagExperiment.R`              | `01_modelScriptsSeedBagExperiment.R` |
| `viabilityData.csv`    | `02_prepDataForViabilityModel.R`              | `jags-viabilityTrials.R`              | `02_modelScriptsViabilityTrials.R` |
| `seedlingFruitingPlantCountsPermanentPlots.csv`    | `03_prepDataForSeedlingSurvivalModel.R`              | `jags-seedlingSurvival.R`              | `03_modelScriptsSeedlingSurvival.R` |
| `countFruitsPerPlantAllPlants.csv`; `countUndamagedDamagedFruitsPerPlantAllPlants.csv`    | `04_prepDataForFruitsPerPlantModel.R`              | `jags-fruitsPerPlant.R`              | `04_modelScriptsFruitsPerPlant.R` |
| `countSeedPerFruit.csv`    | `05_prepDataForSeedsPerFruitModel.R`              | `jags-seedsPerFruit.R`              | `05_modelScriptsSeedsPerFruit.R` |

#### `003_runStatisticalModelDiagnostics`: scripts to run model diagnostics

The folder includes a file with a utility function, `00_tracePlotFunction.R`, to produce trace plots that is used in all other files in the folder. All other files run diagnostic checks on the MCMC samples and produce (1) trace plots and (2) calculate R-hat values and the Heidelberg-Welch diagnostic.

#### `004_checkStatisticalModels`: scripts to perform model checks

The folder includes files that conduct model checks for each model. Model checks vary depending on the nature of the data, and include graphical checks as well as posterior predictive checks.

#### `005_calculatePopulationModelParameters`: scripts to calculate parameters for population model

The folder includes files that use the posterior distribution of parameter estimates to calculate parameters for the population model. The script `00_extractParameters.R` extracts the posterior distributions for parameters used in subsequent scripts. All other scripts in the folder use these posteriors to calculate vital rates. These calculations are described in the section <b>Methods: Computing vital rates</b> in the main text of the manuscript, and in <b>Appendix S4: Computing vital rates</b> in the supplementary materials of the manuscript.

#### `006_testHypotheses`: scripts to test hypotheses in the manuscript

The folder contains scripts to test the hypotheses in the manuscript, and produce results figures in the main text and supplement.

- `00_utilityFunctions.R`: Functions used in other scripts in this folder.
- `01_demographicBetHedgingTest.R`: Script to conduct the demographic test of bet hedging, described in <b>Methods: Analysis: Demographic test of bet hedging</b> in the main text. The script produces panels A-C in Figure 3. Also produces the zoomed in version of Figure 3B included in the appendix.
- `01.1_demographicBetHedgingTest.R`: Script to conduct the demographic test of bet hedging, that uses 'quasi-complete germination' (prob. germination = 0.99) described in <b>Methods: Analysis: Demographic test of bet hedging</b> in the main text. This analysis was suggested by a reviewer during peer review. Results of the analysis are given in the appendix.
- `02_optimalObservedGerminationDensityIndependentModelGridSearch.R`: Script to calculate the optimal germination fractions using a density-independent model for bet hedging. The script uses a grid search to find the optimal germination fractions. NOTE: the output from this script is used to verify the results from the 1D optimization method in (`03_optimalObservedGerminationDensityIndependentModelOptimization.R`) and IS NOT USED to create figures in the main text.
- `03_optimalObservedGerminationDensityIndependentModelOptimization.R`: Script to calculate the optimal germination fractions using a density-independent model for bet hedging, as described in <b>Methods: Analysis: Density-independent model for germination fractions</b> in the main text. The script uses a 1-dimensional optimization routine to find the optimal germination fraction. NOTE: the output from this script IS USED to create Figure 4.
- `04_correlationGerminationSeedSurvival.R`: Script to calculate the correlation between germination fractions and seed survival, described in <b>Methods: Analysis: Correlation between germination and seed survival</b> in the main text. The script produces panel A in Figure 5.
- `05_correlationGerminationReproductiveSuccess.R`: Script to calculate the correlation between germination fractions and variability in per-capita reproductive success, described in <b>Methods: Analysis: Correlation between germination and variability in per-capita reproductive success</b> in the main text. The script produces panel B in Figure 5.
- `06_demographicBetHedgingTestUncertainty.R`: Script to calculate uncertainty in the demographic test of bet hedging, described in <b>S5.1 Accounting for parameter uncertainty in demographic test of bet hedging</b> in the supplementary materials. The script produces panel A-C in Figure S7.
- `07_optimalObservedGerminationDensityIndependentModelUncertainty.R`: Script to calculate uncertainty in the estimates for optimal germination fractions, described in <b>S5.2 Accounting for parameter uncertainty in optimal germination fractions</b> in the supplementary materials. The script produces Figure S8.
- `08_optimalGerminationSensitivity.R`: Script to calculate sensitivity of optimal germination fractions to seed mortality, described in <b>S5.3 Seed mortality before and after germination have opposing effects on
optimal germination</b> in the supplementary materials. The script produces Figures S9 and S10.

#### `007_createFiguresDiagramsTables`: scripts to create additional figures and diagrams for paper

- `00_utilityFunctions.R`: Modifies `ggsn::scalebar` for the map of populations.
- `01_populationMap.R`: Script to produce map of study populations. Produces Figure 1A in the main text.
- `02_reproductiveFailureMissingDataDiagram.R`: Script to create visual summary of years in which populations have observations for per-capita reproductive success, in which no seedlings survived in permanent plots, and in which there is missing data. Produces Figure 1B in the main text.
- `03_abovegroundModelDiagram.R`: Script to produce parts of Figure 2A. Specifically, this script produces the graphs that describe models for aboveground components of demography.
- `04_belowgroundModelDiagram.R`: Script to produce parts of Figure 2B. Specifically, this script produces the graphs that describe models for belowground components of demography.
- `05_summarizeSampleSizes.R`: Script to summarize sample sizes of datasets used in the study. Used to produce Tables S2-S9.
- `06_plotPopulationModelParameters.R`: Script to plot population model parameters (i.e., graphical summary of vital rate parameters).

#### `008_exploratoryAnalysis`: scripts to run exploratory analyses conducted during peer review

- `01_precipitationReproductiveSuccessAnalysis.R`: Script to analyze the relationship between spring precipitation and per-capita reproductive success.
- `02_seedSeedlingDensityDependenceAnalysis.R`: Script to assess evidence for density-dependence in the transition from seed to seedling.
- `03_seedlingFruitingPlantDensityDependenceAnalysis.R`: Script to assess evidence for density-dependence in the transition from seedling to fruiting plant.
