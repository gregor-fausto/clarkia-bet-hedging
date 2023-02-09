####
####
# Primary script to create output directories
####
####

# - Create output directories ----

dir.create("../outputs")

# - +Create directory to hold data in format ready for fitting models with JAGS

dir.create("../outputs/001_prepareDataForModels")

# - +Create directory to hold output of fitting models

dir.create("../outputs/002_fitStatisticalModels")
tmp <- "../outputs/002_statisticalModelFitting"

# - ++Create directories to hold (1) data used to fit models and (2) posterior samples obtained via MCMC

dir.create(paste0(tmp,"/data"))
dir.create(paste0(tmp,"/mcmcSamples"))

# - +Create directory to hold trace plots and histograms of R-hat

dir.create("../outputs/003_runStatisticalModelDiagnostics")
tmp <- "../outputs/003_runStatisticalModelDiagnostics"

# - ++Create directories to hold (1) data used to fit models and (2) posterior samples obtained via MCMC

dir.create(paste0(tmp,"/01_tracePlots"))
dir.create(paste0(tmp,"/02_rhatDistribution"))

# - +Create directory to hold PDFs of model checks

dir.create("../outputs/004_checkStatisticalModels")
tmp <- "../outputs/004_checkStatisticalModels"

# - ++Create directories to hold model checks to include in revised supplement

dir.create(paste0(tmp,"/01_modelChecksSupplement"))

# - +Create directory to hold population model parameters

dir.create("../outputs/005_calculatePopulationModelParameters")
tmp <- "../outputs/005_calculatePopulationModelParameters"

# - ++Create directories to hold (1) posterior samples of the statistical model parameters and derived quantities including (2) parameters of the population model in lists, (3) parameters of the population model in matrices, and (4) computed annual values for per-capita reproductive success

dir.create(paste0(tmp,"/01_parameterPosteriorDistributions"))
dir.create(paste0(tmp,"/02_populationModelParameters"))
dir.create(paste0(tmp,"/03_populationModelParametersMatrix"))
dir.create(paste0(tmp,"/04_reproductiveSuccess"))

# - +Create directory to hold files with optimal germination fractions and results of demographic bet hedging test

dir.create("../outputs/006_hypothesisTesting")

# - +Create directory to hold files used to create figures and diagrams

dir.create("../outputs/007_createFiguresDiagrams")

# - +Create directory to hold text files with sample size summaries

dir.create("../outputs/008_sampleSizeSummaries")

# - Create product directories ----

dir.create("../products")

# - +Create directory to hold figures produced in R scripts

dir.create("../products/figures")
