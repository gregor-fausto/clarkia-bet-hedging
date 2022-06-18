# -------------------------------------------------------------------
# Run scripts for model fitting
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Set JAGS parameters and random seed
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
# -------------------------------------------------------------------
n.adapt = 10000/1
n.update = 15000/1
n.iterations = 15000/1
n.thin = 10

dataDirectory = "~/Dropbox/clarkia-data-processing/analysis-data/"
modelDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/models/"
outDataDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/002_statisticalModelFitting/data/"
outMCMCDirectory = "/Users/Gregor/Dropbox/dataLibrary/seed-banks-reanalysis-2/"
scriptsModelFitting = "/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/scripts/002_statisticalModelFitting/"

# -------------------------------------------------------------------
# Run scripts for model fitting
# -------------------------------------------------------------------

# run models with partial pooling
# source(paste0(scriptsModelFitting,"01_modelScriptsSeedsNonparametric-noTransition.R"))
#Sys.sleep(60*10)
source(paste0(scriptsModelFitting,"01_modelScriptsSeedsNonparametric.R"))
source(paste0(scriptsModelFitting,"01_modelScriptsSeedsNonparametric-broad.R"))
source(paste0(scriptsModelFitting,"01_modelScriptsSeedsNonparametric-uniform.R"))

Sys.sleep(60*10)
source(paste0(scriptsModelFitting,"02_modelScriptsViability.R"))
Sys.sleep(60*10)
source(paste0(scriptsModelFitting,"03_modelScriptsSeedlingSurvival.R"))
Sys.sleep(60*10)
source(paste0(scriptsModelFitting,"04_modelScriptsFruits.R"))
Sys.sleep(60*10)
source(paste0(scriptsModelFitting,"05_modelScriptsSeedsPerFruit.R"))

# run models with no pooling
source(paste0(scriptsModelFitting,"06_modelScriptsSeedlingSurvivalNoPool.R"))
Sys.sleep(60*10)
source(paste0(scriptsModelFitting,"07_modelScriptsFruitsNoPool.R"))
Sys.sleep(60*10)
source(paste0(scriptsModelFitting,"08_modelScriptsSeedsPerFruitNoPool.R"))

# -------------------------------------------------------------------
# Run scripts for evaluating convergence 
# -------------------------------------------------------------------
# scriptConvergenceDirectory = "/Users/Gregor/Desktop/scriptsModelConvergence/"
# fileDirectory = "/Users/Gregor/Desktop/mcmcSamples/"
# outputDirectory = "/Users/Gregor/Desktop/convergence/"

scriptConvergenceDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/scripts/003_statisticalModelDiagnostics/"
fileDirectory = "/Users/Gregor/Dropbox/dataLibrary/seed-banks-reanalysis-2/"
outputDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/convergence/"

# check convergence for belowground models
source(paste0(scriptConvergenceDirectory,"01_modelConvergenceSeeds.R"))
source(paste0(scriptConvergenceDirectory,"02_modelConvergenceViability.R"))

# check convergence for partially pooled aboveground models
source(paste0(scriptConvergenceDirectory,"03_modelConvergenceSeedlingSurvival.R"))
source(paste0(scriptConvergenceDirectory,"04_modelConvergenceFruits.R"))
source(paste0(scriptConvergenceDirectory,"05_modelConvergenceSeedsPerFruit.R"))

# check convergence for no pooling aboveground models
source(paste0(scriptConvergenceDirectory,"06_modelConvergenceSeedlingSurvivalNoPool.R"))
source(paste0(scriptConvergenceDirectory,"07_modelConvergenceFruitsNoPool.R"))
source(paste0(scriptConvergenceDirectory,"08_modelConvergenceSeedsPerFruitNoPool.R"))

# pause script here to check that convergence has been achieved

# -------------------------------------------------------------------
# Run scripts for model checks
# -------------------------------------------------------------------
scriptConvergenceDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/scripts/004_statisticalModelChecking/"
dataDirectory = "~/Dropbox/clarkia-data-processing/analysis-data/"
outDataDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/002_statisticalModelFitting/data/"
fileDirectory = "/Users/Gregor/Dropbox/dataLibrary/seed-banks-reanalysis-2/"
outputDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/modelChecks-2022/"

# source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/modelChecksBelowground.R")
# source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/03_modelChecksSeedlingSurvival.R")
# source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/04_modelChecksFruits.R")
# source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/05_modelChecksSeedsPerFruit.R")

# following scripts not yet written
# source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsModelChecking/modelChecksViability.R")

# -------------------------------------------------------------------
# Run scripts to output summary tables and summary plots
# -------------------------------------------------------------------

# need to consider what the right scripts are for here

# # -------------------------------------------------------------------
# # Run scripts to recover structured model parameters (derived quantities)
# # -------------------------------------------------------------------
scriptConvergenceDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/scripts/004_statisticalModelChecking/"
dataDirectory = "~/Dropbox/clarkia-data-processing/analysis-data/"
outDataDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/002_statisticalModelFitting/data/"
fileDirectory = "/Users/Gregor/Dropbox/dataLibrary/seed-banks-reanalysis-2/"
outputDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/products/figures/modelChecks-2022/"

# 
# source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsParameters/parametersBelowground.R")
# source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsParameters/parametersSeedlingSurvival.R")
# source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsSummaries/parametersFruitsPerPlant.R")
# source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsSummaries/parametersSeedsPerFruit.R")

# -------------------------------------------------------------------
# Run scripts to recover structured model parameters (derived quantities)
# -------------------------------------------------------------------
fileDirectory = "/Users/Gregor/Dropbox/dataLibrary/seed-banks-reanalysis-2/"
outDataDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/002_statisticalModelFitting/data/"

# source("/Users/gregor/Dropbox/clarkiaSeedBanks/analysis/scripts/006_derivedQuantities/01_derivedQuantitiesSeedSurvivalGermination.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/analysis/scripts/006_derivedQuantities/02_derivedQuantitiesSeedlingSurvival.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/analysis/scripts/006_derivedQuantities/03_derivedQuantitiesFruits.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/analysis/scripts/006_derivedQuantities/04_derivedQuantitiesSeedsPerFruit.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/analysis/scripts/006_derivedQuantities/05_calculateReproductiveSuccess.R")

# -------------------------------------------------------------------
# Run scripts to analyze data for paper
# -------------------------------------------------------------------

# create figure summarizing aboveground observations
# generate RDS with populations and years without plants
source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/scripts/zeroFitness.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/scripts/calculateReproductiveSuccess.R")

source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/scripts/estimateCorrelationGerminationSurvival-Reanalysis.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/scripts/estimateCorrelationGerminationReproductiveSuccess-Reanalysis.R")
source("/Users/gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/scripts/simulatePopulationGrowthRateDIModel-Reanalysis.R")
