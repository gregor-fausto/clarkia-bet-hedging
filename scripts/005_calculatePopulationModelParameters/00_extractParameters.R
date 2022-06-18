####
####
# Script to split posterior samples into parameters only
# to speed up working with them in later scripts
# the full MCMC object has both samples for the parameters and model checks
# which makes it very large; this script extracts only parameters
# after model checks have been done
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

mcmcDirectory = "outputs/002_fitStatisticalModels/mcmcSamples/"
outputDirectory = "outputs/005_calculatePopulationModelParameters/01_parameterPosteriorDistributions/"

# - +load libraries ----
library(MCMCvis)
library(tidyverse)
library(magrittr)

# - Read in what's needed ----

# - +Read in MCMC samples ----
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))

# - Read in MCMC samples for seed bag experiment (field) ----

mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedBagExperimentPosteriorSamples.RDS",mcmcSampleDirectory)]])

parameters = c("mu0_s","sigma0_s",
               "mu0_g","sigma0_g",
               "mu0_s0","sigma0_s0",
               "mu_s","mu_g","mu_s0")

mcmcSamples.pars <- MCMCchains(mcmcSamples,params=parameters)

saveRDS(mcmcSamples.pars,paste0(outputDirectory,"seedBagExperimentParameters.RDS"))

rm(mcmcSamples)
rm(mcmcSamples.pars)

# - Read in MCMC samples for viability trials (lab) ----

mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("viabilityPosteriorSamples.RDS",mcmcSampleDirectory)]])

parameters = c("mu0_g","sigma0_g",
               "mu0_v","sigma0_v",
               "mu_g","mu_v")

mcmcSamples.pars <- MCMCchains(mcmcSamples,params=parameters)

saveRDS(mcmcSamples.pars,paste0(outputDirectory,"viabilityTrialsParameters.RDS"))

rm(mcmcSamples)
rm(mcmcSamples.pars)

# - Read in MCMC samples for seedling survival to fruiting ----

mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedlingSurvivalPosteriorSamples.RDS",mcmcSampleDirectory)]])

parameters = c("mu0","sigma0","mu")

mcmcSamples.pars <- MCMCchains(mcmcSamples,params=parameters)

saveRDS(mcmcSamples.pars,paste0(outputDirectory,"seedlingSurvivalParameters.RDS"))

rm(mcmcSamples)
rm(mcmcSamples.pars)

# - Read in MCMC samples for fruits per plant ----

mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("fruitsPerPlantPosteriorSamples.RDS",mcmcSampleDirectory)]])

parameters = c("nu_tfe","sigma0_tfe",
               "sigma_overdisp_tfe",
               "nu_tot","sigma0_tot",
               "sigma_overdisp_tot",
               "mu0","sigma0",
               "mu_tfe","mu_tot","mu")

mcmcSamples.pars <- MCMCchains(mcmcSamples,params=parameters)

saveRDS(mcmcSamples.pars,paste0(outputDirectory,"fruitsPerPlantParameters.RDS"))

rm(mcmcSamples)
rm(mcmcSamples.pars)

# - Read in MCMC samples for seeds per fruit ----

mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedsPerFruitPosteriorSamples.RDS",mcmcSampleDirectory)]])

parameters = c("nu_seeds","sigma0_seeds",
               "sigma_overdisp_seeds",
               "mu0","sigma0",
               "mu_seeds","mu")

mcmcSamples.pars <- MCMCchains(mcmcSamples,params=parameters)

saveRDS(mcmcSamples.pars,paste0(outputDirectory,"seedsPerFruitParameters.RDS"))

rm(mcmcSamples)
rm(mcmcSamples.pars)


# - Read in MCMC samples for seed bag+pot experiment (field) ----

mcmcDirectory = "~/Dropbox/clarkia-bet-hedging/outputs/002_fitStatisticalModels/mcmcSamples/"
outputDirectory = "~/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/01_parameterPosteriorDistributions/"

mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))

mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedBagPotExperimentsPosteriorSamples-s0.RDS",mcmcSampleDirectory)]])

parameters = c("mu0_s","sigma0_s",
               "mu0_g","sigma0_g",
               "mu0_s0","sigma0_s0",
               "mu_s","mu_g","mu_s0")

mcmcSamples.pars <- MCMCchains(mcmcSamples,params=parameters)

saveRDS(mcmcSamples.pars,paste0(outputDirectory,"seedBagPotExperimentParameters.RDS"))

rm(mcmcSamples)
rm(mcmcSamples.pars)