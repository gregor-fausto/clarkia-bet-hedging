####
####
# Script to calculate parameters for population model
# Seedling survival
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

mcmcDirectory = "outputs/005_calculatePopulationModelParameters/01_parameterPosteriorDistributions/"
outputDirectory = "outputs/005_calculatePopulationModelParameters/"
tmpDataDirectory = "outputs/001_prepareDataForModels/"

# - +load libraries ----
library(MCMCvis)
library(tidyverse)
library(magrittr)


# - Read in what's needed ----

# - +Read in MCMC samples ----
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedlingSurvivalParameters.RDS",mcmcSampleDirectory)]])

# - +Read in site names ----
siteAbiotic <- read.csv("data/siteAbioticData.csv",header=TRUE)
siteNames <- siteAbiotic$site

# - +Get population-level estimate ----

mu0Pooled <- MCMCchains(mcmcSamples,params="mu0")
sigma0Pooled <- apply(mu0Pooled,2,boot::inv.logit)

saveRDS(sigma0Pooled,paste0(outputDirectory,"sigma-population-level.RDS"))

# - +Get population- and year-level estimate ----

muPooled <- MCMCchains(mcmcSamples,params="mu")
sigmaPooled <- apply(muPooled,2,boot::inv.logit)

saveRDS(sigmaPooled,paste0(outputDirectory,"sigma-population-year-level.RDS"))

# - +convert to matrix ----

refDf <- readRDS(paste0(tmpDataDirectory,"seedlingRef.RDS"))

refDf<-refDf %>%
  dplyr::select(site,year,siteYearIndex) %>%
  unique

n.iter=dim(sigmaPooled)[1]
index.site<-unique(refDf$site)
dat.list <- list()
years = 2006:2020
for(i in 1:20){
  index.tmp <- refDf[refDf$site==index.site[i],]
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  for(j in 1:15){
    if((years[j] %in% index.tmp$year)){
      index.tmp.tmp = as.numeric(index.tmp$siteYearIndex[index.tmp$year==years[j]])
      mat[,j] = sigmaPooled[,index.tmp.tmp]
    }
  }
  dat.list[[i]] <- mat
  
}

saveRDS(dat.list,paste0(outputDirectory,"sigma-population-year-level-mat.RDS"))

