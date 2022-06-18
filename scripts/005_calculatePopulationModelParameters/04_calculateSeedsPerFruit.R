####
####
# Script to calculate parameters for population model
# Seeds per fruit
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

# - Read in what's needed  ----

# - +Read in MCMC samples ----
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))
mcmcSamplesSeeds <- readRDS(mcmcSampleDirectory[[grep("seedsPerFruitParameters.RDS",mcmcSampleDirectory)]])

# - +Read in site names ----
siteAbiotic <- read.csv("data/siteAbioticData.csv",header=TRUE)
siteNames <- siteAbiotic$site

# - +Get population- and year-level estimate ----

# - ++Seeds per undamaged fruit ----

mu_seeds=MCMCchains(mcmcSamplesSeeds, params=c("mu_seeds"))

# - +Save population- and year-level estimate ----

saveRDS(mu_seeds,paste0(outputDirectory,"phi-population-year-level.RDS"))

# - +convert to matrix ----

refDf <- readRDS(paste0(tmpDataDirectory,"referenceSeedsPerFruit.RDS"))

refDf<-refDf %>%
  dplyr::select(site,year,siteYearIndex_und) %>%
  unique

n.iter=dim(mu_seeds)[1]
index.site<-unique(refDf$site)
dat.list <- list()
years = 2006:2020
for(i in 1:20){
  index.tmp <- refDf[refDf$site==index.site[i],]
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  for(j in 1:15){
    if((years[j] %in% index.tmp$year)){
      index.tmp.tmp = as.numeric(index.tmp$siteYearIndex_und[index.tmp$year==years[j]])
      mat[,j] = mu_seeds[,index.tmp.tmp]
    }
  }
  dat.list[[i]] <- mat
  
}

saveRDS(dat.list,paste0(outputDirectory,"phi-population-year-level-mat.RDS"))


# - +Proportion damage ----

mu_dam=MCMCchains(mcmcSamplesSeeds, params=c("mu"))

# - +Save population- and year-level estimate ----

# saveRDS(mu_dam,paste0(outputDirectory,"phi-damage-population-year-level.RDS"))

# - +convert to matrix ----

refDf <- readRDS(paste0(tmpDataDirectory,"referenceSeedsPerFruit.RDS"))

refDf<-refDf %>%
  dplyr::select(site,year,siteYearIndex_dam) %>%
  unique

n.iter=dim(mu_dam)[1]
index.site<-unique(refDf$site)
dat.list <- list()
years = 2006:2020
for(i in 1:20){
  index.tmp <- refDf[refDf$site==index.site[i],]
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  for(j in 1:15){
    if((years[j] %in% index.tmp$year)){
      index.tmp.tmp = as.numeric(index.tmp$siteYearIndex_dam[index.tmp$year==years[j]])
      mat[,j] = mu_dam[,index.tmp.tmp]
    }
  }
  dat.list[[i]] <- mat
  
}

saveRDS(dat.list,paste0(outputDirectory,"phi-damage-population-year-level-mat.RDS"))
