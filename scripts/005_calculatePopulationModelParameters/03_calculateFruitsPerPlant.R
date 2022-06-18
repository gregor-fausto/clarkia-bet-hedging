####
####
# Script to calculate parameters for population model
# fruits per plant
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

tmpDataDirectory = "outputs/001_prepareDataForModels/"
mcmcDirectory = "outputs/005_calculatePopulationModelParameters/01_parameterPosteriorDistributions/"
outputDirectory = "outputs/005_calculatePopulationModelParameters/"

# - +load libraries ----
library(MCMCvis)
library(tidyverse)
library(magrittr)

# - Read in what's needed ----

# - +Read in data ----
tmpDataDirectory <- paste0(tmpDataDirectory,list.files(tmpDataDirectory))
ref.df <- readRDS(tmpDataDirectory[[grep("referenceSeedsPerFruit.RDS",tmpDataDirectory)]])
ref.df.tot <- readRDS(tmpDataDirectory[[grep("referenceTOTFruitsPerPlant.RDS",tmpDataDirectory)]])
ref.df.tfe <- readRDS(tmpDataDirectory[[grep("referenceTFEFruitsPerPlant.RDS",tmpDataDirectory)]])

# - +Read in MCMC samples ----
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))
mcmcSamplesFruits <- readRDS(mcmcSampleDirectory[[grep("fruitsPerPlantParameters.RDS",mcmcSampleDirectory)]])
mcmcSamplesSeeds <- readRDS(mcmcSampleDirectory[[grep("seedsPerFruitParameters.RDS",mcmcSampleDirectory)]])

# - +Get population- and year-level estimate ----

# - ++Total fruit equivalents ----

# annual estimates for total fruit equivalents (2006-2012)
mu_tfe=MCMCchains(mcmcSamplesFruits, params=c("mu_tfe"))

n.iter=dim(mu_tfe)[1]
index.site<-unique(ref.df.tfe$site)
dat.list <- list()
years = 2006:2012
for(i in 1:20){
  index.tmp <- ref.df.tfe[ref.df.tfe$site==index.site[i],]
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  for(j in 1:length(years)){
    if((years[j] %in% index.tmp$year)){
      index.tmp.tmp = as.numeric(index.tmp$siteYearIndex_tfe[index.tmp$year==years[j]])
      mat[,j] = mu_tfe[,index.tmp.tmp]
    }
  }
  dat.list[[i]] <- mat
  
}

saveRDS(dat.list,paste0(outputDirectory,"tfe-population-year-level-mat.RDS"))


# - ++Total and damaged fruits ----

# annual estimates for total fruits per plant (2013-2020)
mu_tot=MCMCchains(mcmcSamplesFruits, params=c("mu_tot"))
# annual estimates for proportion of fruits damaged (2013-2020)
p_dam=boot::inv.logit(MCMCchains(mcmcSamplesFruits, params=c("mu")))

# product of total fruits per plant and proportion damaged gives average number of damaged fruits per plant
mu_dam = mu_tot*p_dam
# product of total fruits per plant and 1- proportion damaged gives average number of undamaged fruits per plant
mu_und = mu_tot*(1-p_dam)


n.iter=dim(mu_und)[1]
index.site<-unique(ref.df$site)
dat.list_und <- dat.list_dam <- dat.list_tot <- dat.list_pdam <- list()
years = 2013:2020
for(i in 1:20){
  index.tmp <- ref.df.tot[ref.df.tot$site==index.site[i],]
  mat_und = mat_dam = mat_tot = mat_pdam = matrix(NA,nrow = n.iter,ncol = 8)
  for(j in 1:length(years)){
    if((years[j] %in% index.tmp$year)){
      index.tmp.tmp = as.numeric(index.tmp$siteYearIndex_tot[index.tmp$year==years[j]])
      mat_tot[,j] = mu_tot[,index.tmp.tmp]
      mat_und[,j] = mu_und[,index.tmp.tmp]
      mat_pdam[,j] = p_dam[,index.tmp.tmp]
      mat_dam[,j] = mu_dam[,index.tmp.tmp]
    }
  }
  dat.list_tot[[i]] <- mat_tot
  dat.list_und[[i]] <- mat_und
  dat.list_pdam[[i]] <- mat_pdam
  dat.list_dam[[i]] <- mat_dam
  
}

saveRDS(dat.list_tot,paste0(outputDirectory,"f_tot-population-year-level-mat.RDS"))
saveRDS(dat.list_pdam,paste0(outputDirectory,"f_pdam-population-year-level-mat.RDS"))


# - ++Seeds per undamaged and damaged fruits ----

mu_seeds=(MCMCchains(mcmcSamplesSeeds, params=c("mu_seeds")))
# proportion of seeds that are eaten, on average
p_dam_seeds=boot::inv.logit(MCMCchains(mcmcSamplesSeeds, params=c("mu")))

# number of seeds in a damaged fruit
mu_dam_seeds=mu_seeds[,as.numeric(ref.df$siteYearIndex_und)]*p_dam_seeds[,as.numeric(ref.df$siteYearIndex_dam)]

# - +Calculate a composite to recreate TFE ----

# calculate the ratio of seeds in damaged:undamaged fruits
ratio.seeds=mu_dam_seeds/mu_seeds

# convert to matrix
ref.df.tmp <- ref.df %>%
  dplyr::filter(year>2012)

n.iter=dim(mu_und)[1]
index.site<-unique(ref.df$site)
dat.list_prop <- list()
years = 2013:2020
for(i in 1:20){
  index.tmp <- ref.df.tmp[ref.df.tmp$site==index.site[i],]
  mat_prop = matrix(NA,nrow = n.iter,ncol = 8)
  for(j in 1:length(years)){
    if((years[j] %in% index.tmp$year)){
      index.tmp.tmp = as.numeric(index.tmp$siteYearIndex_und[index.tmp$year==years[j]])
      mat_prop[,j] = ratio.seeds[,index.tmp.tmp]
    }
  }
  dat.list_prop[[i]] <- mat_prop

}

# now fill in any missing years

df.test <- data.frame(ref.df,p=ratio.seeds[1,])
df.test2 <- data.frame(ref.df.tot,p2 = mu_dam[1,])

df.joined <- df.test2 %>% dplyr::left_join(df.test,by=c('site','year')) 

# missing undamaged
df.joined[is.na(df.joined$siteYearIndex_und),]

# calculate for LO
mu0_seeds=exp((MCMCchains(mcmcSamplesSeeds, params=c("nu_seeds")))[,13])

# proportion of seeds that are eaten, on average
p0_dam_seeds=boot::inv.logit(MCMCchains(mcmcSamplesSeeds, params=c("mu0")))

# number of seeds in a damaged fruit
mu_lo_dam_seeds=mu0_seeds*p0_dam_seeds[,13]

# - +Calculate a composite to recreate TFE ----

# calculate the ratio of seeds in damaged:undamaged fruits
ratio.seeds.lo=mu_lo_dam_seeds/mu0_seeds

dat.list_prop[[13]][,3] = ratio.seeds.lo

# fill in missing damaged with average numbers

# missing undamaged
dam.missing<-df.joined[is.na(df.joined$siteYearIndex_dam),]

# proportion of seeds that are eaten, on average
p0_dam_seeds=boot::inv.logit(MCMCchains(mcmcSamplesSeeds, params=c("mu0")))

index.site<-unique(ref.df$site)
index.year<-(unique(ref.df$year))[-(1:7)]
index.ref<-(unique(ref.df$siteYearIndex_und))

for(i in 1:dim(dam.missing)[1]){
  site.ref = (index.site) %in% dam.missing$site[i]
  avg.p.dam = p0_dam_seeds[,site.ref]
  year.ref = (index.year) %in% dam.missing$year[i]
  if(is.na(dam.missing$siteYearIndex_und[i])){} else {
  index.ref.tmp = (index.ref) %in% dam.missing$siteYearIndex_und[i]
  
  # number of seeds in a damaged fruit
  mu_dam_seeds.tmp=mu_seeds[,as.numeric(index.ref.tmp)]*avg.p.dam
  
  # - +Calculate a composite compatible with TFE ----
  
  # calculate the ratio of seeds in damaged:undamaged fruits
  ratio.seeds.tmp=mu_dam_seeds.tmp/mu_seeds[,as.numeric(index.ref.tmp)]
  
  dat.list_prop[[(1:20)[site.ref]]][,year.ref] <- ratio.seeds.tmp
  }
}

# - +Calculate a composite to make compatible with TFE ----

# use this ratio to adjust to create a composite estimate for years 2013-2020
# that is comparable to the estimates of total fruit equivalents per plant

years = 2006:2020
for(i in 1:20){

  und.mat <-dat.list_und[[i]]
  dam.mat <-dat.list_dam[[i]]
  prop.mat <- dat.list_prop[[i]]
  tfe.composite <- und.mat + dam.mat*prop.mat
  dat.list[[i]][,8:15] <- tfe.composite
  
}

saveRDS(dat.list,paste0(outputDirectory,"combinedF-population-year-level-mat.RDS"))
