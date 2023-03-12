####
####
# Script to calculate derived quantities for per-capita reproductive success
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

# - Function to calculate mode ----

posterior.mode = function(x){
  if(!is.na(x[1])){ x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]}else if(is.na(x[1])){
    modeParam <- NA
  }
  return(modeParam)
}

# - Directory ----

inputDirectoryPopulation = "outputs/005_calculatePopulationModelParameters/02_populationModelParameters/"
inputDirectoryMatrix = "outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/"
outputDirectory = "outputs/005_calculatePopulationModelParameters/04_reproductiveSuccess/"

# - +load libraries ----
library(MCMCvis)
library(tidyverse)
library(magrittr)

# - Read in what's needed for plotting ----

# - +Read in data ----
inputDir <- paste0(inputDirectoryMatrix,list.files(inputDirectoryMatrix))
inputDirPop <- paste0(inputDirectoryPopulation,list.files(inputDirectoryPopulation))

# - Model with population-level pooling ----

# - +Get derived quantities ----
# - ++seedling survival ----
sigma = readRDS(inputDir[grep("sigma-population-year-level-mat.RDS",inputDir)])
# - ++total fruit equivalents and composite fruit equivalents per plant ----
fruits <- readRDS(inputDir[grep("combinedF-population-year-level-mat.RDS",inputDir)])
# - ++seeds per undamaged fruit ----
seeds <- readRDS(inputDir[grep("phi-population-year-level-mat.RDS",inputDir)])

# - population-level estimates----
sigma0 <- readRDS(inputDirPop[grep("sigma-population-level.RDS",inputDirPop)])
phi0 <- readRDS(inputDirPop[grep("phi-population-level.RDS",inputDirPop)])

# - site information----

siteAbiotic <- read.csv("data/siteAbioticData.csv",header=TRUE)

position<-siteAbiotic %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

# get mode of estimates to summarize
sigma.mode <- unlist(lapply(sigma,apply,2,posterior.mode))
fec.mode <- unlist(lapply(fruits,apply,2,posterior.mode))
phi.mode <- unlist(lapply(seeds,apply,2,posterior.mode))

tmp.df=data.frame(site=rep(position$site,each=15),
                  year=rep(2006:2020,20),
                  cbind(sigma.mode,fec.mode,phi.mode) ) %>%
  dplyr::mutate(perCapitaReproductiveSuccess.mode = sigma.mode*fec.mode*phi.mode)

write.csv(tmp.df,paste0(outputDirectory,"reproductiveSuccess.csv"),row.names = FALSE)

# - Calculate per-capita reproductive success as derived quantity, after adustments ----

perCapitaRS.list <- list()
for(i in 1:20){
  perCapitaRS.mat = matrix(NA,nrow=dim(sigma[[1]])[1],ncol=15)
  for(j in 1:15){
    perCapitaRS.mat[,j] <- sigma[[i]][,j]*fruits[[i]][,j]*seeds[[i]][,j]
  }
  perCapitaRS.list[[i]] <- perCapitaRS.mat
}

saveRDS(perCapitaRS.list,paste0(outputDirectory,"reproductiveSuccess-populationYear-mat.RDS"))

# - Handling missing data ----
# I followed rules 1-4 below to deal with missing observations
# if
# 1. if there are no seedlings or fruiting plants at all (n=3 population-years)
# 2. if there are seedlings but no fruiting plants at all (n=9 population-years)
# I assumed sigma (seedlings survival) was the population mean from all years with observations
# I assumed fruits per plant (F) and seeds per fruit (phi) were 0
# if
# 3. in 1 year at 1 population (LO, 2015) there was 1 plant with 1 fruit and the fruit was not collected
# in this case I substituted the population average seeds/fruit
# if
# 4. if there are no seedlings in permanent plots but plants elsewhere in the populaton (n=13 population-years)
# I estimated seedling survival to fruiting as the population average from all years with observations

# - Adjust missing data ----
# the loop below applies the rules above
# and assigns 0s or averages to missing data

for( i in 1:20 ){
 
    # take the first row of each object
    sigma.tmp = sigma[[i]][1,]
    f.tmp = fruits[[i]][1,]
    phi.tmp = seeds[[i]][1,]
    vr.tmp <-cbind(sigma=sigma.tmp,f=f.tmp,phi=phi.tmp)
    
    # Case 1: all values are missing
    allMissing <- apply(is.na(vr.tmp),1,sum)==3
    
    # assign population-level estimate for seedling survival
    sigma[[i]][,allMissing] <-  sigma0[,i]
    
    # assign 0 to fruits per plant and seeds per fruit
    fruits[[i]][,allMissing] <- 0
    seeds[[i]][,allMissing] <- 0
    
    # updated object
    sigma.tmp = sigma[[i]][1,]
    f.tmp = fruits[[i]][1,]
    phi.tmp = seeds[[i]][1,]
    vr.tmp <-cbind(sigma=sigma.tmp,f=f.tmp,phi=phi.tmp)
    
    # Case 2: fruits per plant and seeds per fruit are missing
    zeroSurvival <- !is.na(vr.tmp[,1])&(apply(is.na(vr.tmp[,2:3]),1,sum)==2)
    
    # assign 0 to fruits per plant and seeds per fruit
    fruits[[i]][,zeroSurvival] <- 0
    seeds[[i]][,zeroSurvival] <- 0
    
    # take the first row of updated object
    sigma.tmp = sigma[[i]][1,]
    f.tmp = fruits[[i]][1,]
    phi.tmp = seeds[[i]][1,]
    vr.tmp <-cbind(sigma=sigma.tmp,f=f.tmp,phi=phi.tmp)
    
    # Case 3: seeds per fruit is missing
    missingSeeds <- is.na(vr.tmp[,3])&(apply(is.na(vr.tmp[,1:2]),1,sum)==0)
    
    # assign population-level estimate to seeds per fruit
    seeds[[i]][,missingSeeds] <- phi0[,i]
    
    # take the first row of updated object
    sigma.tmp = sigma[[i]][1,]
    f.tmp = fruits[[i]][1,]
    phi.tmp = seeds[[i]][1,]
    vr.tmp <-cbind(sigma=sigma.tmp,f=f.tmp,phi=phi.tmp)
    
    # Case 4: seedling survival to fruiting is missing
    missingSurvival <- is.na(vr.tmp[,1])&(apply(is.na(vr.tmp[,2:3]),1,sum)==0)
    
    # assign population-level estimate to seedling survival
    sigma[[i]][,missingSurvival] <- sigma0[,i]
    
}

saveRDS(sigma,paste0(outputDirectory,"sigmaWithCorrectionForMissingness-populationYear-mat.RDS"))
saveRDS(fruits,paste0(outputDirectory,"fecWithCorrectionForMissingness-populationYear-mat.RDS"))
saveRDS(seeds,paste0(outputDirectory,"phiWithCorrectionForMissingness-populationYear-mat.RDS"))

# - write out estimates with adjustments for missing data ----

sigma.mode <- unlist(lapply(sigma,apply,2,posterior.mode))
fec.mode <- unlist(lapply(fruits,apply,2,posterior.mode))
phi.mode <- unlist(lapply(seeds,apply,2,posterior.mode))

tmp.df=data.frame(site=rep(position$site,each=15),
                  year=rep(2006:2020,20),
                  cbind(sigma.mode,fec.mode,phi.mode) ) %>%
  dplyr::mutate(perCapitaReproductiveSuccess.mode = sigma.mode*fec.mode*phi.mode)

tmp.df[,3:6] <- apply(tmp.df[,3:6],2,signif,4)

write.csv(tmp.df,paste0(outputDirectory,"reproductiveSuccessWithCorrectionForMissingness.csv"),row.names = FALSE)

# - Calculate per-capita reproductive success as derived quantity, after adustments ----

perCapitaRS.list <- list()
for(i in 1:20){
  perCapitaRS.mat = matrix(NA,nrow=dim(sigma[[1]])[1],ncol=15)
  for(j in 1:15){
    perCapitaRS.mat[,j] <- sigma[[i]][,j]*fruits[[i]][,j]*seeds[[i]][,j]
  }
  perCapitaRS.list[[i]] <- perCapitaRS.mat
}


saveRDS(perCapitaRS.list,paste0(outputDirectory,"reproductiveSuccessWithCorrectionForMissingness-populationYear-mat.RDS"))
