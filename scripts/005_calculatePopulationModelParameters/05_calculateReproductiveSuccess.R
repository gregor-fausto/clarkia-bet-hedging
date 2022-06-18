####
####
# Script to calculate derived quantities for per-capita reproductive success
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

outputDirectory = "outputs/005_calculatePopulationModelParameters/"

# - +load libraries ----
library(MCMCvis)
library(tidyverse)
library(magrittr)

# - Read in what's needed for plotting ----

# - +Read in data ----
outputDirectory <- paste0(outputDirectory,list.files(outputDirectory))

# - Model with population-level pooling ----

# - +Get derived quantities ----
# - ++seedling survival ----
sigma = readRDS(outputDirectory[grep("sigma-population-year-level-mat.RDS",outputDirectory)])
# - ++total fruit equivalents and composite fruit equivalents per plant ----
fruits <- readRDS(outputDirectory[grep("combinedF-population-year-level-mat.RDS",outputDirectory)])
# - ++seeds per undamaged fruit ----
seeds <- readRDS(outputDirectory[grep("phi-population-year-level-mat.RDS",outputDirectory)])

# - +Check matrix dimensions ----

sigma0 <- readRDS(outputDirectory[grep("sigma-population-level.RDS",outputDirectory)])


# get mode of estimates to summarize

posterior.mode = function(x){
  if(!is.na(x[1])){ x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]}else if(is.na(x[1])){
    modeParam <- NA
  }
  return(modeParam)
}

sigma.mode <- unlist(lapply(sigma,apply,2,posterior.mode))
fec.mode <- unlist(lapply(fruits,apply,2,posterior.mode))
phi.mode <- unlist(lapply(seeds,apply,2,posterior.mode))

siteAbiotic <- read.csv("data/siteAbioticData.csv",header=TRUE)

position<-siteAbiotic %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteIndex <- order(position$easting,decreasing=FALSE)
siteNames = unique(position$site)[siteIndex]

siteNames <- unique(position$site)

tmp.df=data.frame(site=rep(position$site,each=15),
                  year=rep(2006:2020,20),
                  cbind(sigma.mode,fec.mode,phi.mode) ) %>%
  dplyr::mutate(perCapitaReproductiveSuccess = sigma.mode*fec.mode*phi.mode)

write.csv(tmp.df,paste0(outputDirectory,"reproductiveSuccess.csv"),row.names = FALSE)


rs.list <- list()
for( i in 1:20 ){
  mat = matrix(NA,nrow=dim(sigma[[1]])[1],ncol=15)
  for(j in 1:15){
    seed.production<-ifelse(is.na(fruits[[i]][,j]*seeds[[i]][,j]),0,fruits[[i]][,j]*seeds[[i]][,j])
    seedling.survival <- ifelse(is.na(sigma[[i]][,j])&!is.na(fruits[[i]][,j]),sigma0[,i],sigma[[i]][,j])
    mat[,j] <- seedling.survival*seed.production
  }
  rs.list[[i]] <- mat
}

# - +Calculate per-capita reproductive success ----
saveRDS(rs.list,paste0(outputDirectory,"reproductiveSuccess-population-year-level-mat.RDS"))
