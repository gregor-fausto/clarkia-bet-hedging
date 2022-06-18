####
####
# Script to fit model for belowground vital rates
# germination
####
####

# - Environment ----
outDataDirectory = "outputs/002_fitStatisticalModels/data/"

# - Libraries ----
library(tidyverse)
library(parallel)
library(stringr)
library(bbmle)

# - Likelihood function ----

source("scripts/009_compareBayesianLikelihood/01_functions.R")

# - Optimization ----
# with bbmle::mle2 

siteNumber = 1:20
modelList = list()
modelList.ci = list()
dataList = list()
for(i in 1:20){

  # Prep/filter data 
  siteNumbers = siteNumber[i]
  
  data <- readRDS(file=paste0(outDataDirectory,"viabilityData.RDS"))
  
  dat = list()
  st = data$siteViab%in%siteNumbers
  dat$siteViab = data$siteViab[st]; dat$germStart = data$germStart[st]; dat$germCount = data$germCount[st]
  dat$indexViab = data$indexViab[st];
  st_v = data$siteViab_v%in%siteNumbers
  dat$viabStart = data$viabStart[st_v]; dat$viabStain = data$viabStain[st_v]; 
  dat$indexViab_v = data$indexViab_v[st_v];
 
  dat$n1 = length(data$germStart)
  dat$n2 = length(data$viabStart)

  # - +create data frame ----
  
  f = function(x){
    tmp = sum(is.na(x))>0
    return(tmp)
  }
  
  if(f(dat$germStart)|f(dat$germCount)){
    t.vec = !is.na(data$germCount)
    dat$germStart = dat$germStart[t.vec]
    dat$germCount = dat$germCount[t.vec]
  }
  
  if(f(dat$viabStart)|f(dat$viabStain)){
    t.vec = !is.na(data$viabStain)
    dat$viabStain = dat$viabStain[t.vec]
    dat$viabStart = dat$viabStart[t.vec]
  }
  
  dataTmp = list(dat$germStart,dat$germCount,
                 dat$viabStart,dat$viabStain,
                 dat$indexViab,dat$indexViab_v)
  
  dataList[[i]] = dataTmp
  
  # set up optimization
  parnames(likViabilityModel) =c(paste0("g",1:6),paste0("v",1:6))
  
  # fit model
  m1 = mle2(likViabilityModel, start = c(g1 = .5, g2 = .5, g3 = .5, 
                                       g4 = .5, g5 = .5, g6 = .5,
                                       v1=.5,v2=.5,v3=.5,v4=.5,
                                       v5=.5,v6=.5 ) ,
            data = list( dat = dataTmp ),
            method="L-BFGS-B", lower=rep(0.00001,12),upper=rep(.9999999,12))
  
  modelList[[i]] = m1

}


saveRDS(modelList,file=paste0("outputs/009_compareBayesianLikelihood/viabilityTrialsLikelihoodModels.rds"))
saveRDS(dataList,file=paste0("outputs/009_compareBayesianLikelihood/viabilityTrialsLikelihoodData.rds"))
