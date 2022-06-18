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
  
  data <- readRDS(file=paste0(outDataDirectory,"seedData.RDS"))
  
  st = data$siteSurvival%in%siteNumbers
  data$siteSurvival = data$siteSurvival[st]; data$seedStart = data$seedStart[st]; data$y = data$y[st]
  data$months = data$months[st]; data$yearSurvival = data$yearSurvival[st]; data$compIndex = data$compIndex[st]
  data$n_siteSurvival = length(siteNumbers)
  
  # Recode variables
  data$roundSurvival = ifelse(data$compIndex[st]%in%c(1:6),1,ifelse(data$compIndex[st]%in%c(7:10),2,3))
  
  # Filter data
  st = data$siteGermination%in%siteNumbers
  data$siteGermination = data$siteGermination[st]; data$totalJan = data$totalJan[st]; 
  data$seedlingJan = data$seedlingJan[st]; data$germinationIndex = data$germinationIndex[st]
  data$n_siteGermination = length(siteNumbers)
  
  data$round = ifelse(data$germinationIndex[st]%in%c(1,2,3),1,ifelse(data$germinationIndex[st]%in%c(4,5),2,3))
  
  data$binary = ifelse(data$germinationIndex==1,1,
                       ifelse(data$germinationIndex==2,2,
                              ifelse(data$germinationIndex==3,3,
                                     ifelse(data$germinationIndex==4,4,
                                            ifelse(data$germinationIndex==5,3,5)))))
  
  st = data$sitePlot%in%siteNumbers
  data$sitePlot = data$sitePlot[st]; data$yearPlot = data$yearPlot[st]; data$fec = data$fec[st]
  data$plotSeedlings = data$plotSeedlings[st]; data$fecIndex = data$fecIndex[st]; data$fecCompIndex = data$fecCompIndex[st]
  
  data$n1 = length(data$siteSurvival)
  data$n2 = length(data$siteGermination)
  data$n3 = length(data$sitePlot)
  
  # - +create data frame ----
  
  f = function(x){
    tmp = sum(is.na(x))>0
    return(tmp)
  }
  
  if(f(data$seedStart)|f(data$y)){
    t.vec = !is.na(data$y)
    data$seedStart = data$seedStart[t.vec]
    data$y = data$y[t.vec]
    data$compIndex = data$compIndex[t.vec]
  }
  
  if(f(data$totalJan)|f(data$seedlingJan)){
    t.vec = !is.na(data$seedlingJan)
    data$seedlingJan = data$seedlingJan[t.vec]
    data$totalJan = data$totalJan[t.vec]
    data$germinationIndex = data$germinationIndex[t.vec]
  }
  
  dataTmp = list(data$seedStart,data$y,
                 data$totalJan,data$seedlingJan,
                 data$compIndex,
                 data$germinationIndex)
  
  dataList[[i]] = dataTmp
  
  # set up optimization
  parnames(likSeedBagModel) =c(paste0("g",1:6),paste0("s",1:12))
  
  # fit model
  m1 = mle2(likSeedBagModel, start = c(g1 = .5, g2 = .5, g3 = .5, 
                                       g4 = .5, g5 = .5, g6 = .5,
                                       s1=.5,s2=.5,s3=.5,s4=.5,
                                       s5=.5,s6=.5,s7=.5,s8=.5,
                                       s9=.5,s10=.5,s11=.5,s12=.5) ,
            data = list( dat = dataTmp ),
            method="L-BFGS-B", lower=rep(0.00001,18),upper=rep(.9999999,18))
  
  modelList[[i]] = m1

}

saveRDS(modelList,file=paste0("outputs/009_compareBayesianLikelihood/seedSamplesLikelihoodModels.rds"))
saveRDS(dataList,file=paste0("outputs/009_compareBayesianLikelihood/seedSamplesLikelihoodData.rds"))

# - Clear graphics device ----

dev.off()

# - Extract parameter estimates ----
empty.list = list()
for(i in 1:20){
  m.tmp = modelList[[i]]
  empty.list[[i]] <- rbind(data.frame(site = i,parameter = paste0("g",1:6),estimate = m.tmp@coef[1:6]),
  data.frame(site = i,parameter = paste0("s",1:12),estimate = m.tmp@coef[7:18]))
}

estimates=do.call(rbind,empty.list)

# - Join by site name ----

climate <- readRDS("~/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
position<-climate %>% 
  dplyr::ungroup() %>%
  dplyr::select(site,easting,intenseDemography) %>%
  dplyr::mutate(easting=easting/1000)

position = position %>% dplyr::filter(intenseDemography==1) %>% unique
estimates <- estimates %>% dplyr::left_join(data.frame(siteName=position$site,site = 1:20),by='site')

saveRDS(estimates,file=paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","seedSamplesLikelihood.rds"))

