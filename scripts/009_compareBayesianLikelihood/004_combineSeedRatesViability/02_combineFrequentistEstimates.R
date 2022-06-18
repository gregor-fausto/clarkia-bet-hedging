
# - Libraries ----
library(tidyverse)
library(parallel)
library(stringr)
library(bbmle)

# - Read in estimates ----

viabilityList <- readRDS(paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","viabilityEstimates.rds"))
modelList <- readRDS(paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","seedSamplesLikelihoodModels.rds"))

# - Read in additional data ----

climate <- readRDS("~/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
position<-climate %>% 
  dplyr::ungroup() %>%
  dplyr::select(site,easting,intenseDemography) %>%
  dplyr::mutate(easting=easting/1000)

position = position %>% dplyr::filter(intenseDemography==1) %>% unique

orderedIndex = order(position$easting,decreasing=FALSE)
siteNames = unique(position$site)

# - Calculate rates ----

germinationList <- list()
survivalList <- list()

for(i in 1:20){

  # observations
  model.tmp = modelList[i][[1]]
  viab.tmp = viabilityList[[i]]
  
  # viability
  # 1, 4, 6 are round 1
  # 2, 5 are round 2
  # 3 is round 3
  
  # get parameter estimates
  g_est = model.tmp@coef[1:6]
  theta_est = model.tmp@coef[7:18]
  nu_est = viab.tmp[c(1,4,6,2,5,3)]
  
  # calculations for round 1
  
  nu_ratio1 = ifelse(nu_est[2]/nu_est[1]>1,1,nu_est[2]/nu_est[1])
  nu_ratio2 = ifelse(nu_est[3]/nu_est[2]>1,1,nu_est[3]/nu_est[2])
  
  s1.1 = theta_est[1] * (g_est[1] + (1-g_est[1])*(nu_est[1]^(1/3)))
  g1.1 = g_est[1]/(1-(1-(nu_est[1]^(1/3)))*(1-g_est[1]))
  s2.1 = theta_est[2] * ((nu_est[1]^(2/3)))
  s3.1 = ((theta_est[1] * (1-g_est[1]) *theta_est[2]*theta_est[3]) * (g_est[2] + (1-g_est[2])*(nu_est[1]*(nu_ratio1)^(1/3))))/(s1.1*(1-g1.1)*s2.1)
  
  g2.1 = g_est[2]/(1-(1-(nu_est[1]*(nu_ratio1)^(1/3)))*(1-g_est[2]))
  s4.1 = theta_est[4] * (nu_est[1]*(nu_ratio1)^(2/3))
  s5.1 = ((theta_est[1] * (1-g_est[1]) *theta_est[2]*theta_est[3]* (1-g_est[2])*theta_est[4]*theta_est[5]) * (g_est[3] + (1-g_est[3])*(nu_est[2]*(nu_ratio2)^(1/3))))/(s1.1*(1-g1.1)*s2.1*s3.1*(1-g2.1)*s4.1)
  
  g3.1 = g_est[3]/(1-(1-(nu_est[2]*(nu_ratio2)^(1/3)))*(1-g_est[3]))
  s6.1 = theta_est[6] * (nu_est[2]*(nu_ratio2)^(2/3))
  
  # calculations for round 2
  
  nu_ratio3 = ifelse(nu_est[5]/nu_est[4]>1,1,nu_est[5]/nu_est[4])

  s1.2 = theta_est[7] * (g_est[4] + (1-g_est[4])*(nu_est[4]^(1/3)))
  g1.2 = g_est[4]/(1-(1-(nu_est[4]^(1/3)))*(1-g_est[4]))
  s2.2 = theta_est[8] * ((nu_est[4]^(2/3)))
  s3.2 = ((theta_est[7] * (1-g_est[4]) *theta_est[8]*theta_est[9]) * (g_est[5] + (1-g_est[5])*(nu_est[4]*(nu_ratio3)^(1/3))))/(s1.2*(1-g1.2)*s2.2)
  
  g2.2 = g_est[5]/(1-(1-(nu_est[4]*(nu_ratio3)^(1/3)))*(1-g_est[5]))
  s4.2 = theta_est[10] * (nu_est[4]*(nu_ratio3)^(2/3))

  # calculations for round 3
  s1.3 = theta_est[11] * (g_est[6] + (1-g_est[6])*(nu_est[6]^(1/3)))
  g1.3 = g_est[6]/(1-(1-(nu_est[6]^(1/3)))*(1-g_est[6]))
  s2.3 = theta_est[12] * ((nu_est[6]^(2/3)))

  index = c(1,2,3,4,5,6,1,2,3,4,1,2)
  round = c(1,1,1,1,1,1,2,2,2,2,3,3)
  age = c(1,1,2,2,3,3,1,1,2,2,1,1)
  survival = c(s1.1,s2.1,s3.1,s4.1,s5.1,s6.1,s1.2,s2.2,s3.2,s4.2,s1.3,s2.3)
  germination = c(g1.1,g2.1,g3.1,g1.2,g2.2,g1.3)
  
  survivalList[[i]] <- survival
  germinationList[[i]] <- germination
}

# bind into matrix
germination.dat<-do.call(rbind,germinationList)
survival.dat<-do.call(rbind,survivalList)


survival.dat[survival.dat > 1] <- 1
colnames(survival.dat) = c("s1","s2","s3","s4","s5","s6","s1","s2","s3","s4","s1","s2")
colnames(germination.dat) = c("g1","g2","g3","g1","g2","g1")

likelihoodEstimates <- list(germination.dat,survival.dat)

saveRDS(likelihoodEstimates,paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","seedRateEstimatesFrequentistCombined.rds"))

s1 = apply(survival.dat[,c(1,7,11)],1,mean)
s2 = apply(survival.dat[,c(2,8,12)],1,mean)
s3 = apply(survival.dat[,c(3,9)],1,mean)
s4 = apply(survival.dat[,c(4,10)],1,mean)
s5 = survival.dat[,5]
s6 = survival.dat[,6]

# calculate persistence, on average
g1 = apply(germination.dat[,c(1,4,6)],1,mean)
g2 = apply(germination.dat[,c(2,5)],1,mean)
g3 = germination.dat[,3]

germination.dat <- cbind(g1,g2,g3)
survival.dat <- cbind(s1,s2,s3,s4,s5,s6)

likelihoodEstimates <- list(germination.dat,survival.dat)

saveRDS(likelihoodEstimates,paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","seedRateEstimatesFrequentistCombinedPopulationLevel.rds"))
