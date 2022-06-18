
# - Libraries ----
library(tidyverse)
library(parallel)
library(stringr)
library(bbmle)

# - Read in estimates ----

dataList <- readRDS(paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","viabilityTrialsLikelihoodData.rds"))
modelList <- readRDS(paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","viabilityTrialsLikelihoodModels.rds"))

# - Read in additional data ----

climate <- readRDS("~/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
position<-climate %>% 
  dplyr::ungroup() %>%
  dplyr::select(site,easting,intenseDemography) %>%
  dplyr::mutate(easting=easting/1000)

position = position %>% dplyr::filter(intenseDemography==1) %>% unique

orderedIndex = order(position$easting,decreasing=FALSE)
siteNames = unique(position$site)

# - Obs. vs predicted for germination trials ----
# purple: observations for age 1 seeds
# red: observations for age 2 seeds
# pink: observations for age 3 seeds

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  
  index = orderedIndex[i]
  
  # observations
  data.tmp = dataList[[index]]
  model.tmp = modelList[index][[1]]
  
  # obs vs. predicted for germination
  plot(x=data.tmp[[5]]+rnorm(length(data.tmp[[5]]),0,.05),
       y=data.tmp[[2]]/data.tmp[[1]],
       pch=16,xlim=c(.5,6.5), main='',
       ylab='',xlab='',xaxt='n',yaxt='n',ylim=c(0,1))
  points(1:6,model.tmp@coef[1:6],pch=16,cex=2,
         col=c("purple","purple","purple","red","red","pink"))
  abline(v=c(3.5,5.5),lty='dotted')
  
  mtext(siteNames[index], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}


# - Obs. vs predicted for viability trials ----
# purple: observations for age 1 seeds
# red: observations for age 2 seeds
# pink: observations for age 3 seeds

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  
  index = orderedIndex[i]
  
  # observations
  data.tmp = dataList[[index]]
  model.tmp = modelList[index][[1]]
  
  # obs vs. predicted for germination
  plot(x=data.tmp[[5]]+rnorm(length(data.tmp[[5]]),0,.05),
       y=data.tmp[[4]]/data.tmp[[3]],
       pch=16,xlim=c(.5,6.5), main='',
       ylab='',xlab='',xaxt='n',yaxt='n',ylim=c(0,1))
  points(1:6,model.tmp@coef[7:12],pch=16,cex=2,
         col=c("purple","purple","purple","red","red","pink"))
  abline(v=c(3.5,5.5),lty='dotted')
  
  mtext(siteNames[index], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}


# - Calculate overall probability of viability ----
# purple: observations for age 1 seeds
# red: observations for age 2 seeds
# pink: observations for age 3 seeds

viabilityList <- list()

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  
  # models
  model.tmp = modelList[i][[1]]
  
  germ_estimates <- model.tmp@coef[1:6];
  viab_estimates <- model.tmp@coef[7:12];
  viab_tot <- germ_estimates+viab_estimates*(1-germ_estimates)
  
  viabilityList[[i]] <- viab_tot

  
}




par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  
  # obs vs. predicted for germination
  plot(NA,
       pch=16,xlim=c(.5,6.5), main='',
       ylab='',xlab='',xaxt='n',yaxt='n',ylim=c(0,1))
  
  points(1:3, viabilityEstimate[[i]][c(1,4,6)],type='b')
  points(1:2, viabilityEstimate[[i]][c(2,5)],type='b')
  points(1, viabilityEstimate[[i]][c(3)],type='b')
  
  abline(v=c(3.5,5.5),lty='dotted')
  
  mtext(siteNames[i], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

saveRDS(viabilityList,file=paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","viabilityEstimates.rds"))
