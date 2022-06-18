
# - Libraries ----
library(tidyverse)
library(parallel)
library(stringr)
library(bbmle)

# - Read in estimates ----

dataList <- readRDS(paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","seedSamplesLikelihoodData-constantSurvival.rds"))
modelList <- readRDS(paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","seedSamplesLikelihoodModels-constantSurvival.rds"))

# - Read in additional data ----

climate <- readRDS("~/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
position<-climate %>% 
  dplyr::ungroup() %>%
  dplyr::select(site,easting,intenseDemography) %>%
  dplyr::mutate(easting=easting/1000)

position = position %>% dplyr::filter(intenseDemography==1) %>% unique

orderedIndex = order(position$easting,decreasing=FALSE)
siteNames = unique(position$site)

# - Obs. vs predicted for germination ----

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  
  index = orderedIndex[i]
  
  # observations
  data.tmp = dataList[[index]]
  model.tmp = modelList[index][[1]]
  
  # obs vs. predicted for germination
  plot(x=data.tmp[[6]]+rnorm(length(data.tmp[[6]]),0,.05),
       y=data.tmp[[4]]/data.tmp[[3]],
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



# - Obs. vs predicted for survival ----

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
       pch=16,xlim=c(.5,12.5), main='',
       ylab='',xlab='',xaxt='n',yaxt='n',ylim=c(0,1))
  
  
  surv=model.tmp@coef[7:8]
  surv=c(surv[1],surv[2]^(8/12),surv[2]^(4/12),surv[2]^(8/12),surv[2]^(4/12),surv[2]^(8/12))
  r1=cumprod(surv[1:6])
  r2=cumprod(surv[1:4])
  r3=cumprod(surv[1:2])
  points(1:12,c(r1,r2,r3),pch=16,cex=2,
         col=c("purple","purple","purple","purple","purple","purple",
               "red","red","red","red",
               "pink","pink"))

  abline(v=c(6.5,10.5),lty='dotted')
  
  mtext(siteNames[index], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}


# - Cumulative loss and persistence ----

# vector to hold values after 3 years
persistence.32 = c()

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  
  index = orderedIndex[i]
  
  # model
  model.tmp = modelList[index][[1]]
  
  # obs vs. predicted for germination
  plot(NA,
       pch=16,xlim=c(0,3), main='',
       ylab='',xlab='',xaxt='n',yaxt='n',ylim=c(0,1))
  
  surv=model.tmp@coef[7:8]
  surv=c(surv[1],surv[2]^(8/12),surv[2]^(4/12),surv[2]^(8/12),surv[2]^(4/12),surv[2]^(8/12))
  germ = model.tmp@coef[1:6]
  
  # round 1
  survival.1=(surv[1:6])
  germination.1=germ[1:3]
  notgermination.1 = 1-germination.1
  
  cumPersistence = c(1,survival.1[1],notgermination.1[1],survival.1[2:3],notgermination.1[2],survival.1[4:5],notgermination.1[3],survival.1[6])
  cumPersistence.1 = cumprod(cumPersistence)
  time = c(0,4,4,12,16,16,24,28,28,36)/12
  points(time,cumPersistence.1,type='o',pch=16,
         col=c(1,1,2,1,1,2,1,1,2,1))
  
  persistence.32[i] = cumPersistence.1[10]
  
  # round 2
  survival.2=(surv[1:4])
  germination.2=germ[4:5]
  notgermination.2 = 1-germination.2
  
  cumPersistence = c(1,survival.2[1],notgermination.2[1],survival.2[2:3],notgermination.2[2],survival.2[4])
  cumPersistence.2 = cumprod(cumPersistence)
  time = c(0,4,4,12,16,16,24)/12
  points(time,cumPersistence.2,type='o',pch=16,
         col=c(1,1,2,1,1,2,1,1,2,1),lty='dotted')
  
  # round 2
  survival.3=(surv[1:2])
  germination.3=germ[6]
  notgermination.3 = 1-germination.3
  
  cumPersistence = c(1,survival.3[1],notgermination.3[1],survival.2[2])
  cumPersistence.3 = cumprod(cumPersistence)
  time = c(0,4,4,12)/12
  points(time,cumPersistence.3,type='o',pch=16,
         col=c(1,1,2,1,1,2,1,1,2,1),lty='dashed')
  
  
  abline(v=c(6.5,10.5),lty='dotted')
  
  legend("topright",legend=siteNames[index],bty='n')
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}


plot(1:20,persistence.32)

plot(persistence.3,persistence.32,xlim=c(0,.5),ylim=c(0,.5))
abline(a=0,b=1)
