
rm(list=(ls(all=TRUE))) # if using in source(script)

# - Libraries ----
library(tidyverse)
library(parallel)
library(stringr)
library(bbmle)

# - Read in estimates ----

viabilityList <- readRDS(paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","viabilityEstimates.rds"))
seedBagList <- readRDS(paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","seedSamplesLikelihoodModels.rds"))

# - +Read seed bag trials & posterior ----

fileDirectory = "/Users/Gregor/Dropbox/dataLibrary/seed-banks-reanalysis-2/"

mcmcSampleDirectory <- paste0(fileDirectory,list.files(fileDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedSamplesNonparametric.rds",mcmcSampleDirectory)]])

# - +Read viability trials & posterior ----

mcmcSamplesViability <- readRDS(mcmcSampleDirectory[[grep("viabilityTrialSamples.rds",mcmcSampleDirectory)]])

# - Read in additional data ----

climate <- readRDS("~/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
position<-climate %>% 
  dplyr::ungroup() %>%
  dplyr::select(site,easting,intenseDemography) %>%
  dplyr::mutate(easting=easting/1000)

position = position %>% dplyr::filter(intenseDemography==1) %>% unique

orderedIndex = order(position$easting,decreasing=FALSE)
siteNames = unique(position$site)

# - Functions ----

posterior.mode = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]
  return(modeParam)
}


# - Calculate rates ----

germinationList <- survivalList <- viabList <- list()

for(i in 1:20){
  
  # observations
  model.tmp = seedBagList[i][[1]]
  viab.tmp = viabilityList[[i]]
  
  # viability
  # 1, 4, 6 are round 1
  # 2, 5 are round 2
  # 3 is round 3
  
  # get parameter estimates
  g_est = model.tmp@coef[1:6]
  theta_est = model.tmp@coef[7:18]
  nu_est = viab.tmp[c(1,4,6,2,5,3)]
  
  index = c(1,2,3,4,5,6,1,2,3,4,1,2)
  round = c(1,1,1,1,1,1,2,2,2,2,3,3)
  age = c(1,1,2,2,3,3,1,1,2,2,1,1)
  survival = theta_est
  germination = g_est
  viability = nu_est
  
  survivalList[[i]] <- survival
  germinationList[[i]] <- germination
  viabList[[i]] <- viability
}

# bind into matrix
germination.dat<-do.call(rbind,germinationList)
survival.dat<-do.call(rbind,survivalList)
viability.dat<-do.call(rbind,viabList)

likelihoodEstimates <- list(germination.dat,survival.dat,viability.dat)

survival.dat[survival.dat > 1] <- 1
colnames(survival.dat) = c("s1","s2","s3","s4","s5","s6","s1","s2","s3","s4","s1","s2")
colnames(germination.dat) = c("g1","g2","g3","g1","g2","g1")
colnames(viability.dat) = c("v1","v2","v3","v1","v2","v1")

# calculate persistence, on average
survival.dat = cbind(apply(survival.dat[,c(1,7,11)],1,mean),
                     apply(survival.dat[,c(2,8,12)],1,mean),
                     apply(survival.dat[,c(3,9)],1,mean),
                     apply(survival.dat[,c(4,10)],1,mean),
                     survival.dat[,5],
                     survival.dat[,6])

# calculate germination, on average
germination.dat = cbind(apply(germination.dat[,c(1,4,6)],1,mean),
                        apply(germination.dat[,c(2,5)],1,mean),
                        germination.dat[,3])


# ---
# - Seed bag burial posteriors ----
# ---

# - +Germination: population- and year-level ----

mu0_g=MCMCchains(mcmcSamples,params="mu0_g")

# - ++index for age 1, 2, 3 ----
index.g1=grep(c(",1\\]"),colnames(mu0_g))
index.g2=grep(",2\\]",colnames(mu0_g))
index.g3=grep(",3\\]",colnames(mu0_g))

# - ++transform to probability scale and summarize ----
gamma1 = boot::inv.logit(mu0_g[,index.g1])
gamma1.sum=apply(gamma1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

gamma2 = boot::inv.logit(mu0_g[,index.g2])
gamma2.sum=apply(gamma2,2,quantile,probs=c(0.025,.25,.5,.75,.975))

gamma3 = boot::inv.logit(mu0_g[,index.g3])
gamma3.sum=apply(gamma3,2,quantile,probs=c(0.025,.25,.5,.75,.975))

# - +Survival: population-level ----

mu0_s=MCMCchains(mcmcSamples,params="mu0_s")

index.s1=grep(",1\\]",colnames(mu0_s))
index.s2=grep(",2\\]",colnames(mu0_s))
index.s3=grep(",3\\]",colnames(mu0_s))
index.s4=grep(",4\\]",colnames(mu0_s))
index.s5=grep(",5\\]",colnames(mu0_s))
index.s6=grep(",6\\]",colnames(mu0_s))

sigma1 = boot::inv.logit(mu0_s[,index.s1])
sigma1.sum=apply(sigma1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma2 = boot::inv.logit(mu0_s[,index.s2])
sigma2.sum=apply(sigma2,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma3 = boot::inv.logit(mu0_s[,index.s3])
sigma3.sum=apply(sigma3,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma4 = boot::inv.logit(mu0_s[,index.s4])
sigma4.sum=apply(sigma4,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma5 = boot::inv.logit(mu0_s[,index.s5])
sigma5.sum=apply(sigma5,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma6 = boot::inv.logit(mu0_s[,index.s6])
sigma6.sum=apply(sigma6,2,quantile,probs=c(0.025,.25,.5,.75,.975))


# - + Summarize estimates by medians ----

germinationBayesian <- cbind(gamma1.sum[3,1:20],gamma2.sum[3,1:20],gamma3.sum[3,1:20])

survivalBayesian <- cbind(sigma1.sum[3,1:20],sigma2.sum[3,1:20],sigma3.sum[3,1:20],
                          sigma4.sum[3,1:20],sigma5.sum[3,1:20],sigma6.sum[3,1:20])

# - Compare estimates ----

par(mfrow=c(1,2))
plot(germination.dat,germinationBayesian,xlim=c(0,.6),ylim=c(0,.6))
abline(a=0,b=1)

plot(survival.dat,survivalBayesian,xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)

# - + Compare germination ----

par(mfrow=c(1,3),mar=c(3,3,2,2))

plot(germination.dat[,1],germinationBayesian[,1],xlim=c(0,.6),ylim=c(0,.6),
     main="Round 1: g1")
abline(a=0,b=1)
segments(x0=germination.dat[,1],y0=gamma1.sum[1,1:20],y1=gamma1.sum[5,1:20])

plot(germination.dat[,2],germinationBayesian[,2],xlim=c(0,.6),ylim=c(0,.6),
     main="Round 1: g2")
abline(a=0,b=1)
segments(x0=germination.dat[,2],y0=gamma2.sum[1,1:20],y1=gamma2.sum[5,1:20])

plot(germination.dat[,3],germinationBayesian[,3],xlim=c(0,.8),ylim=c(0,.8),
     main="Round 1: g3")
abline(a=0,b=1)
segments(x0=germination.dat[,3],y0=gamma3.sum[1,1:20],y1=gamma3.sum[5,1:20])

# - Calculate probability of persistence ----

# - +persistence & viability: round 1 ----

viableAfterThreeYears.freq <-survival.dat[,1] * (1-germination.dat[,1]) * survival.dat[,2]  * survival.dat[,3] * (1-germination.dat[2]) * survival.dat[,4]  * survival.dat[,5] * (1-germination.dat[3]) * survival.dat[,6]
viableAfterThreeYears.bayes <- s1[,1:20] * (1-g1[,1:20]) * s2[,1:20]  * s3[,1:20] * (1-g2[,1:20]) * s4[,1:20]  * s5[,1:20] * (1-g3[,1:20]) * s6[,1:20]
viableAfterThreeYears.bayes.summary <- apply(viableAfterThreeYears.bayes,2,quantile,c(.025,.5,.975))

par(mfrow=c(1,1),mar=c(4,4,4,4))

plot(viableAfterThreeYears.freq,viableAfterThreeYears.bayes.summary[2,],xlim=c(0,.3),ylim=c(0,.3),
     main="Persistent & viable after 3 years",xlab="Frequentist estimate",ylab="Bayesian estimate")
abline(a=0,b=1)
segments(x0=viableAfterThreeYears.freq,y0=viableAfterThreeYears.bayes.summary[1,],y1=viableAfterThreeYears.bayes.summary[3,])

plot(position$easting, viableAfterThreeYears.freq,ylim=c(0,.3),col='black',pch=16,type='n')
points(position$easting, viableAfterThreeYears.bayes.summary[2,],col='red',pch=16)
segments(x0=position$easting,y0=viableAfterThreeYears.bayes.summary[1,],y1=viableAfterThreeYears.bayes.summary[3,],col='red')
points(position$easting, viableAfterThreeYears.freq,col='black',pch=16)
