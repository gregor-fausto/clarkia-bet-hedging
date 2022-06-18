
rm(list=(ls(all=TRUE))) # if using in source(script)

# - Libraries ----
library(tidyverse)
library(parallel)
library(stringr)
library(bbmle)
library(MCMCvis)

# - Read in estimates ----

viabilityList <- readRDS(paste0("outputs/009_compareBayesianLikelihood/viabilityTrialsLikelihoodModels.rds"))
seedBagList <- readRDS(paste0("outputs/009_compareBayesianLikelihood/seedSamplesLikelihoodModels.rds"))

# - +Read in posterior for seed bag trials & viabiility trials ----

mcmcDirectory = "outputs/002_fitStatisticalModels/mcmcSamples/"
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedBagExperimentPosteriorSamples.RDS",mcmcSampleDirectory)]])
mcmcSamplesViability <- readRDS(mcmcSampleDirectory[[grep("viabilityPosteriorSamples.RDS",mcmcSampleDirectory)]])

# - Read in additional data ----

# - +Read in site names ----
siteAbiotic <- read.csv("data/siteAbiotic.csv",header=TRUE)

position<-siteAbiotic %>% 
  dplyr::ungroup() %>%
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

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

germinationList <- survivalList <- viab_gList <- viab_vList <- list()

for(i in 1:20){
  
  # observations
  model.tmp = seedBagList[[i]]
  viab.tmp = viabilityList[[i]]
  
  # viability
  # 1, 4, 6 are round 1
  # 2, 5 are round 2
  # 3 is round 3
  
  # get parameter estimates
  g_est = model.tmp@coef[1:6]
  theta_est = model.tmp@coef[7:18]
  nu_g_est = viab.tmp@coef[1:6]
  nu_v_est = viab.tmp@coef[7:12]
  
  index = c(1,2,3,4,5,6,1,2,3,4,1,2)
  round = c(1,1,1,1,1,1,2,2,2,2,3,3)
  age = c(1,1,2,2,3,3,1,1,2,2,1,1)
  survival = theta_est
  germination = g_est
  viability_g = nu_g_est
  viability_v = nu_v_est
  
  survivalList[[i]] <- survival
  germinationList[[i]] <- germination
  viab_gList[[i]] <- viability_g
  viab_vList[[i]] <- viability_v
}

# bind into matrix
germination.dat<-do.call(rbind,germinationList)
survival.dat<-do.call(rbind,survivalList)
viability_g.dat<-do.call(rbind,viab_gList)
viability_v.dat<-do.call(rbind,viab_vList)

likelihoodEstimates <- list(germination.dat,survival.dat,viability_g.dat,viability_v.dat)

# survival.dat[survival.dat > 1] <- 1
# colnames(survival.dat) = c("s1","s2","s3","s4","s5","s6","s1","s2","s3","s4","s1","s2")
# colnames(germination.dat) = c("g1","g2","g3","g1","g2","g1")
# colnames(viability.dat) = c("v1","v2","v3","v1","v2","v1")


# ---
# - Seed bag burial posteriors ----
# ---

# - +Germination: population- and year-level ----

mu_g=MCMCchains(mcmcSamples,params="mu_g")

# - ++index for age 1, 2, 3 ----
index.g1=grep(c(",1\\]|4\\]|6\\]"),colnames(mu_g))
index.g2=grep(",2\\]|5\\]",colnames(mu_g))
index.g3=grep(",3\\]",colnames(mu_g))

# - ++transform to probability scale and summarize ----
gamma1 = boot::inv.logit(mu_g[,index.g1])
gamma1.sum=apply(gamma1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

gamma2 = boot::inv.logit(mu_g[,index.g2])
gamma2.sum=apply(gamma2,2,quantile,probs=c(0.025,.25,.5,.75,.975))

gamma3 = boot::inv.logit(mu_g[,index.g3])
gamma3.sum=apply(gamma3,2,quantile,probs=c(0.025,.25,.5,.75,.975))

# - +Survival: population-level ----

mu_s=MCMCchains(mcmcSamples,params="mu_s")

index.s1=grep(",1\\]|,7\\]|,11\\]",colnames(mu_s))
index.s2=grep(",2\\]|,8\\]|,12\\]",colnames(mu_s))
index.s3=grep(",3\\]|,9\\]",colnames(mu_s))
index.s4=grep(",4\\]|,10\\]",colnames(mu_s))
index.s5=grep(",5\\]",colnames(mu_s))
index.s6=grep(",6\\]",colnames(mu_s))

sigma1 = boot::inv.logit(mu_s[,index.s1])
sigma1.sum=apply(sigma1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma2 = boot::inv.logit(mu_s[,index.s2])
sigma2.sum=apply(sigma2,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma3 = boot::inv.logit(mu_s[,index.s3])
sigma3.sum=apply(sigma3,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma4 = boot::inv.logit(mu_s[,index.s4])
sigma4.sum=apply(sigma4,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma5 = boot::inv.logit(mu_s[,index.s5])
sigma5.sum=apply(sigma5,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma6 = boot::inv.logit(mu_s[,index.s6])
sigma6.sum=apply(sigma6,2,quantile,probs=c(0.025,.25,.5,.75,.975))

gamma1.mode <- apply(gamma1,2,posterior.mode)
gamma2.mode <- apply(gamma2,2,posterior.mode)
gamma3.mode <- apply(gamma3,2,posterior.mode)
sigma1.mode <- apply(sigma1,2,posterior.mode)
sigma2.mode <- apply(sigma2,2,posterior.mode)
sigma3.mode <- apply(sigma3,2,posterior.mode)
sigma4.mode <- apply(sigma4,2,posterior.mode)
sigma5.mode <- apply(sigma5,2,posterior.mode)
sigma6.mode <- apply(sigma6,2,posterior.mode)

# - + Summarize estimates by medians ----

germinationBayesian <- cbind(gamma1.sum[3,1:20],gamma2.sum[3,1:20],gamma3.sum[3,],
                             gamma1.sum[3,21:40],gamma2.sum[3,21:40],
                             gamma1.sum[3,41:60])

survivalBayesian <- cbind(sigma1.sum[3,1:20],sigma2.sum[3,1:20],sigma3.sum[3,1:20],
                          sigma4.sum[3,1:20],sigma5.sum[3,1:20],sigma6.sum[3,1:20],
                          sigma1.sum[3,21:40],sigma2.sum[3,21:40],sigma3.sum[3,21:40],
                          sigma4.sum[3,21:40],
                          sigma1.sum[3,41:60],sigma2.sum[3,41:60])


survivalBayesian.mode <- cbind(sigma1.mode[1:20],sigma2.mode[1:20],sigma3.mode[1:20],
                               sigma4.mode[1:20],sigma5.mode[1:20],sigma6.mode[1:20],
                               sigma1.mode[21:40],sigma2.mode[21:40],sigma3.mode[21:40],
                               sigma4.mode[21:40],
                               sigma1.mode[41:60],sigma2.mode[41:60])


# - Compare estimates ----

par(mfrow=c(1,2))
plot(germination.dat,germinationBayesian,xlim=c(0,.6),ylim=c(0,.6))
abline(a=0,b=1)

plot(survival.dat,survivalBayesian,xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)

# - + Compare germination ----

par(mfrow=c(3,3),mar=c(3,3,2,2))

plot(germination.dat[,1],germinationBayesian[,1],xlim=c(0,.4),ylim=c(0,.4),
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

plot(germination.dat[,4],germinationBayesian[,4],xlim=c(0,.3),ylim=c(0,.3),
     main="Round 2: g1")
abline(a=0,b=1)
segments(x0=germination.dat[,4],y0=gamma1.sum[1,21:40],y1=gamma1.sum[5,21:40])

plot(germination.dat[,5],germinationBayesian[,5],xlim=c(0,.45),ylim=c(0,.45),
     main="Round 2: g2")
abline(a=0,b=1)
segments(x0=germination.dat[,5],y0=gamma2.sum[1,21:40],y1=gamma2.sum[5,21:40])

plot(NA,xlim=c(0,.6),ylim=c(0,.6),type='n',bty='n',xaxt='n',yaxt='n',
     main="")

plot(germination.dat[,6],germinationBayesian[,6],xlim=c(0,.7),ylim=c(0,.7),
     main="Round 3: g1")
abline(a=0,b=1)
segments(x0=germination.dat[,6],y0=gamma1.sum[1,41:60],y1=gamma1.sum[5,41:60])

plot(NA,xlim=c(0,.6),ylim=c(0,.6),type='n',bty='n',xaxt='n',yaxt='n',
     main="")

plot(NA,xlim=c(0,.6),ylim=c(0,.6),type='n',bty='n',xaxt='n',yaxt='n',
     main="")

# - + Compare survival ----

par(mfrow=c(2,3),mar=c(3,3,2,2))

plot(survival.dat[,1],survivalBayesian[,1],xlim=c(0,1),ylim=c(0,1),
     main="Round 1: s1")
abline(a=0,b=1)
segments(x0=survival.dat[,1],y0=sigma1.sum[1,1:20],y1=sigma1.sum[5,1:20])

plot(survival.dat[,2],survivalBayesian[,2],xlim=c(0,1),ylim=c(0,1),
     main="Round 1: s2")
abline(a=0,b=1)
segments(x0=survival.dat[,2],y0=sigma2.sum[1,1:20],y1=sigma2.sum[5,1:20])

plot(survival.dat[,3],survivalBayesian[,3],xlim=c(0,1),ylim=c(0,1),
     main="Round 1: s3")
abline(a=0,b=1)
segments(x0=survival.dat[,3],y0=sigma3.sum[1,1:20],y1=sigma3.sum[5,1:20])

plot(survival.dat[,4],survivalBayesian[,4],xlim=c(0,1),ylim=c(0,1),
     main="Round 1: s4")
abline(a=0,b=1)
segments(x0=survival.dat[,4],y0=sigma4.sum[1,1:20],y1=sigma4.sum[5,1:20])

plot(survival.dat[,5],survivalBayesian[,5],xlim=c(0,1),ylim=c(0,1),
     main="Round 1: s5")
abline(a=0,b=1)
segments(x0=survival.dat[,5],y0=sigma5.sum[1,1:20],y1=sigma5.sum[5,1:20])

plot(survival.dat[,6],survivalBayesian[,6],xlim=c(0,1),ylim=c(0,1),
     main="Round 1: s6")
abline(a=0,b=1)
segments(x0=survival.dat[,6],y0=sigma6.sum[1,1:20],y1=sigma6.sum[5,1:20])

par(mfrow=c(2,2),mar=c(3,3,2,2))

plot(survival.dat[,7],survivalBayesian[,7],xlim=c(0,1),ylim=c(0,1),
     main="Round 2: s1")
abline(a=0,b=1)
segments(x0=survival.dat[,7],y0=sigma1.sum[1,21:40],y1=sigma1.sum[5,21:40])

plot(survival.dat[,8],survivalBayesian[,8],xlim=c(0,1),ylim=c(0,1),
     main="Round 2: s2")
abline(a=0,b=1)
segments(x0=survival.dat[,8],y0=sigma2.sum[1,21:40],y1=sigma2.sum[5,21:40])

plot(survival.dat[,9],survivalBayesian[,9],xlim=c(0,1),ylim=c(0,1),
     main="Round 2: s3")
abline(a=0,b=1)
segments(x0=survival.dat[,9],y0=sigma3.sum[1,21:40],y1=sigma3.sum[5,21:40])

plot(survival.dat[,10],survivalBayesian[,10],xlim=c(0,1),ylim=c(0,1),
     main="Round 2: s4")
abline(a=0,b=1)
segments(x0=survival.dat[,10],y0=sigma4.sum[1,21:40],y1=sigma4.sum[5,21:40])


par(mfrow=c(1,2),mar=c(3,3,2,2))

plot(survival.dat[,11],survivalBayesian[,11],xlim=c(0,1),ylim=c(0,1),
     main="Round 3: s1")
abline(a=0,b=1)
segments(x0=survival.dat[,11],y0=sigma1.sum[1,41:60],y1=sigma1.sum[5,41:60])

plot(survival.dat[,12],survivalBayesian[,12],xlim=c(0,1),ylim=c(0,1),
     main="Round 3: s2")
abline(a=0,b=1)
segments(x0=survival.dat[,12],y0=sigma2.sum[1,41:60],y1=sigma2.sum[5,41:60])




# ---
# - Viability experiment posteriors ----
# ---

# - +Population estimates ----
# germination tests
mu_g<-MCMCchains(mcmcSamplesViability, params = "mu_g")
# viability assays
mu_v<-MCMCchains(mcmcSamplesViability, params = "mu_v")

# - ++index for age 1, 2, 3 ----
index.g1=grep(c(",1\\]|4\\]|6\\]"),colnames(mu_g))
index.g2=grep(",2\\]|5\\]",colnames(mu_g))
index.g3=grep(",3\\]",colnames(mu_g))

# - ++transform to probability scale and summarize ----
mu_g = boot::inv.logit(mu_g)
mu_g.sum=apply(mu_g,2,quantile,probs=c(0.025,.25,.5,.75,.975))

mu_v = boot::inv.logit(mu_v)
mu_v.sum=apply(mu_v,2,quantile,probs=c(0.025,.25,.5,.75,.975))


viability_gBayesian <- cbind(mu_g.sum[3,1:20],mu_g.sum[3,21:40],
                             mu_g.sum[3,41:60],mu_g.sum[3,61:80],
                             mu_g.sum[3,81:100],mu_g.sum[3,101:120])

viability_vBayesian <- cbind(mu_v.sum[3,1:20],mu_v.sum[3,21:40],
                             mu_v.sum[3,41:60],mu_v.sum[3,61:80],
                             mu_v.sum[3,81:100],mu_v.sum[3,101:120])


# - Compare viability estimates ----

par(mfrow=c(1,2))
plot(viability_g.dat,viability_gBayesian,xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)

plot(viability_v.dat,viability_vBayesian,xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)

# - + Compare germination ----

par(mfrow=c(2,3),mar=c(3,3,2,2))

plot(viability_g.dat[,1],viability_gBayesian[,1],xlim=c(0,1),ylim=c(0,1),
     main="Round 1: age 1")
abline(a=0,b=1)
segments(x0=viability_g.dat[,1],y0=mu_g.sum[1,1:20],y1=mu_g.sum[5,1:20])

plot(viability_g.dat[,2],viability_gBayesian[,2],xlim=c(0,1),ylim=c(0,1),
     main="Round 2: age 1")
abline(a=0,b=1)
segments(x0=viability_g.dat[,2],y0=mu_g.sum[1,21:40],y1=mu_g.sum[5,21:40])

plot(viability_g.dat[,3],viability_gBayesian[,3],xlim=c(0,1),ylim=c(0,1),
     main="Round 3: age 1")
abline(a=0,b=1)
segments(x0=viability_g.dat[,3],y0=mu_g.sum[1,41:60],y1=mu_g.sum[5,41:60])


plot(viability_g.dat[,4],viability_gBayesian[,4],xlim=c(0,1),ylim=c(0,1),
     main="Round 1: age 2")
abline(a=0,b=1)
segments(x0=viability_g.dat[,4],y0=mu_g.sum[1,61:80],y1=mu_g.sum[5,61:80])

plot(viability_g.dat[,5],viability_gBayesian[,5],xlim=c(0,1),ylim=c(0,1),
     main="Round 2: age 2")
abline(a=0,b=1)
segments(x0=viability_g.dat[,5],y0=mu_g.sum[1,81:100],y1=mu_g.sum[5,81:100])

plot(viability_g.dat[,6],viability_gBayesian[,6],xlim=c(0,1),ylim=c(0,1),
     main="Round 1: age 3")
abline(a=0,b=1)
segments(x0=viability_g.dat[,6],y0=mu_g.sum[1,101:120],y1=mu_g.sum[5,101:120])


# - + Compare viability ----

par(mfrow=c(2,3),mar=c(3,3,2,2))

plot(viability_v.dat[,1],viability_vBayesian[,1],xlim=c(0,1),ylim=c(0,1),
     main="Round 1: age 1")
abline(a=0,b=1)
segments(x0=viability_v.dat[,1],y0=mu_v.sum[1,1:20],y1=mu_v.sum[5,1:20])

plot(viability_v.dat[,2],viability_vBayesian[,2],xlim=c(0,1),ylim=c(0,1),
     main="Round 2: age 1")
abline(a=0,b=1)
segments(x0=viability_v.dat[,2],y0=mu_v.sum[1,21:40],y1=mu_v.sum[5,21:40])

plot(viability_v.dat[,3],viability_vBayesian[,3],xlim=c(0,1),ylim=c(0,1),
     main="Round 3: age 1")
abline(a=0,b=1)
segments(x0=viability_v.dat[,3],y0=mu_v.sum[1,41:60],y1=mu_v.sum[5,41:60])


plot(viability_v.dat[,4],viability_vBayesian[,4],xlim=c(0,1),ylim=c(0,1),
     main="Round 1: age 2")
abline(a=0,b=1)
segments(x0=viability_v.dat[,4],y0=mu_v.sum[1,61:80],y1=mu_v.sum[5,61:80])

plot(viability_v.dat[,5],viability_vBayesian[,5],xlim=c(0,1),ylim=c(0,1),
     main="Round 2: age 2")
abline(a=0,b=1)
segments(x0=viability_v.dat[,5],y0=mu_v.sum[1,81:100],y1=mu_v.sum[5,81:100])

plot(viability_v.dat[,6],viability_vBayesian[,6],xlim=c(0,1),ylim=c(0,1),
     main="Round 1: age 3")
abline(a=0,b=1)
segments(x0=viability_v.dat[,6],y0=mu_v.sum[1,101:120],y1=mu_v.sum[5,101:120])



# - Look at germination data ----

dataList <- readRDS(file=paste0("outputs/009_compareBayesianLikelihood/seedSamplesLikelihoodData.rds"))

observationTally <- observationNumber <- observationVar <- list()

for(i in 1:20){
  
  # observations
  obs.number = dataList[[i]][[6]]
  tmp.tally <-  table(obs.number)
  
  # sample size: start
  obs.samplesize = dataList[[i]][[3]]
  tmp.samplesize = c()
  for(j in 1:6){tmp.samplesize[j] = sum(obs.samplesize[obs.number==j])}

  # sample: variability in outcome
  obs.samplesize = dataList[[i]][[4]]
  tmp.var = c()
  for(j in 1:6){tmp.var[j] = var(obs.samplesize[obs.number==j])}
  
  
  observationTally[[i]] <- tmp.tally
  observationNumber[[i]] <- tmp.samplesize
  observationVar[[i]] <- tmp.var
  
}

tally.dat <-do.call(rbind,observationTally)
sample.dat <-do.call(rbind,observationNumber)
var.dat <-do.call(rbind,observationVar)

# mle - bayesian estimate
germination.delta=(germination.dat-germinationBayesian)

par(mfrow=c(1,3))

# difference in estimates vs. number of observations
plot(tally.dat,germination.delta)
abline(h=0)

# more total seed numbers decreases difference between estimates 
plot(sample.dat,germination.delta)
abline(h=0)

# but main effect of pooling is due to 
# greater variability among bags in the germination estimate 
# increases the difference between estimate from frequentist vs. bayesian estimate
plot(var.dat,germination.delta,ylab="MLE-Posterior median")
abline(h=0)

par(mfrow=c(2,3))
for(i in 1:6){
plot(var.dat[,i],germination.delta[,i],ylab="MLE-Posterior median")
abline(h=0)
}

df=data.frame(position[,1:2],variance=var.dat,germ.delta=germination.delta)


# - Look at survival data ----


dataList <- readRDS(file=paste0("outputs/009_compareBayesianLikelihood/seedSamplesLikelihoodData.rds"))

observationTally <- observationNumber <- observationVar <- list()

for(i in 1:20){
  
  # observations
  obs.number = dataList[[i]][[5]]
  tmp.tally <-  table(obs.number)
  
  # sample size: start
  obs.samplesize = dataList[[i]][[1]]
  tmp.samplesize = c()
  for(j in 1:12){tmp.samplesize[j] = sum(obs.samplesize[obs.number==j])}
  
  # sample: variability in outcome
  obs.samplesize = dataList[[i]][[2]]
  tmp.var = c()
  for(j in 1:12){tmp.var[j] = var(obs.samplesize[obs.number==j])}
  
  
  observationTally[[i]] <- tmp.tally
  observationNumber[[i]] <- tmp.samplesize
  observationVar[[i]] <- tmp.var
  
}

tally.dat <-do.call(rbind,observationTally)
sample.dat <-do.call(rbind,observationNumber)
var.dat <-do.call(rbind,observationVar)

# mle - bayesian estimate
survival.delta=(survival.dat-survivalBayesian)

par(mfrow=c(1,3))

# difference in estimates vs. number of observations
plot(tally.dat[,c(1,7,11)],survival.delta[,c(1,7,11)])
abline(h=0)

plot(sample.dat[,c(1,7,11)],survival.delta[,c(1,7,11)])
abline(h=0)

#  main effect of pooling is due to 
# greater variability among bags in the estimates
# increases the difference between estimate from frequentist vs. bayesian estimate
# but note that in the survival data this is mostly obvious in the 1st observation
plot(var.dat[,c(1,7,11)],survival.delta[,c(1,7,11)],ylab="MLE-Posterior median")
abline(h=0)

# the effect disappears when looking at everthing in aggregate
plot(var.dat,survival.delta,ylab="MLE-Posterior median")
points(var.dat[,c(1,7,11)],survival.delta[,c(1,7,11)],pch=16,col='red')
abline(h=0)

# this is because the bias is compensating
# low estimates in s2 become high estimates the following time point, and then continue to oscillate
par(mfrow=c(4,5))
#plot(1:6,survival.delta[1,1:6],type='n',ylim=c(-.4,.3))
for(i in 1:20){
  plot(1:6,survival.delta[i,1:6],type='b')
  abline(h=0,lty='dotted')
}

par(mfrow=c(1,2))
plot(1:4,survival.delta[1,7:10],type='n',ylim=c(-.4,.3))
abline(h=0,lty='dotted')
for(i in 1:20){
  points(1:4,survival.delta[i,7:10],type='b')
}

plot(1:2,survival.delta[1,11:12],type='n',ylim=c(-.4,.3))
abline(h=0,lty='dotted')
for(i in 1:20){
  points(1:2,survival.delta[i,11:12],type='b')
}

