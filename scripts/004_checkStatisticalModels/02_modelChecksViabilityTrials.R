####
####
# Script for model checks for viability
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

mcmcDirectory = "outputs/002_fitStatisticalModels/mcmcSamples/"
modelDataDirectory = "outputs/002_fitStatisticalModels/data/"
outputDirectory = "outputs/004_checkStatisticalModels/"

# - Libraries ----
library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)

# - Read in what's needed for plotting ----

# - +Read in data ----
modelDataDirectory <- paste0(modelDataDirectory,list.files(modelDataDirectory))
data <- readRDS(modelDataDirectory[[grep("viabilityData.RDS",modelDataDirectory)]])

# - +Read in MCMC samples ----
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("viabilityPosteriorSamples.RDS",mcmcSampleDirectory)]])

# - +Read in site names ----
siteAbiotic <- read.csv("data/siteAbioticData.csv",header=TRUE)
# create variable that arranges populations by easting
siteIndex <- order(siteAbiotic$easting,decreasing=FALSE)
siteNames = unique(siteAbiotic$site)[siteIndex]


# - Function to plot p-vals ----

f.plot = function(diagnostic.test.mat,years){
  
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(min(years)-.5,max(years)+.5),
       xlab="",ylab="",
       main=NULL,type='n',
       xaxt='n',yaxt='n',frame=FALSE)
  
  start = years[1]
  offset = seq(0,length(years),by=1)
  
  for(j in 1:length(years)){
    
    points(x=rep(start+offset[j],20)+rnorm(20,0,.05),
           pch=21,cex=.8,
           bg=ifelse(diagnostic.test.mat[,j]>0.95|
                       diagnostic.test.mat[,j]<0.05,'orange','gray95'),
           y=diagnostic.test.mat[,j])
    
  }
  
  abline(h=.05,lty='dotted')
  abline(h=.95,lty='dotted')
  box(which="plot",bty="l",col='black')
  
  axis(1, seq(min(years),max(years),by=1),
       labels = seq(min(years),max(years),by=1), 
       las = 1,
       col = NA, col.ticks = 1, cex.axis = 1)
  axis(2, seq(0,1,by=.2),
       seq(0,1,by=.2), las = 1,
       col = NA, col.ticks = 1, cex.axis = 1)
}

# ---
# Germination ----
# ---

# ---
# *Posterior predictive checks ----
# ---

# extract chi-2 values
chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.germCount.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.germCount.sim"))

n.iter = dim(chi2.obs)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=6)
for(j in 1:20){
  for(i in 1:6){
    tmp.chi2.obs=chi2.obs[,data$siteViab==j&data$indexViab==i]
    tmp.chi2.sim=chi2.sim[,data$siteViab==j&data$indexViab==i]
    fit.obs=apply(tmp.chi2.obs,1,sum)
    fit.sim=apply(tmp.chi2.sim,1,sum)
    p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
    p.pop[,i] = p.chi2.calc
  }
  p.chi.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
}
p.chi.mat=do.call(rbind,p.chi.list)

f=function(y.sim=chains,y.obs=data,n.obs=data2,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=6)
  for(j in 1:20){
    for(i in 1:6){
      index=data$siteViab==j&data$indexViab==i
      tmp.obs=y.obs[index]/n.obs[index]
      tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
      test.obs=model.fun(tmp.obs)
      test.sim=apply(tmp.sim,1,model.fun)
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      p.pop[,i] = p.test.calc
    }
    p.test.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
  }
  p.test.mat=do.call(rbind,p.test.list)
  return(p.test.mat)
}

sims=MCMCchains(mcmcSamples,params=c("germCount_sim"))
df=data$germCount
df2=data$germStart

germ.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)

# convert lists to matrix
f.convert = function(test.list){
  out.list <- list()
  for(i in 1:20){
    tmpYear = data$indexViab[data$siteViab==i]
    
    out.list[[i]] <- ifelse((1:6 %in% tmpYear), test.list[[i]],
                            ifelse(!(1:6 %in% tmpYear), NA, 0))
  }
  out.mat <- do.call(rbind,out.list)
  return(out.mat)
}

# return matrix objects
germ.mean.mat<-germ.mean
p.chi.mat<-f.convert(p.chi.list)


pdf(file=paste0(outputDirectory,"pvals-germinationSeedBagExperiment.pdf"),height=3,width=6)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

# germination g1 in 2006, 2007, and 2008 
# correspond to rounds 1, 2, 3
f.plot(diagnostic.test.mat = germ.mean.mat[,1:3],years=2006:2008)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat[,1:3],years=2006:2008)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

# germination g1 in  2007, and 2008 
# correspond to rounds 2
f.plot(diagnostic.test.mat = germ.mean.mat[,4:5],years=2007:2008)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat[,4:5],years=2007:2008)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)


dev.off()


# ---
# Graphical checks ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "germCount_sim")
n.iter=dim(y.sim)[1]

f = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density(x, from= x.min, to = x.max)
  return(dres)
}

pdf(file=paste0(outputDirectory,"ppc-viabilityTrialsGermination.pdf"),height=6,width=6)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,.25,.25) + 0.1)
for(i in 1:20){
  
  tmp.i = siteIndex[i]
  
  # for each site get index for year 1 germination
  index=data$siteViab==tmp.i&data$indexViab==1
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$germStart[index], FUN = '/')
  
  sample.index = sample(1:n.iter,50)
  tmp=tmp[sample.index,]

  for(j in 1:50){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,
          col=ifelse(p.chi.mat[i,1]>.95|p.chi.mat[i,1]<.05,
                     'orange','gray75'))
  }

  dres = f(data$germCount[index]/data$germStart[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  mtext(siteNames[i], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination in lab trial\n (1 year after seed production, round 1)", 
      side = 1, outer = TRUE, line = 3.5)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,.25,.25) + 0.1)
for(i in 1:20){
  
  tmp.i = siteIndex[i]
  
  index=data$siteViab==tmp.i&data$indexViab==2
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$germStart[index], FUN = '/')
  
  sample.index = sample(1:n.iter,50)
  tmp=tmp[sample.index,]
  
  for(j in 1:50){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,
          col=ifelse(p.chi.mat[i,2]>.95|p.chi.mat[i,2]<.05,
                     'orange','gray75'))
  }
  
  dres = f(data$germCount[index]/data$germStart[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  mtext(siteNames[i], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination in lab trial\n (2 years after seed production, round 1)", 
      side = 1, outer = TRUE, line = 3.5)
mtext("Density", side = 2, outer = TRUE, line = 2.2)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,.25,.25) + 0.1)
for(i in 1:20){
  
  tmp.i = siteIndex[i]
  
  # for each site get index for year 1 germination
  index=data$siteViab==tmp.i&data$indexViab==4
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$germStart[index], FUN = '/')
  
  sample.index = sample(1:n.iter,50)
  tmp=tmp[sample.index,]
  
  for(j in 1:50){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,
          col=ifelse(p.chi.mat[i,4]>.95|p.chi.mat[i,4]<.05,
                     'orange','gray75'))
  }
  
  dres = f(data$germCount[index]/data$germStart[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  mtext(siteNames[i], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination in lab trial\n (1 year after seed production, round 2)", 
      side = 1, outer = TRUE, line = 3.5)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,.25,.25) + 0.1)
for(i in 1:20){
  
  tmp.i = siteIndex[i]
  
  index=data$siteViab==tmp.i&data$indexViab==5
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$germStart[index], FUN = '/')
  
  sample.index = sample(1:n.iter,50)
  tmp=tmp[sample.index,]
  
  for(j in 1:50){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,
          col=ifelse(p.chi.mat[i,5]>.95|p.chi.mat[i,5]<.05,
                     'orange','gray75'))
  }
  
  dres = f(data$germCount[index]/data$germStart[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  mtext(siteNames[i], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination in lab trial\n (2 years after seed production, round 2)", 
      side = 1, outer = TRUE, line = 3.5)
mtext("Density", side = 2, outer = TRUE, line = 2.2)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,.25,.25) + 0.1)
for(i in 1:20){
  
  tmp.i = siteIndex[i]
  
  index=data$siteViab==tmp.i&data$indexViab==6
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$germStart[index], FUN = '/')
  
  sample.index = sample(1:n.iter,50)
  tmp=tmp[sample.index,]
  
  for(j in 1:50){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,
          col=ifelse(p.chi.mat[i,6]>.95|p.chi.mat[i,6]<.05,
                     'orange','gray75'))
  }
  
  dres = f(data$germCount[index]/data$germStart[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  mtext(siteNames[i], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination in lab trial\n (1 years after seed production, round 3)", 
      side = 1, outer = TRUE, line = 3.5)
mtext("Density", side = 2, outer = TRUE, line = 2.2)


dev.off()

# ---
# Viable, staining seeds count ----
# ---

# ---
# *Posterior predictive checks ----
# ---

# get posterior for chi-2
chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.viabStain.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.viabStain.sim"))

# get length of chains
n.iter = dim(chi2.obs)[1]

# construct list to hold p-value for chi-square
p.chi.list = list()

# put together the p-values by age of seeds
p.pop = matrix(NA,nrow=n.iter,ncol=6)
for(j in 1:20){
    
  for(i in 1:6){
    index = data$siteViab_v==j&data$indexViab_v==i
    tmp.chi2.obs=chi2.obs[,index]
    tmp.chi2.sim=chi2.sim[,index]
    fit.obs=apply(tmp.chi2.obs,1,sum)
    fit.sim=apply(tmp.chi2.sim,1,sum)
    p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
    p.pop[,i] = p.chi2.calc
  }
  p.chi.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
}
p.chi.mat=do.call(rbind,p.chi.list)

# write a function to calculate the p-value for different functions
f=function(y.sim=chains,y.obs=data,n.obs=data2,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=6)
  for(j in 1:20){
    
    for(i in 1:6){
      index = data$siteViab_v==j&data$indexViab_v==i
      tmp.obs=y.obs[index]/n.obs[index]
      tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
      test.obs=model.fun(tmp.obs)
      test.sim=apply(tmp.sim,1,model.fun)
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      p.pop[,i] = p.test.calc
    }
    p.test.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
  }
  p.test.mat=do.call(rbind,p.test.list)
  return(p.test.mat)
}

sims=MCMCchains(mcmcSamples,params=c("viabStain_sim"))
df=data$viabStain
df2=data$viabStart

seeds.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)

# convert lists to matrix
f.convert = function(test.list){
  out.list <- list()
  for(i in 1:20){
    tmpYear = data$indexViab_v[data$siteViab_v==i]
    
    out.list[[i]] <- ifelse((1:6 %in% tmpYear), test.list[[i]],
                            ifelse(!(1:6 %in% tmpYear), NA, 0))
  }
  out.mat <- do.call(rbind,out.list)
  return(out.mat)
}

# return matrix objects
seeds.mean.mat<-seeds.mean
p.chi.mat<-f.convert(p.chi.list)

pdf(file=paste0(outputDirectory,"pvals-viabilityTrialsViability.pdf"),height=3,width=6)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

# correspond to rounds 1
f.plot(diagnostic.test.mat = seeds.mean.mat[,1:3],years=2006:2008)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat[,1:3],years=2006:2008)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

# correspond to rounds 2
f.plot(diagnostic.test.mat = seeds.mean.mat[,4:5],years=2007:2008)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat[,4:5],years=2007:2008)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)


dev.off()

# ---
# Graphical checks ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "viabStain_sim")
n.iter=dim(y.sim)[1]

f = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density(x, from= x.min, to = x.max)
  return(dres)
}

pdf(file=paste0(outputDirectory,"ppc-viabilityTrialsViability.pdf"),height=6,width=6)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,.25,.25) + 0.1)
for(i in 1:20){
  
  tmp.i = siteIndex[i]
  
  # for each site get index for year 1 germination
  index=data$siteViab_v==tmp.i&data$indexViab_v==1
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$viabStart[index], FUN = '/')
  
  sample.index = sample(1:n.iter,50)
  tmp=tmp[sample.index,]
  
  for(j in 1:50){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,
          col=ifelse(p.chi.mat[i,1]>.95|p.chi.mat[i,1]<.05,
                     'orange','gray75'))
  }
  
  dres = f(data$viabStain[index]/data$viabStart[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  mtext(siteNames[i], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of viability, conditional on not germinating, in lab trial\n (1 year after seed production, round 1)", 
      side = 1, outer = TRUE, line = 3.5,cex=.8)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,.25,.25) + 0.1)
for(i in 1:20){
  
  tmp.i = siteIndex[i]
  
  index=data$siteViab_v==tmp.i&data$indexViab_v==2
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$viabStart[index], FUN = '/')
  
  sample.index = sample(1:n.iter,50)
  tmp=tmp[sample.index,]
  
  for(j in 1:50){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,
          col=ifelse(p.chi.mat[i,2]>.95|p.chi.mat[i,2]<.05,
                     'orange','gray75'))
  }
  
  dres = f(data$viabStain[index]/data$viabStart[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  mtext(siteNames[i], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of viability, conditional on not germinating, in lab trial\n (2 years after seed production, round 1)", 
      side = 1, outer = TRUE, line = 3.5,cex=.8)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,.25,.25) + 0.1)
for(i in 1:20){
  
  tmp.i = siteIndex[i]
  
  # for each site get index for year 1 germination
  index=data$siteViab_v==tmp.i&data$indexViab_v==4
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$viabStart[index], FUN = '/')
  
  sample.index = sample(1:n.iter,50)
  tmp=tmp[sample.index,]
  
  for(j in 1:50){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,
          col=ifelse(p.chi.mat[i,4]>.95|p.chi.mat[i,4]<.05,
                     'orange','gray75'))
  }
  
  dres = f(data$viabStain[index]/data$viabStart[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  mtext(siteNames[i], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of viability, conditional on not germinating, in lab trial\n (1 year after seed production, round 2)", 
      side = 1, outer = TRUE, line = 3.5,cex=.8)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,.25,.25) + 0.1)
for(i in 1:20){
  
  tmp.i = siteIndex[i]
  
  index=data$siteViab_v==tmp.i&data$indexViab_v==5
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$viabStart[index], FUN = '/')
  
  sample.index = sample(1:n.iter,50)
  tmp=tmp[sample.index,]
  
  for(j in 1:50){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,
          col=ifelse(p.chi.mat[i,5]>.95|p.chi.mat[i,5]<.05,
                     'orange','gray75'))
  }
  
  dres = f(data$viabStain[index]/data$viabStart[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  mtext(siteNames[i], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of viability, conditional on not germinating, in lab trial\n (2 years after seed production, round 2)", 
      side = 1, outer = TRUE, line = 3.5,cex=.8)
mtext("Density", side = 2, outer = TRUE, line = 2.2)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,.25,.25) + 0.1)
for(i in 1:20){
  
  tmp.i = siteIndex[i]
  
  index=data$siteViab_v==tmp.i&data$indexViab_v==6
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$viabStart[index], FUN = '/')
  
  sample.index = sample(1:n.iter,50)
  tmp=tmp[sample.index,]
  
  for(j in 1:50){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,
          col=ifelse(p.chi.mat[i,6]>.95|p.chi.mat[i,6]<.05,
                     'orange','gray75'))
  }
  
  dres = f(data$viabStain[index]/data$viabStart[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  mtext(siteNames[i], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of viability, conditional on not germinating, in lab trial\n (1 year after seed production, round 3)", 
      side = 1, outer = TRUE, line = 3.5,cex=.8)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

dev.off()