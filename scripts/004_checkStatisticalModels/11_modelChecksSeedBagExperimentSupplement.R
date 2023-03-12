####
####
# Script for model checks for seed bag experiment
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

mcmcDirectory = "outputs/002_fitStatisticalModels/mcmcSamples/"
modelDataDirectory = "outputs/002_fitStatisticalModels/data/"
outputDirectory = "outputs/004_checkStatisticalModels/01_modelChecksSupplement/"

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
data <- readRDS(modelDataDirectory[[grep("seedData.RDS",modelDataDirectory)]])

# - +Read in MCMC samples ----
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedBagExperimentPosteriorSamples.RDS",mcmcSampleDirectory)]])

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

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sim"))

# calculations are rowwise
fit.obs=apply(chi2.obs,1,sum)
fit.sim=apply(chi2.sim,1,sum)
p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
mean(p.chi2.calc)

# gIndex = list(c(1,4,6),c(2,5),c(3))
p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=6)
for(j in 1:20){
  for(i in 1:6){
    tmp.chi2.obs=chi2.obs[,data$siteGermination==j&data$germinationIndex %in% i]
    tmp.chi2.sim=chi2.sim[,data$siteGermination==j&data$germinationIndex %in% i]
    fit.obs=apply(tmp.chi2.obs,1,sum)
    fit.sim=apply(tmp.chi2.sim,1,sum)
    #   exclude=fit.obs!=fit.sim
    #   p.chi2.calc=ifelse(fit.sim[exclude]-fit.obs[exclude]>=0,1,0)
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
      index=data$siteGermination==j&data$germinationIndex %in% i
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

sims=MCMCchains(mcmcSamples,params=c("seedlingJan_sim"))
df=data$seedlingJan
df2=data$totalJan

germ.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)

# convert lists to matrix
f.convert = function(test.list){
  out.list <- list()
  for(i in 1:20){
    tmpYear = data$germinationIndex[data$siteGermination==i]
    
    out.list[[i]] <- ifelse((1:6 %in% tmpYear), test.list[[i]],
                            ifelse(!(1:6 %in% tmpYear), NA, 0))
  }
  out.mat <- do.call(rbind,out.list)
  return(out.mat)
}

# return matrix objects
germ.mean.mat<-germ.mean
p.chi.mat<-f.convert(p.chi.list)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(3)

pdf(file=paste0(outputDirectory,"pvals-germinationSeedBagExperiment.pdf"),height=3,width=6)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

# germination g1 in 2006, 2007, and 2008 
# correspond to germinationIndex 1, 2, 3
f.plot(diagnostic.test.mat = germ.mean.mat[,1:3],years=2006:2008)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat[,1:3],years=2006:2008)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)

dev.off()

# ---
# *Posterior predictive distribution ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "seedlingJan_sim")
n.iter=dim(y.sim)[1]

f = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density(x, from= x.min, to = x.max)
  return(dres)
}

pdf(file=paste0(outputDirectory,"ppc-germinationSeedBagExperiment.pdf"),height=6,width=6)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

for(i in 1:20){
  
  i.tmp = siteIndex[i]
  
  # for each site get index for year 1 germination
  index=data$siteGermination==i.tmp&data$germinationIndex%in%c(1)
  
  tmp=sweep(y.sim[,index], 2, data$totalJan[index], FUN = '/')
  sample.index = sample(1:n.iter,25)
  tmp=tmp[sample.index,]
  list.dens=apply(tmp,1,density,na.rm=TRUE)
  all.max.y=max(unlist(lapply(list.dens, "[", "y")))
  all.max.x=max(unlist(lapply(list.dens, "[", "x")))
  
  plot(NA,NA,
       ylim=c(0,all.max.y),xlim=c(0,1),
       xaxt='n',xlab='',ylab='',yaxt='n')
  
  for(j in 1:25){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,
          col=ifelse(p.chi.mat[i,1]>.95|p.chi.mat[i,1]<.05,
                     'orange','gray75'))
  }
  
  dres = f(data$seedlingJan[index]/data$totalJan[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  legend("topright",paste0(siteNames[i],"\n n=",length(data$totalJan[index])),bty='n')
  
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination in 2006 (year t+1 after seed production)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

dev.off()

# ---
# Intact seed counts ----
# ---

# ---
# *Posterior predictive checks ----
# ---

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.yobs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.ysim"))

n.iter = dim(chi2.obs)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=12)
for(j in 1:20){
  
  # classes = list(c(1,7,11),c(2,8,12),c(3,9),c(4,10),c(5),c(6))
  
  for(i in 1:12){
    #cIndex=classes[[i]]
    index = data$siteSurvival==j&data$compIndex%in%i
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


f=function(y.sim=chains,y.obs=data,n.obs=data2,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=12)
  for(j in 1:20){
    
    # classes = list(c(1,7,11),c(2,8,12),c(3,9),c(4,10),c(5),c(6))
    
    for(i in 1:12){
      #  cIndex=classes[[i]]
      index = data$siteSurvival==j&data$compIndex%in%i
      
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

sims=MCMCchains(mcmcSamples,params=c("y_sim"))
df=data$y
df2=data$seedStart

seeds.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)

# convert lists to matrix
f.convert = function(test.list){
  out.list <- list()
  for(i in 1:20){
    tmpYear = data$compIndex[data$siteSurvival==i]
    
    out.list[[i]] <- ifelse((1:12 %in% tmpYear), test.list[[i]],
                            ifelse(!(1:12 %in% tmpYear), NA, 0))
  }
  out.mat <- do.call(rbind,out.list)
  return(out.mat)
}

# return matrix objects
seeds.mean.mat<-seeds.mean
p.chi.mat<-f.convert(p.chi.list)

pdf(file=paste0(outputDirectory,"pvals-seedSurvivalSeedBagExperiment.pdf"),height=3,width=6)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

# seed survival s1 in 2006, 2007, and 2008 
# correspond to compIndex  c(1,7,11)

f.plot(diagnostic.test.mat = seeds.mean.mat[,c(1,7,11)],years=2006:2008)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat[,c(1,7,11)],years=2006:2008)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

# seed survival s2 in 2006, 2007, and 2008 
# correspond to compIndex  c(2,8,12)
# classes = list(c(1,7,11),c(2,8,12),c(3,9),c(4,10),c(5),c(6))

f.plot(diagnostic.test.mat = seeds.mean.mat[,c(2,8,12)],years=2006:2008)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat[,c(2,8,12)],years=2006:2008)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)


par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

# seed survival s3 in 2006, 2007
# correspond to compIndex  c(3,9)

f.plot(diagnostic.test.mat = seeds.mean.mat[,c(3,9)],years=2006:2007)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat[,c(3,9)],years=2006:2007)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)

dev.off()

# ---
# *Posterior predictive distribution ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "y_sim")
n.iter=dim(y.sim)[1]

f = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density(x, from= x.min, to = x.max)
  return(dres)
}

pdf(file=paste0(outputDirectory,"ppc-seedSurvivalSeedBagExperiment.pdf"),height=6,width=6)
par(mfrow=c(3,3),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

time.sample = 1:3
months.names = c("Four months", "Twelve months", "Sixteen months")

group.tmp=sample(1:20,3)

for(j in 1:3){
  
  j.tmp=group.tmp[j]
  
  n.samples = 25
  iter.ind = sample(1:n.iter,n.samples)
  
  for(i in 1:3){
    
    index = data$siteSurvival==j.tmp&data$compIndex==i
    y_sim=y.sim[,index]
    
    tmp=sweep(y_sim, 2, data$seedStart[index], FUN = '/')
    tmp=tmp[iter.ind,]
    
    list.dens=apply(tmp,1,density)
    all.max=max(unlist(lapply(list.dens, "[", "y")))
    
    plot(NA,NA,
         xlim=c(0,all.max),ylim=c(0,1),
         main=ifelse(j==1,months.names[i],''),
         xaxt='n',xlab='',ylab='',yaxt='n')
    for(h in 1:n.samples){
      m.max=max(tmp[h,],na.rm=TRUE)
      m.min=min(tmp[h,],na.rm=TRUE)
      
      dens.x = density(tmp[h,],from=m.min,to=m.max)
      
      lines(y=dens.x$x,x=dens.x$y,lwd=0.25,
            col=ifelse(p.chi.mat[i,1]>.95|p.chi.mat[i,1]<.05,'orange','gray75'))
    }
    
    p = data$y[index]/data$seedStart[index]
    m.max=max(p,na.rm=TRUE)
    m.min=min(p,na.rm=TRUE)
    
    dens.x = density(p,from=m.min,to=m.max, na.rm=TRUE)
    
    lines(y=dens.x$x,x=dens.x$y)
    
    ifelse(i%in%c(1),axis(2L),NA)
    ifelse(j%in%c(3),axis(1L),NA)
    if(i==1){
      legend("topright",
                     paste0(siteNames[group.tmp[j]],"\n n=",
                            length(data$seedStart[index])),bty='n') 
    } else {
      legend("topright",
             paste0("n=",
                    length(data$seedStart[index])),bty='n') 
      }
    
  }
  
  
}
mtext("Probability of a seed being intact (first round seed cohort)", side = 2, outer = TRUE, line = 2.2)
mtext("Density", side = 1, outer = TRUE, line = 2.2)

dev.off()

# ---
# *s0 ----
# ---

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.plot.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.plot.sim"))

n.iter = dim(chi2.obs)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=1)
for(j in 1:20){
  
  index = data$sitePlot==j
  tmp.chi2.obs=chi2.obs[,index]
  tmp.chi2.sim=chi2.sim[,index]
  fit.obs=apply(as.matrix(tmp.chi2.obs),1,sum)
  fit.sim=apply(as.matrix(tmp.chi2.sim),1,sum)
  p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
  p.pop[,1] = p.chi2.calc
  p.chi.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
}
p.chi.mat=do.call(rbind,p.chi.list)

data$n_sitePlot <- length(data$sitePlot)

f=function(y.sim=chains,y.obs=data,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  #p.pop = matrix(NA,nrow=n.iter,ncol=length(datasitePlot))
  for(j in 1:20){
    
      index = data$sitePlot==j
      
      tmp.obs=as.matrix(y.obs[index])
      tmp.sim=as.matrix(y.sim[,index])
      test.obs=model.fun(tmp.obs,na.rm=TRUE)
      test.sim=apply(tmp.sim,1,model.fun,na.rm=TRUE)
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)

    p.test.list[[j]] = mean(p.test.calc,na.rm=TRUE)
  }
  p.test.mat=do.call(rbind,p.test.list)
  return(p.test.mat)
}

sims=MCMCchains(mcmcSamples,params=c("y_sim"))
df=data$plotSeedlings

seedlings.mean=f(y.sim=sims,y.obs=df,model.fun=mean)

# return matrix objects
seedlings.mean.mat<-seedlings.mean
p.chi.mat<-do.call(rbind,p.chi.list)


pdf(file=paste0(outputDirectory,"pvals-seedlingsEmergingPlots.pdf"),height=6,width=6)


par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))


f.plot(diagnostic.test.mat = seedlings.mean.mat,years=2008)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat,years=2008)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)

dev.off()

# ---
# *s0 ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "plotSeedlings_sim")

pdf(file=paste0(outputDirectory,"ppc-seedlingsEmergingPlots.pdf"),height=6,width=6)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  n.samples = 25
  iter.ind = sample(1:n.iter,n.samples)
  
  i.tmp = siteIndex[i]
  
  index=data$sitePlot==i.tmp
  y_sim=y.sim[,index]
  
  tmp=y_sim
  tmp=tmp[iter.ind,]
  
  list.dens=apply(tmp,1,density)
  
  # also density of observed
  p.obs = data$plotSeedlings[index]
  m.max=max(p.obs)
  m.min=min(p.obs)
  dens.x.obs = density(p.obs,from=m.min,to=m.max)
  list.dens[[26]] <- dens.x.obs
  
  all.max.y=max(unlist(lapply(list.dens, "[", "y")))
  all.max.x=max(unlist(lapply(list.dens, "[", "x")))
  
  plot(NA,NA,
       ylim=c(0,all.max.y),xlim=c(0,all.max.x),
       main='',xaxt='n',xlab='',ylab='',yaxt='n')
  
  for(h in 1:n.samples){
    m.max=max(tmp[h,])
    m.min=min(tmp[h,])
    
    dens.x = density(tmp[h,],from=m.min,to=m.max)
    
    lines(x=dens.x$x,y=dens.x$y,lwd=0.25,col='orange')
  }
  
  p = data$plotSeedlings[index]
  m.max=max(p)
  m.min=min(p)
  
  dens.x = density(p,from=m.min,to=m.max)
  
  lines(x=dens.x$x,y=dens.x$y)
  
  legend("topright",paste0(siteNames[i],"\n n=",sum(data$plotSeedlings[index])),bty='n')
  
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Seeds emerging in permanent plots in January 2008, \n from seeds produced in year t-1 and t-2", side = 1, outer = TRUE, line = 3.5)
mtext("Density", side = 2, outer = TRUE, line = 2.2)
dev.off()

