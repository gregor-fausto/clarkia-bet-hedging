####
####
# Script for model checks for seedling survival to fruiting
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

# - +load libraries ----
library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)

# - Read in what's needed ----

mcmcDirectory = "outputs/002_fitStatisticalModels/mcmcSamples/"
fullDataDirectory = "outputs/001_prepareDataForModels/"
outputDirectory = "outputs/004_checkStatisticalModels/"

# - +Read in data ----
fullDataDirectory <- paste0(fullDataDirectory,list.files(fullDataDirectory))
data <- readRDS(fullDataDirectory[[grep("seedlingFruitingPlantCountsPermanentPlots.RDS",fullDataDirectory)]])

# - +Read in MCMC samples ----
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedlingSurvivalPosteriorSamples.RDS",mcmcSampleDirectory)]])

# - +Read in site names ----
siteAbiotic <- read.csv("data/siteAbioticData.csv",header=TRUE) 
siteNames = unique(siteAbiotic$site)

# create variable that arranges populations by easting
siteIndex <- order(siteAbiotic$easting,decreasing=FALSE)


# - Function to plot Bayesian p-values ----
f.plot = function(diagnostic.test.mat,years){
  
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(min(years),max(years)),
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
  
  axis(1, seq(min(years),max(years),by=2),
       labels = seq(min(years),max(years),by=2), 
       las = 1,
       col = NA, col.ticks = 1, cex.axis = 1)
  axis(2, seq(0,1,by=.2),
       seq(0,1,by=.2), las = 1,
       col = NA, col.ticks = 1, cex.axis = 1)
}

# - Posterior Predictive Checks ---
# - ++Get chi-squared calculated in JAGS ----
chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sim"))
n.iter=dim(chi2.obs)[1]

# - ++Calculate Bayesian p-value per population/year ----
# chi-squared
p.pop = matrix(NA,nrow=n.iter,ncol=data$n_siteYearIndex)
for(j in 1:20){
  
  tmpYear = data$year_observed[data$site_observed==j]
  
  for(i in 1:length(tmpYear)){
    index = data$site==j&data$year==tmpYear[i]
    tmp.chi2.obs=chi2.obs[,index]
    tmp.chi2.sim=chi2.sim[,index]
    
    if(is.matrix(tmp.chi2.obs)){
      fit.obs=apply(tmp.chi2.obs,1,sum)
      fit.sim=apply(tmp.chi2.sim,1,sum)
    } else if(is.vector(tmp.chi2.obs)){
      fit.obs=sum(tmp.chi2.obs)
      fit.sim=sum(tmp.chi2.sim)
    }
    
    p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
    index = data$siteYearIndex_observed[data$site_observed==j&data$year_observed==tmpYear[i]] 
    p.pop[,index] = p.chi2.calc
  }
}

p.chi.mat = apply(p.pop,2,mean,na.rm=TRUE)

# put in list
p.chi.list = list()
for(j in 1:20){
  tmpIndex = data$siteYearIndex_observed[data$site_observed==j]
  p.chi.list[[j]] = p.chi.mat[tmpIndex]
}

# - +++function to calculate Bayesian p-value for specific test statistic ----
f=function(y.sim=chains,y.obs=data,n.obs=data2,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=data$n_siteYearIndex)
  for(j in 1:20){
    
    tmpYear = data$year_observed[data$site_observed==j]
    
    for(i in 1:length(tmpYear)){
      
      index = data$site==j&data$year==tmpYear[i]
      
      tmp.obs=as.matrix(y.obs[index])
      tmp.sim=as.matrix(y.sim[,index])
      test.obs=model.fun(tmp.obs,na.rm=TRUE)
      test.sim=apply(tmp.sim,1,model.fun,na.rm=TRUE)
      
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      index = data$siteYearIndex_observed[data$site_observed==j&data$year_observed==tmpYear[i]] 
      p.pop[,index] = p.test.calc
      
    }
  }
  
  p.test.mat = apply(p.pop,2,mean,na.rm=TRUE)
  
  p.test.list = list()
  for(j in 1:20){
    tmpIndex = data$siteYearIndex_observed[data$site_observed==j]
    p.test.list[[j]] = p.test.mat[tmpIndex]
  }
  
  return(p.test.list)
}

# - ++get simulated observations and observed data ----
sims=MCMCchains(mcmcSamples,params=c("fruitplNumber_sim"))
df=data$fruitplNumber
df2=data$seedlingNumber

# - ++Calculate Bayesian p-value ----
# mean
seedlings.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)

# - +++function to convert list to matrix ----
f.convert = function(test.list){
  out.list <- list()
  for(i in 1:20){
    tmpYear = data$year_observed[data$site_observed==i]
    
    out.list[[i]] <- ifelse((1:15 %in% tmpYear), test.list[[i]], 
                            ifelse(!(1:15 %in% tmpYear), NA, 0))
  }
  out.mat <- do.call(rbind,out.list)
  return(out.mat)
}

# return matrix objects
p.chi.mat<-f.convert(p.chi.list)
seedlings.mean.mat<-f.convert(seedlings.mean)

# - +Plot Bayesian p-vals ----

pdf(file=paste0(outputDirectory,"pvals-seedlingSurvivalFruiting-chi.pdf"),height=4,width=6)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

f.plot(diagnostic.test.mat = seedlings.mean.mat,years=2006:2020)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat,years=2006:2020)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)

dev.off()

# *Posterior predictive distribution ----
# - ++get simulated values and observations----
y.sim=MCMCchains(mcmcSamples, params = "fruitplNumber_sim")
y.obs=data$fruitplNumber
n.obs=data$seedlingNumber
years=2006:2020
n.samples = 25
n.sim=5000

# - +Plot posterior predictive distributionin all years----
pdf(file=paste0(outputDirectory,"ppc-seedlingSurvivalFruiting.pdf"),height=6,width=6)

# for each year
for(i in 1:length(years)){
  
  par(mfrow = c(4,5),
      oma = c(4,5,0,0) + 0.1,
      mar = c(1,0,1,1) + 0.1)
  
  n.samples = 25
  iter.ind = sample(1:n.iter,n.samples)
  
  tmpSite = data$site_observed[data$year_observed==i]
  
  # for each population
  for(j in 1:20){
    
    # order populations geographically
    j.tmp = siteIndex[j]
    if(j.tmp %in% tmpSite){
      
      # get index
      index=data$site==j.tmp&data$year==i
      
      # get simulated and observed values
      p.obs=y.obs[index]
      p.sim=y.sim[,index]
      if(is.matrix(p.sim)){
        # calculate density for simulated values 
        p.sim=p.sim[iter.ind,]
        list.dens=apply(p.sim,1,density,na.rm=TRUE)
        # get max values for plot dimensions
        all.max.y=max(unlist(lapply(list.dens, "[", "y")))
        all.max.x=max(unlist(lapply(list.dens, "[", "x")))
        
        plot(NA,NA,
             ylim=c(0,all.max.y),xlim=c(0,all.max.x),
             xaxt='n',xlab='',ylab='',yaxt='n')
        # for each dataset, plot the posterior predictive distribution
        # coloring it according to whether or not the p-value is extreme
        # use the p-values for the mean test statistic
        for(h in 1:n.samples){
          m.max=max(p.sim[h,],na.rm=TRUE)
          m.min=min(p.sim[h,],na.rm=TRUE)
          dens.x = density(p.sim[h,],from=m.min,to=m.max,na.rm=TRUE)
          lines(x=dens.x$x,y=dens.x$y,lwd=0.25,
                col=if(p.chi.mat[j.tmp,i]>.95|p.chi.mat[j.tmp,i]<.05&seedlings.mean.mat[j.tmp,i]>.95|seedlings.mean.mat[j.tmp,i]<.05){'orange'}
                else if(p.chi.mat[j.tmp,i]>.95|p.chi.mat[j.tmp,i]<.05&seedlings.mean.mat[j.tmp,i]<.95&seedlings.mean.mat[j.tmp,i]>.05){'purple'}
                else{ 'gray75'})
        }
        
        # get the distribution for the observed data
        # plot it as a black line overlaying the ppd
        p = data$fruitplNumber[index]
        m.max=max(p,na.rm=TRUE)
        m.min=min(p,na.rm=TRUE)
        dens.x = density(p,from=m.min,to=m.max,na.rm=TRUE)
        lines(x=dens.x$x,y=dens.x$y)
        
        # plot single observations as a dotted line
        if(length(unique(p))==1){ abline(v=unique(p),lty='dotted',col='black')}
        
        # lines below are for the case where there is only 1 observation
        # plot a single dotted line and make a discrete density plot
      } else if(is.vector(p.sim)){
        
        p.sim = p.sim[iter.ind]
        all.max.x= max(density(p.sim,na.rm=TRUE)$x)
        all.max.y= max(density(p.sim,na.rm=TRUE)$y)
        
        plot(NA,NA,
             ylim=c(0,all.max.y),xlim=c(0,all.max.x),
             xaxt='n',xlab='',ylab='',yaxt='n')
        
        # plot observed value
        p = data$fruitplNumber[index]
        abline(v=p,lty='dotted',col='black')
        
        # overlay simulated values
        tb = table(p.sim)
        tb.x = as.numeric(names(tb))
        segments(x0=tb.x,y0=0,y1=tb/sum(tb),
                 col=if(p.chi.mat[j.tmp,i]>.95|p.chi.mat[j.tmp,i]<.05&seedlings.mean.mat[j.tmp,i]>.95|seedlings.mean.mat[j.tmp,i]<.05){'orange'}
                 else if(p.chi.mat[j.tmp,i]>.95|p.chi.mat[j.tmp,i]<.05&seedlings.mean.mat[j.tmp,i]<.95&seedlings.mean.mat[j.tmp,i]>.05){'purple'}
                 else{ 'gray75'})
      }
      
      axis(2,cex=.25,tick=FALSE,line=-1)
      axis(1,cex=.25,tick=FALSE,line=-1)
      
      legend("topright",paste0(siteNames[j.tmp],"\n n=",length(p)),bty='n')
    } else {
      plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
    }
  }
  
  mtext(paste0("Fruiting plant counts (",years[i],")"), side = 1, outer = TRUE, line = 1.5)
  mtext("Density", side = 2, outer = TRUE, line = 2.2)
}

dev.off()

# - ++ Plot for 2007 ----
# plot 4 sites
# pick 3 ahead of time (site 1, 3, 16, and then randomly sample 1 more)
# here graph the plot-level simulated values against the observed values

sims=MCMCchains(mcmcSamples,params=c("fruitplNumber_sim"))
df=data$fruitplNumber
df2=data$seedlingNumber

pdf(file=paste0(outputDirectory,"ppc-seedlingSurvivalFruiting-perPlot.pdf"),height=6,width=6)

par(mfrow = c(2,2),
    oma = c(2,2,0,0) + 0.1,
    mar = c(1,0,1,1) + 0.1)

for(j in c(1,3,16,sample(c(2,4:15,17:20),1))){
  
  j.tmp = siteIndex[j]
  if(j.tmp %in% tmpSite){
    
    index=data$site==j.tmp&data$year==2
    
    p.obs=df[index]
    p.sim=sims[,index]
    
    plot(NA,NA,
         xlim=c(1,length(p.obs)+1),
         ylim=c(0,max(p.sim)),
         xaxt='n',xlab='',ylab='',yaxt='n')
    
    rect(xleft=seq(1,50,by=2),
         xright=seq(1,50,by=2)+1,
         ybottom=-1000,ytop=1000,col='gray99',border='gray99')
    
    for(i in 1:length(p.obs)){
      index.rand = sample(1:45000,50)
      tmp<-table(p.sim[index.rand,i])
      segments(y0=as.numeric(names(tmp)),y1=as.numeric(names(tmp)),
               x0=i, x1=i+(tmp/max(tmp))*.8)
    }
    box()
    points(1:length(p.obs),p.obs,pch=16)
    
    axis(2,cex=.25,tick=FALSE,line=-1)
    axis(1,cex=.25,tick=FALSE,line=-1)
    
    legend("topright",paste0(siteNames[j.tmp],"\n n=",length(p.obs)),bty='n')
    
  }
  
}

mtext(paste0("Fruiting plant counts"), side = 2, outer = TRUE, line = 1)
mtext("Observation (index)", side = 1, outer = TRUE, line = 0.2)

dev.off()

# - ++ Plot for LCW ----
# here graph the plot-level simulated values against the observed values
# for all years for LCW

pdf(file=paste0(outputDirectory,"ppc-seedlingSurvivalFruiting-populationLevel.pdf"),height=6,width=6)

par(mfrow = c(3,5),
    oma = c(2,2,0,0) + 0.1,
    mar = c(1,0,1,1) + 0.1)

# get LCW
for(j in 3){
  
  j.tmp = siteIndex[j]
  if(j.tmp %in% tmpSite){
    
    for(i in 1:15){
      index=data$site==j.tmp&data$year==i
      
      p.obs=df[index]
      p.sim=as.matrix(sims[,index])
      
      if(sum(index)>0){
        plot(NA,NA,
             xlim=c(1,length(p.obs)+1),
             ylim=c(0,max(p.sim)),
             xaxt='n',xlab='',ylab='',yaxt='n')
        
        rect(xleft=seq(1,50,by=2),
             xright=seq(1,50,by=2)+1,
             ybottom=-1000,ytop=1000,col='gray99',border='gray99')
        
        for(h in 1:length(p.obs)){ 
          index.rand = sample(1:45000,50)
          tmp<-table(p.sim[index.rand,h])
          segments(y0=as.numeric(names(tmp)),y1=as.numeric(names(tmp)),
                   x0=h, x1=h+(tmp/max(tmp))*.8)
        }
        box()
        points(1:length(p.obs),p.obs,pch=16)
        
        axis(2,cex=.25,tick=FALSE,line=-1)
        axis(1,cex=.25,tick=FALSE,line=-1)
        
        legend("topright",paste0(years[i],"\n",siteNames[j.tmp],"\n n=",length(p.obs)),bty='n')
        
      }
    } 
  }
  else {
    plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
  }
}

mtext(paste0("Fruiting plant counts"), side = 2, outer = TRUE, line = 1)
mtext("Observation (index)", side = 1, outer = TRUE, line = 0.2)

dev.off()