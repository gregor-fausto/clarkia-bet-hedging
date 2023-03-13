####
####
# Script for model checks for seedling survival to fruiting
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

mcmcDirectory = "outputs/002_fitStatisticalModels/mcmcSamples/"
fullDataDirectory = "outputs/001_prepareDataForModels/"
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


# ---
# ---
# - Posterior Predictive Checks ---
# ---
# ---

# ---
# Test statistics ----
# ---
chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.obs"))

n.iter=dim(chi2.obs)[1]
subsample = sample(n.iter,5000)

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.obs"))[subsample,]
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sim"))[subsample,]

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

p.chi.list = list()
for(j in 1:20){
  tmpIndex = data$siteYearIndex_observed[data$site_observed==j]
  p.chi.list[[j]] = p.chi.mat[tmpIndex]
}


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

sims=MCMCchains(mcmcSamples,params=c("fruitplNumber_sim"))[subsample,]
df=data$fruitplNumber
df2=data$seedlingNumber

seedlings.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)

# convert lists to matrix
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

var.mat = matrix(NA,nrow=20,ncol=15)
for(i in 1:20){
  for(j in 1:15){
    index=data$site==i&data$year==j
    n_sdlg=data$seedlingNumber[index]
    n_frts=data$fruitplNumber[index]
    #p=n_frts/n_sdlg
    #  var.mat[i,j]=var(n_sdlg,na.rm=TRUE)
    var.mat[i,j]=sum(n_sdlg,na.rm=TRUE)
  }
}
dev.off()


# ---
# Graphical checks ----
# ---


y.sim=MCMCchains(mcmcSamples, params = "fruitplNumber_sim")
y.obs=data$fruitplNumber
n.obs=data$seedlingNumber
years=2006:2020
n.samples = 50
n.sim=5000


pdf(file=paste0(outputDirectory,"ppc-seedlingSurvivalFruiting.pdf"),height=6,width=6)

for(i in 1:length(years)){
  
  par(mfrow = c(4,5),
      oma = c(4,5,0,0) + 0.1,
      mar = c(1,0,1,1) + 0.1)
  
  n.samples = 50
  iter.ind = sample(1:length(subsample),n.samples)
  
  tmpSite = data$site_observed[data$year_observed==i]
  
  for(j in 1:20){
    
    j.tmp = siteIndex[j]
    if(j.tmp %in% tmpSite){
      
      index=data$site==j.tmp&data$year==i
      
      p.obs=y.obs[index]
      p.sim=y.sim[,index]
      if(is.matrix(p.sim)){
        p.sim=p.sim[iter.ind,]
        list.dens=apply(p.sim,1,density,na.rm=TRUE)
        all.max.y=max(unlist(lapply(list.dens, "[", "y")))
        all.max.x=max(unlist(lapply(list.dens, "[", "x")))
        
        plot(NA,NA,
             ylim=c(0,all.max.y),xlim=c(0,all.max.x),
             xaxt='n',xlab='',ylab='',yaxt='n')
        for(h in 1:n.samples){
          m.max=max(p.sim[h,],na.rm=TRUE)
          m.min=min(p.sim[h,],na.rm=TRUE)
          
          dens.x = density(p.sim[h,],from=m.min,to=m.max,na.rm=TRUE)
          lines(x=dens.x$x,y=dens.x$y,lwd=0.25,
                col=ifelse(p.chi.mat[j.tmp,i]>.95|p.chi.mat[j.tmp,i]<.05,
                           'orange','gray75'))
        }
        
        p = data$fruitplNumber[index]
        m.max=max(p,na.rm=TRUE)
        m.min=min(p,na.rm=TRUE)
        
        dens.x = density(p,from=m.min,to=m.max,na.rm=TRUE)
        
        lines(x=dens.x$x,y=dens.x$y)
        
        if(length(unique(p))==1){ abline(v=unique(p),lty='dotted',col='black')}
        
      } else if(is.vector(p.sim)){
        # lines below are for the case where there is only 1 observation
        
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
                 col=ifelse(p.chi.mat[j.tmp,i]>.95|p.chi.mat[j.tmp,i]<.05,
                            'orange','gray75'))
        
      }
      
      # ifelse(i%in%c(1),axis(2L),NA)
      
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


