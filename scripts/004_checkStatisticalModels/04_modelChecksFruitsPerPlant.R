####
####
# Script for model checks for seeds per fruit
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

fullDataDirectory = "outputs/001_prepareDataForModels/"
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
fullDataDirectory <- paste0(fullDataDirectory,list.files(fullDataDirectory))
data <- readRDS(fullDataDirectory[[grep("fruitsPerPlantAllPlants.RDS",fullDataDirectory)]])

# - +Read in MCMC samples ----
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("fruitsPerPlantPosteriorSamples.RDS",mcmcSampleDirectory)]])

# - +Read in site names ----
siteAbiotic <- read.csv("data/siteAbioticData.csv",header=TRUE)
# create variable that arranges populations by easting
siteIndex <- order(siteAbiotic$easting,decreasing=FALSE)
siteNames = unique(siteAbiotic$site)[siteIndex]

# - Function to plot p-vals ----

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

# ---
# - Total fruit equivalents ----
# - +Simulate observations ----
# ---

z_tfe=MCMCchains(mcmcSamples, params = "z_tfe")

n.sim = 5000
draws.iter=sample(1:dim(z_tfe)[1],n.sim)
z_tfe.sub = z_tfe[draws.iter,]
n.obs = dim(z_tfe.sub)[2]
# for each observation
# simulate observations from a Poisson
y_tfe.sim = matrix(NA,nrow=n.sim,ncol=dim(z_tfe.sub)[2])
for(i in 1:n.sim){
  y_tfe.sim[i,]<-rpois(n.obs,z_tfe.sub[i,])
}
y.sim = y_tfe.sim

# --
# - Posterior Predictive Checks  ----
# --

# ---
# +Test statistics: Total fruit equivalents per plant ----
# ---

y.obs=data$y_tfe
years=2006:2012

chi2.obs=(sweep(z_tfe.sub, 2, y.obs, FUN = '-')^2)/z_tfe.sub
chi2.sim=((y_tfe.sim-z_tfe.sub)^2)/z_tfe.sub

n.iter = dim(chi2.obs)[1]

p.pop = matrix(NA,nrow=n.iter,ncol=data$n_siteYearIndex_tfe)
for(j in 1:20){
  
  tmpYear = data$year_tfe_observed[data$site_tfe_observed==j]
  
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
    index = data$siteYearIndex_tfe_observed[data$site_tfe_observed==j&data$year_tfe_observed==tmpYear[i]] 
    p.pop[,index] = p.chi2.calc
  }
}

p.chi.mat = apply(p.pop,2,mean,na.rm=TRUE)

p.chi.list = list()
for(j in 1:20){
  tmpIndex = data$siteYearIndex_tfe_observed[data$site_tfe_observed==j]
  p.chi.list[[j]] = p.chi.mat[tmpIndex]
}

f=function(y.sim=chains,y.obs=data,n.obs=data2,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=data$n_siteYearIndex_tfe)
  for(j in 1:20){
    
    tmpYear = data$year_tfe_observed[data$site_tfe_observed==j]
    
    for(i in 1:length(tmpYear)){
      
      index = data$site==j&data$year==tmpYear[i]
      tmp.obs=as.matrix(y.obs[index])
      tmp.sim=as.matrix(y.sim[,index])
      test.obs=model.fun(tmp.obs,na.rm=TRUE)
      test.sim=apply(tmp.sim,1,model.fun,na.rm=TRUE)
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      index = data$siteYearIndex_tfe_observed[data$site_tfe_observed==j&data$year_tfe_observed==tmpYear[i]] 
      p.pop[,index] = p.test.calc
      
    }
  }
  
  p.test.mat = apply(p.pop,2,mean,na.rm=TRUE)
  
  p.test.list = list()
  for(j in 1:20){
    tmpIndex = data$siteYearIndex_tfe_observed[data$site_tfe_observed==j]
    p.test.list[[j]] = p.test.mat[tmpIndex]
  }
  
  return(p.test.list)
}

sims=y.sim
df=data$y_tfe

tfe.mean=f(y.sim=sims,y.obs=df,model.fun=mean)

# convert lists to matrix
f.convert = function(test.list){
  out.list <- list()
  for(i in 1:20){
    tmpYear = data$year_tfe_observed[data$site_tfe_observed==i]
    
    out.list[[i]] <- ifelse((1:7 %in% tmpYear), test.list[[i]],
                            ifelse(!(1:7 %in% tmpYear), NA, 0))
  }
  out.mat <- do.call(rbind,out.list)
  return(out.mat)
}

# return matrix objects
p.chi.mat<-f.convert(p.chi.list)
tfe.mean.mat<-f.convert(tfe.mean)


pdf(file=paste0(outputDirectory,"pvals-totalFruitEquivalentsPerPlant.pdf"),height=3,width=6)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

f.plot(diagnostic.test.mat = tfe.mean.mat,years=2006:2012)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat,years=2006:2012)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)

dev.off()

# ---
# Graphical checks: Total fruit equivalents per plant ----
# use 2009 as the example
# ---

y.obs=data$y_tfe
n.samples = 50

years=2006:2012
plot.yr = 2009

pdf(file=paste0(outputDirectory,"ppc-totalFruitEquivalentsPerPlant.pdf"),height=4,width=6)

par(mfrow = c(2,4),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(2,1,0))

for(i in (1:length(years))[ years %in% plot.yr]){
  if(colSums(p.chi.mat>.95,na.rm=TRUE)[i]==0){
    
  } else {
    
    iter.ind = sample(1:n.sim,n.samples)
    
    # get the sites with poor fit
    site.index<-(1:20)[(p.chi.mat>.95&!is.na(p.chi.mat))[,i]]
    # add two random sites
    add.index <- c(sample((1:20)[!(1:20%in%site.index)],4))
    for(j in c(site.index,add.index)){
      
      index=data$site==j&data$year==i
      
      p.obs=y.obs[index]
      p.sim=as.matrix(y.sim[,index])
      p.sim=as.matrix(p.sim[iter.ind,])
      
      if (dim(p.sim)[2]<2) {NA} else {
        
        list.dens=apply(p.sim,1,density,na.rm=TRUE)
        all.max.y=max(unlist(lapply(list.dens, "[", "y")))
        all.max.x=max(unlist(lapply(list.dens, "[", "x")))
        
        m.max=max(p.obs,na.rm=TRUE)
        m.min=min(p.obs,na.rm=TRUE)
        dens.obs = density(p.obs,from=m.min,to=m.max,na.rm=TRUE)
        
        all.max.y = max(all.max.y,max(dens.obs$y))
        all.max.x = max(all.max.x,max(dens.obs$x))
        
        plot(NA,NA,
             ylim=c(0,all.max.y),xlim=c(0,all.max.x),
             xaxt='n',xlab='',ylab='',yaxt='n')
        for(h in 1:n.samples){
          m.max=max(p.sim[h,],na.rm=TRUE)
          m.min=min(p.sim[h,],na.rm=TRUE)
          
          dens.x = density(p.sim[h,],from=m.min,to=m.max,na.rm=TRUE)
          lines(x=dens.x$x,y=dens.x$y,lwd=0.25,
                col=ifelse(j %in% site.index,'orange','gray75'))
        }
        
        
        lines(x=dens.obs$x,y=dens.obs$y)
        
        ifelse(i%in%c(1),axis(2L),NA)
        
        axis(2,cex=.25,tick=FALSE,line=-1)
        axis(1,cex=.25,tick=FALSE,line=-1)
        
        legend("topright",paste0(siteAbiotic$site[j],': ',years[i],"\n n=",length(p.obs)),bty='n')
        
      }
    }
    
    mtext(paste0("Counts of total fruit equivalents per plant"), side = 1, 
          outer = TRUE, line = .75)
    mtext("Density", side = 2, outer = TRUE, line = 1.2)
  }
}

dev.off()


# ---
# Simulate observations: Total fruits per plant ----
# ---

z_tot=MCMCchains(mcmcSamples, params = "z_tot")

n.sim = 5000
draws.iter=sample(1:dim(z_tot)[1],n.sim)
z_tot.sub = z_tot[draws.iter,]
n.obs = dim(z_tot.sub)[2]
# for each observation
# simulate observations from a Poisson
y_tot.sim = matrix(NA,nrow=n.sim,ncol=dim(z_tot.sub)[2])
for(i in 1:n.sim){
  y_tot.sim[i,]<-rpois(n.obs,z_tot.sub[i,])
}
y_tot.sim = y_tot.sim

# ---
# Test statistics: Total fruits ----
# ---

y.obs=data$y_tot
years=2013:2020

chi2.obs=(sweep(z_tot.sub, 2, y.obs, FUN = '-')^2)/z_tot.sub
chi2.sim=((y_tot.sim-z_tot.sub)^2)/z_tot.sub

n.iter = dim(chi2.obs)[1]

p.pop = matrix(NA,nrow=n.iter,ncol=data$n_siteYearIndex_tot)
for(j in 1:20){
  
  tmpYear = data$year_tot_observed[data$site_tot_observed==j]
  
  for(i in 1:length(tmpYear)){
    index = data$site2==j&data$year2==tmpYear[i]
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
    index = data$siteYearIndex_tot_observed[data$site_tot_observed==j&data$year_tot_observed==tmpYear[i]] 
    p.pop[,index] = p.chi2.calc
  }
}

p.chi.mat = apply(p.pop,2,mean,na.rm=TRUE)

p.chi.list = list()
for(j in 1:20){
  tmpIndex = data$siteYearIndex_tot_observed[data$site_tot_observed==j]
  p.chi.list[[j]] = p.chi.mat[tmpIndex]
}



f=function(y.sim=chains,y.obs=data,n.obs=data2,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=data$n_siteYearIndex_tot)
  for(j in 1:20){
    
    tmpYear = data$year_tot_observed[data$site_tot_observed==j]
    
    for(i in 1:length(tmpYear)){
      
      index = data$site2==j&data$year2==tmpYear[i]
      tmp.obs=as.matrix(y.obs[index])
      tmp.sim=as.matrix(y.sim[,index])
      test.obs=model.fun(tmp.obs,na.rm=TRUE)
      test.sim=apply(tmp.sim,1,model.fun,na.rm=TRUE)
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      index = data$siteYearIndex_tot_observed[data$site_tot_observed==j&data$year_tot_observed==tmpYear[i]] 
      p.pop[,index] = p.test.calc
      
    }
  }
  
  p.test.mat = apply(p.pop,2,mean,na.rm=TRUE)
  
  p.test.list = list()
  for(j in 1:20){
    tmpIndex = data$siteYearIndex_tot_observed[data$site_tot_observed==j]
    p.test.list[[j]] = p.test.mat[tmpIndex]
  }
  
  return(p.test.list)
}

# convert lists to matrix
f.convert = function(test.list){
  out.list <- list()
  for(i in 1:20){
    tmpYear = data$year_tot_observed[data$site_tot_observed==i]
    
    out.list[[i]] <- ifelse((1:8 %in% tmpYear), test.list[[i]],
                            ifelse(!(1:8 %in% tmpYear), NA, 0))
  }
  out.mat <- do.call(rbind,out.list)
  return(out.mat)
}

sims=y_tot.sim
df=data$y_tot

tot.mean=f(y.sim=sims,y.obs=df,model.fun=mean)

p.chi.mat<-f.convert(p.chi.list)
tot.mean.mat<-f.convert(tot.mean)

pdf(file=paste0(outputDirectory,"pvals-totalFruitsPerPlant.pdf"),height=3,width=6)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

f.plot(diagnostic.test.mat = tot.mean.mat,years=2013:2020)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat,years=2013:2020)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)

dev.off()


# ---
# Graphical checks: Total fruits per plant ----
# ---


y.obs=data$y_tot
n.samples = 50

years=2013:2020
plot.yr = 2013

pdf(file=paste0(outputDirectory,"ppc-totalFruitsPerPlant.pdf"),height=2,width=4)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(2,1,0))

for(i in (1:length(years))[ years %in% plot.yr]){
  if(colSums(p.chi.mat>.95,na.rm=TRUE)[i]==0){
    
  } else {
    
    iter.ind = sample(1:n.iter,n.samples)
    
    # get the sites with poor fit
    site.index<-(1:20)[(p.chi.mat>.95&!is.na(p.chi.mat))[,i]]
    # add two random sites
    add.index <- c(sample((1:20)[!(1:20%in%site.index)],1))
    for(j in c(site.index,add.index)){
      
      index=data$site2==j&data$year2==i
      
      p.obs=y.obs[index]
      p.sim=as.matrix(y_tot.sim[,index])
      p.sim=as.matrix(p.sim[iter.ind,])
      
      if(is.matrix(p.sim)&dim(p.sim)[2]>1){
        
        list.dens=apply(p.sim,1,density,na.rm=TRUE)
        all.max.y=max(unlist(lapply(list.dens, "[", "y")))
        all.max.x=max(unlist(lapply(list.dens, "[", "x")))
        
        m.max=max(p.obs,na.rm=TRUE)
        m.min=min(p.obs,na.rm=TRUE)
        dens.obs = density(p.obs,from=m.min,to=m.max,na.rm=TRUE)
        
        all.max.y = max(all.max.y,max(dens.obs$y))
        all.max.x = max(all.max.x,max(dens.obs$x))
        
        plot(NA,NA,
             ylim=c(0,all.max.y),xlim=c(0,all.max.x),
             xaxt='n',xlab='',ylab='',yaxt='n')
        for(h in 1:n.samples){
          m.max=max(p.sim[h,],na.rm=TRUE)
          m.min=min(p.sim[h,],na.rm=TRUE)
          
          dens.x = density(p.sim[h,],from=m.min,to=m.max,na.rm=TRUE)
          lines(x=dens.x$x,y=dens.x$y,lwd=0.25,
                col=ifelse(j %in% site.index,'orange','gray75'))
        }
        
        p = data$y_tot[index]
        
        lines(x=dens.obs$x,y=dens.obs$y)
        
        if(length(unique(p))==1){ abline(v=unique(p),lty='dotted',col='black')}
        
      } else if(is.matrix(p.sim)&dim(p.sim)[2]==1){
        # lines below are for the case where there is only 1 observation
        
        p.sim = p.sim
        all.max.x= max(density(p.sim,na.rm=TRUE)$x)
        all.max.y= max(density(p.sim,na.rm=TRUE)$y)
        
        plot(NA,NA,
             ylim=c(0,all.max.y),xlim=c(0,all.max.x),
             xaxt='n',xlab='',ylab='',yaxt='n')
        
        # plot observed value
        p = data$y_tot[index]
        abline(v=p,lty='dotted',col='black')
        
        # overlay simulated values
        tb = table(p.sim)
        tb.x = as.numeric(names(tb))
        segments(x0=tb.x,y0=0,y1=tb/sum(tb),
                 col=ifelse(j %in% site.index,'orange','gray75'))
        
      }
      
      axis(2,cex=.25,tick=FALSE,line=-1)
      axis(1,cex=.25,tick=FALSE,line=-1)
      
      legend("topright",paste0(siteAbiotic$site[j],': ',years[i],"\n n=",length(p)),bty='n')
      
    }
    
    mtext(paste0("Counts of total fruits per plant"), side = 1,           
          outer = TRUE, line = .75)
    mtext("Density", side = 2, outer = TRUE, line = 1.2)
  }
}
dev.off()



# ---
# Simulate observations: Damaged fruits per plant ----
# ---

p_dam=MCMCchains(mcmcSamples, params = "prop_dam")

n.sim = 5000
draws.iter=sample(1:dim(p_dam)[1],n.sim)
p_dam.sub = p_dam[draws.iter,]
n.obs = dim(p_dam.sub)[2]
# for each observation
# simulate observations from a Poisson
y_dam.sim = matrix(NA,nrow=n.sim,ncol=dim(p_dam.sub)[2])
for(i in 1:n.obs){
  y_dam.sim[,i]<-rbinom(n.sim,data$y_tot[i],p_dam.sub[,i])
}
y_dam.sim = y_dam.sim


# ---
# Test statistics: Damaged fruits ----
# ---

y.obs=data$y_dam
years = 2013:2020

chi2.obs=(sweep(data$y_tot*p_dam.sub, 2, y.obs, FUN = '-')^2)/data$y_tot*p_dam.sub
chi2.sim=((y_dam.sim-data$y_tot*p_dam.sub)^2)/data$y_tot*p_dam.sub

n.iter = dim(chi2.obs)[1]

p.pop = matrix(NA,nrow=n.iter,ncol=data$n_siteYearIndex_tot)
for(j in 1:20){
  
  tmpYear = data$year_tot_observed[data$site_tot_observed==j]
  
  for(i in 1:length(tmpYear)){
    index = data$site2==j&data$year2==tmpYear[i]
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
    index = data$siteYearIndex_tot_observed[data$site_tot_observed==j&data$year_tot_observed==tmpYear[i]] 
    p.pop[,index] = p.chi2.calc
  }
}

p.chi.mat = apply(p.pop,2,mean,na.rm=TRUE)

p.chi.list = list()
for(j in 1:20){
  tmpIndex = data$siteYearIndex_tot_observed[data$site_tot_observed==j]
  p.chi.list[[j]] = p.chi.mat[tmpIndex]
}



f=function(y.sim=chains,y.obs=data,n.obs=data2,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=data$n_siteYearIndex_tot)
  for(j in 1:20){
    
    tmpYear = data$year_tot_observed[data$site_tot_observed==j]
    
    for(i in 1:length(tmpYear)){
      
      index = data$site2==j&data$year2==tmpYear[i]
      tmp.obs=as.matrix(y.obs[index])
      tmp.sim=as.matrix(y.sim[,index])
      test.obs=model.fun(tmp.obs,na.rm=TRUE)
      test.sim=apply(tmp.sim,1,model.fun,na.rm=TRUE)
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      index = data$siteYearIndex_tot_observed[data$site_tot_observed==j&data$year_tot_observed==tmpYear[i]] 
      p.pop[,index] = p.test.calc
      
    }
  }
  
  p.test.mat = apply(p.pop,2,mean,na.rm=TRUE)
  
  p.test.list = list()
  for(j in 1:20){
    tmpIndex = data$siteYearIndex_tot_observed[data$site_tot_observed==j]
    p.test.list[[j]] = p.test.mat[tmpIndex]
  }
  
  return(p.test.list)
}

sims=y_dam.sim
df=data$y_dam

dam.mean=f(y.sim=sims,y.obs=df,model.fun=mean)

p.chi.mat <- f.convert(p.chi.list)
dam.mean.mat<-f.convert(dam.mean)

pdf(file=paste0(outputDirectory,"pvals-damagedFruitsPerPlant.pdf"),height=3,width=6)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

f.plot(diagnostic.test.mat = dam.mean.mat,years=2013:2020)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat,years=2013:2020)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)

dev.off()


# ---
# Graphical checks: Damaged fruits per plant ----
# ---

y.obs=data$y_dam
n.iter=dim(y_dam.sim)[1]
years=2013:2020
n.samples = 50
plot.yr = 2013

pdf(file=paste0(outputDirectory,"ppc-damagedFruitsPerPlant.pdf"),height=4,width=6)

par(mfrow = c(2,4),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(2,1,0))

for(i in (1:length(years))[ years %in% plot.yr]){
  if(colSums(p.chi.mat>.95,na.rm=TRUE)[i]==0){
    
  } else {
    
    iter.ind = sample(1:n.iter,n.samples)
    
    # get the sites with poor fit
    site.index<-(1:20)[(p.chi.mat>.95&!is.na(p.chi.mat))[,i]]
    site.index<-c(site.index,(1:20)[(p.chi.mat<.05&!is.na(p.chi.mat))[,i]])
    # add two random sites
    add.index <- c(sample((1:20)[!(1:20%in%site.index)],4))
    
    for(j in c(site.index,add.index)){
      
      index=data$site2==j&data$year2==i
      
      p.obs=y.obs[index]
      if(length(p.obs)>0){
      p.sim=as.matrix(y_dam.sim[,index])
      p.sim=as.matrix(p.sim[iter.ind,])
      
      if(is.matrix(p.sim)&dim(p.sim)[2]>1){
        
        list.dens=apply(p.sim,1,density,na.rm=TRUE)
        all.max.y=max(unlist(lapply(list.dens, "[", "y")))
        all.max.x=max(unlist(lapply(list.dens, "[", "x")))
        
        m.max=max(p.obs,na.rm=TRUE)
        m.min=min(p.obs,na.rm=TRUE)
        dens.obs = density(p.obs,from=m.min,to=m.max,na.rm=TRUE)
        
        all.max.y = max(all.max.y,max(dens.obs$y))
        all.max.x = max(all.max.x,max(dens.obs$x))
        
        plot(NA,NA,
             ylim=c(0,all.max.y),xlim=c(0,all.max.x),
             xaxt='n',xlab='',ylab='',yaxt='n')
        for(h in 1:n.samples){
          m.max=max(p.sim[h,],na.rm=TRUE)
          m.min=min(p.sim[h,],na.rm=TRUE)
          
          dens.x = density(p.sim[h,],from=m.min,to=m.max,na.rm=TRUE)
          lines(x=dens.x$x,y=dens.x$y,lwd=0.25,
                col=ifelse(j %in% site.index,'orange','gray75'))
        }
        
        p = data$y_dam[index]
        
        lines(x=dens.obs$x,y=dens.obs$y)
        
        if(length(unique(p))==1){ abline(v=unique(p),lty='dotted',col='black')}
        
      } else if(is.matrix(p.sim)&dim(p.sim)[2]==1){
        # lines below are for the case where there is only 1 observation
        
        p.sim = p.sim
        all.max.x= max(density(p.sim,na.rm=TRUE)$x)
        all.max.y= max(density(p.sim,na.rm=TRUE)$y)
        
        plot(NA,NA,
             ylim=c(0,all.max.y),xlim=c(0,all.max.x),
             xaxt='n',xlab='',ylab='',yaxt='n')
        
        # plot observed value
        p = data$y_dam[index]
        abline(v=p,lty='dotted',col='black')
        
        # overlay simulated values
        tb = table(p.sim)
        tb.x = as.numeric(names(tb))
        segments(x0=tb.x,y0=0,y1=tb/sum(tb),
                 col=ifelse(j %in% site.index,'orange','gray75'))
        
      }
      
      axis(2,cex=.25,tick=FALSE,line=-1)
      axis(1,cex=.25,tick=FALSE,line=-1)
      
      legend("topright",paste0(siteAbiotic$site[j],': ',years[i],"\n n=",length(p)),bty='n')
      
    }
    }
    mtext(paste0("Counts of damaged fruits per plant"), side = 1,           
          outer = TRUE, line = .75)
    mtext("Density", side = 2, outer = TRUE, line = 1.2)
  }
}
dev.off()
