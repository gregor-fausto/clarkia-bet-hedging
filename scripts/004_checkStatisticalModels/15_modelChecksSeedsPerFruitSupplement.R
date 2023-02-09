####
####
# Script for model checks for seeds per fruit
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

tmpDataDirectory = "outputs/001_prepareDataForModels/"
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
tmpDataDirectory <- paste0(tmpDataDirectory,list.files(tmpDataDirectory))
data <- readRDS(tmpDataDirectory[[grep("seedsPerFruit.RDS",tmpDataDirectory)]])

# - +Read in MCMC samples ----
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedsPerFruitPosteriorSamples.RDS",mcmcSampleDirectory)]])

# - +Read in site names ----
siteAbiotic <- read.csv("data/siteAbioticData.csv",header=TRUE)
# create variable that arranges populations by easting
siteIndex <- order(siteAbiotic$easting,decreasing=FALSE)
siteNames = unique(siteAbiotic$site)[siteIndex]

# ---
# ---
# Graphical Posterior Predictive Checks
# ---
# ---

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ---
# Test statistics: Seeds per undamaged fruit ----
# ---
y.sim=MCMCchains(mcmcSamples, params = "y_sd.sim")
y.obs=data$sdno
n.iter=dim(y.sim)[1]
years=2006:2020
sample.index = sample(1:n.iter,5000)
y.sim=y.sim[sample.index,]

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.sd.obs"))[sample.index,]
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sd.sim"))[sample.index,]

n.iter=dim(y.sim)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=data$n_siteYearIndex_und)
for(j in 1:20){
  
  tmpYear = data$year_und_observed[data$site_und_observed==j]
  
  for(i in 1:length(tmpYear)){
    index = data$site3==j&data$year3==tmpYear[i]
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
    index = data$siteYearIndex_und_observed[data$site_und_observed==j&data$year_und_observed==tmpYear[i]] 
    p.pop[,index] = p.chi2.calc
  }
}

p.chi.mat = apply(p.pop,2,mean,na.rm=TRUE)

p.chi.list = list()
for(j in 1:20){
  tmpIndex = data$siteYearIndex_und_observed[data$site_und_observed==j]
  p.chi.list[[j]] = p.chi.mat[tmpIndex]
}


f=function(y.sim=chains,y.obs=data,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=data$n_siteYearIndex_und)
  for(j in 1:20){
    
    tmpYear = data$year_und_observed[data$site_und_observed==j]
    
    for(i in 1:length(tmpYear)){
      
      index = data$site3==j&data$year3==tmpYear[i]
      tmp.chi2.obs=chi2.obs[,index]
      tmp.chi2.sim=chi2.sim[,index]
      
      tmp.obs=as.matrix(y.obs[index])
      test.obs=model.fun(tmp.obs,na.rm=TRUE)
      tmp.sim=as.matrix(y.sim[,index])
      test.sim=apply(tmp.sim,1,model.fun,na.rm=TRUE)
      
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      index = data$siteYearIndex_und_observed[data$site_und_observed==j&data$year_und_observed==tmpYear[i]] 
      p.pop[,index] = p.test.calc
      
    }
  }
  p.test.mat = apply(p.pop,2,mean,na.rm=TRUE)
  
  p.test.list = list()
  for(j in 1:20){
    tmpIndex = data$siteYearIndex_und_observed[data$site_und_observed==j]
    p.test.list[[j]] = p.test.mat[tmpIndex]
  }
  
  return(p.test.list)
}

sims=y.sim
df=data$sdno

#sdno.min=f(y.sim=sims,y.obs=df,model.fun=min)
#sdno.max=f(y.sim=sims,y.obs=df,model.fun=max)
sdno.mean=f(y.sim=sims,y.obs=df,model.fun=mean)
#sdno.sd=f(y.sim=sims,y.obs=df,model.fun=sd)

# convert lists to matrix

f.convert = function(test.list){
  out.list <- list()
  for(i in 1:20){
    tmpYear = data$year_und_observed[data$site_und_observed==i]
    
    out.list[[i]] <- ifelse((1:15 %in% tmpYear), test.list[[i]], 
                            ifelse(!(1:15 %in% tmpYear), NA, 0))
  }
  out.mat <- do.call(rbind,out.list)
  return(out.mat)
}

# return matrix objects
p.chi.sdno.mat<-f.convert(p.chi.list)
# sdno.min.mat<-f.convert(sdno.min)
# sdno.max.mat<-f.convert(sdno.max)
sdno.mean.mat<-f.convert(sdno.mean)
# sdno.sd.mat<-f.convert(sdno.sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(length(years))


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


pdf(file=paste0(outputDirectory,"pvals-seedsPerUndamagedFruit.pdf"),height=4,width=6)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

f.plot(diagnostic.test.mat = sdno.mean.mat,years=2006:2020)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.sdno.mat,years=2006:2020)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)

dev.off()


# ---
# Test statistics: Seeds per damaged fruit ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "y_sd_dam.sim")
y.obs=data$sdno_dam
n.iter=dim(y.sim)[1]
years=2013:2020
sample.index = sample(1:n.iter,5000)
y.sim=y.sim[sample.index,]

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.sd_dam.obs"))[sample.index,]
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sd_dam.sim"))[sample.index,]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=data$n_siteYearIndex_dam)
for(j in 1:20){
  
  tmpYear = data$year_dam_observed[data$site_dam_observed==j]
  
  for(i in 1:length(tmpYear)){
    index = data$site4==j&data$year4==tmpYear[i]
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
    index = data$siteYearIndex_dam_observed[data$site_dam_observed==j&data$year_dam_observed==tmpYear[i]] 
    p.pop[,index] = p.chi2.calc
  }
}

p.chi.mat = apply(p.pop,2,mean,na.rm=TRUE)

p.chi.list = list()
for(j in 1:20){
  tmpIndex = data$siteYearIndex_dam_observed[data$site_dam_observed==j]
  p.chi.list[[j]] = p.chi.mat[tmpIndex]
}


f=function(y.sim=chains,y.obs=data,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=data$n_siteYearIndex_dam)
  for(j in 1:20){
    
    tmpYear = data$year_dam_observed[data$site_dam_observed==j]
    
    for(i in 1:length(tmpYear)){
      
      index = data$site4==j&data$year4==tmpYear[i]
      tmp.chi2.obs=chi2.obs[,index]
      tmp.chi2.sim=chi2.sim[,index]
      
      tmp.obs=as.matrix(y.obs[index])
      test.obs=model.fun(tmp.obs,na.rm=TRUE)
      tmp.sim=as.matrix(y.sim[,index])
      test.sim=apply(tmp.sim,1,model.fun,na.rm=TRUE)
      
      p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
      index = data$siteYearIndex_dam_observed[data$site_dam_observed==j&data$year_dam_observed==tmpYear[i]] 
      p.pop[,index] = p.test.calc
      
    }
  }
  p.test.mat = apply(p.pop,2,mean,na.rm=TRUE)
  
  p.test.list = list()
  for(j in 1:20){
    tmpIndex = data$siteYearIndex_dam_observed[data$site_dam_observed==j]
    p.test.list[[j]] = p.test.mat[tmpIndex]
  }
  
  return(p.test.list)
}

sims=y.sim
df=data$sdno_dam

# sdno.min=f(y.sim=sims,y.obs=df,model.fun=min)
# sdno.max=f(y.sim=sims,y.obs=df,model.fun=max)
sdno.mean=f(y.sim=sims,y.obs=df,model.fun=mean)
# sdno.sd=f(y.sim=sims,y.obs=df,model.fun=sd)

# convert lists to matrix

f.convert = function(test.list){
  out.list <- list()
  for(i in 1:20){
    tmpYear = data$year_dam_observed[data$site_dam_observed==i]
    
    out.list[[i]] <- ifelse((1:8 %in% tmpYear), test.list[[i]], 
                            ifelse(!(1:8 %in% tmpYear), NA, 0))
  }
  out.mat <- do.call(rbind,out.list)
  return(out.mat)
}

# return matrix objects
p.chi.sdnodam.mat<-f.convert(p.chi.list)
# sdno.min.mat<-f.convert(sdno.min)
# sdno.max.mat<-f.convert(sdno.max)
sdnodam.mean.mat<-f.convert(sdno.mean)
# sdno.sd.mat<-f.convert(sdno.sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(length(years))

pdf(file=paste0(outputDirectory,"pvals-seedsPerDamagedFruit.pdf"),height=4,width=6)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

f.plot(diagnostic.test.mat = sdnodam.mean.mat,years=2013:2020)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.sdnodam.mat,years=2013:2020)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)

dev.off()


# ---
# Graphical checks: Seeds per undamaged fruit ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "y_sd.sim")
y.obs=data$sdno
#n.iter=dim(y.sim)[1]
years=2006:2020
n.samples = 50
#y.sim = y.sim[sample(1:n.iter,n.samples),]
plot.yr = c(2010)

pdf(file=paste0(outputDirectory,"ppc-seedsPerUndamagedFruit-2010.pdf"),height=2.5,width=4)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(2,1,0))

for(i in (1:length(years))[ years %in% plot.yr]){
  if(colSums(p.chi.sdno.mat<.05,na.rm=TRUE)[i]==0){
    
  } else {
    
    iter.ind = sample(1:n.iter,n.samples)
    
    # get the sites with poor fit
    site.index<-(1:20)[(p.chi.sdno.mat<.05&!is.na(p.chi.sdno.mat))[,i]]
    # add two random sites
    add.index <- c(sample((1:20)[!(1:20%in%site.index)],1))
    for(j in c(site.index,add.index)){
      
      index=data$site3==j&data$year3==i
      
      p.obs=y.obs[index]
      p.sim=(y.sim[,index])
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
        
        p = data$sdno[index]
        
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
        p = data$sdno[index]
        abline(v=p,lty='dotted',col='black')
        
        # overlay simulated values
        tb = table(p.sim)
        tb.x = as.numeric(names(tb))
        segments(x0=tb.x,y0=0,y1=tb/sum(tb),
                 col=ifelse(j %in% site.index,'orange','gray75'))
        
      }
      
      axis(2,cex=.25,tick=FALSE,line=-1)
      axis(1,cex=.25,tick=FALSE,line=-1)
      
      legend("topright",paste0(siteNames[j],': ',years[i],"\n n=",length(p)),bty='n')
      
    } 
  }
  
  mtext(paste0("Counts of seeds per undamaged frut (",years[i],")"),
        side = 1, outer = TRUE, line = .75)
  mtext("Density", side = 2, outer = TRUE, line = 1.2)
} 
dev.off()

y.sim=MCMCchains(mcmcSamples, params = "y_sd.sim")
y.obs=data$sdno
#n.iter=dim(y.sim)[1]
years=2006:2020
n.samples = 50
plot.yr = c(2014)

pdf(file=paste0(outputDirectory,"ppc-seedsPerUndamagedFruit-2014.pdf"),height=2.5,width=4)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(2,1,0))

for(i in (1:length(years))[ years %in% plot.yr]){
  if(colSums(p.chi.sdno.mat>.95,na.rm=TRUE)[i]==0){
    
  } else {
    
    iter.ind = sample(1:n.iter,n.samples)
    
    # get the sites with poor fit
    site.index<-(1:20)[(p.chi.sdno.mat>.95&!is.na(p.chi.sdno.mat))[,i]]
    # add two random sites
    add.index <- c(sample((1:20)[!(1:20%in%site.index)],1))
    for(j in c(site.index,add.index)){
      
      index=data$site3==j&data$year3==i
      
      p.obs=y.obs[index]
      p.sim=as.matrix(y.sim[,index])
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
        
        p = data$sdno[index]
        
        lines(x=dens.obs$x,y=dens.obs$y)
        
        if(length(unique(p))==1){ abline(v=unique(p),lty='dotted',col='black')}
        
      } else if(is.matrix(p.sim)&dim(p.sim)[2]==1){
        # lines below are for the case where there is only 1 observation
        
        p.sim = p.sim
        all.max.x= max(density(p.sim,na.rm=TRUE)$x)
        all.max.y= max(density(p.sim,na.rm=TRUE)$y)
        
        plot(NA,NA,
             ylim=c(0,all.max.y*1.5),xlim=c(0,all.max.x),
             xaxt='n',xlab='',ylab='',yaxt='n')
        
        # plot observed value
        p = data$sdno[index]
        abline(v=p,lty='dotted',col='black')
        
        # overlay simulated values
        tb = table(p.sim)
        tb.x = as.numeric(names(tb))
        segments(x0=tb.x,y0=0,y1=tb/sum(tb),
                 col=ifelse(j %in% site.index,'orange','gray75'))
        
      }
      
      axis(2,cex=.25,tick=FALSE,line=-1)
      axis(1,cex=.25,tick=FALSE,line=-1)
      
      legend("topright",paste0(siteNames[j],': ',years[i],"\n n=",length(p)),bty='n')
      
    } 
  }
  
  mtext(paste0("Counts of seeds per undamaged frut (",years[i],")"),
        side = 1, outer = TRUE, line = .75)
  mtext("Density", side = 2, outer = TRUE, line = 1.2)
} 
dev.off()

# ---
# Graphical checks: Seeds per damaged fruit ----
# ---

y_dam.sim=MCMCchains(mcmcSamples, params = "y_sd_dam.sim")
y_dam.obs=data$sdno_dam
n.iter=dim(y_dam.sim)[1]
years=2013:2020
n.samples = 50
plot.yr = 2013

pdf(file=paste0(outputDirectory,"ppc-seedsPerDamagedFruit.pdf"),height=4,width=6)


par(mfrow = c(2,5),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(2,1,0))

for(i in (1:length(years))[ years %in% plot.yr]){
  if(colSums(sdnodam.mean.mat>.95,na.rm=TRUE)[i]==0){
    
  } else {
    
    iter.ind = sample(1:n.iter,n.samples)
    
    # get the sites with poor fit
    site.index<-(1:20)[(sdnodam.mean.mat>.95&!is.na(sdnodam.mean.mat))[,i]]
    # add two random sites
    add.index <- c(sample((1:20)[!(1:20%in%site.index)],5))
    for(j in c(site.index,add.index)){
      
      index=data$site4==j&data$year4==i
      
      p.obs=y_dam.obs[index]
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
          
          p = data$sdno_dam[index]
          
          lines(x=dens.obs$x,y=dens.obs$y)
          
          if(length(unique(p))==1){ abline(v=unique(p),lty='dotted',col='black')}
          
        } else if(is.matrix(p.sim)&dim(p.sim)[2]==1){
          # lines below are for the case where there is only 1 observation
          
          p.sim = p.sim
          all.max.x= max(density(p.sim,na.rm=TRUE)$x)
          all.max.y= max(density(p.sim,na.rm=TRUE)$y)
          
          plot(NA,NA,
               ylim=c(0,all.max.y*1.5),xlim=c(0,all.max.x),
               xaxt='n',xlab='',ylab='',yaxt='n')
          
          # plot observed value
          p = data$sdno_dam[index]
          abline(v=p,lty='dotted',col='black')
          
          # overlay simulated values
          tb = table(p.sim)
          tb.x = as.numeric(names(tb))
          segments(x0=tb.x,y0=0,y1=tb/sum(tb),
                   col=ifelse(j %in% site.index,'orange','gray75'))        
        }
        
        axis(2,cex=.25,tick=FALSE,line=-1)
        axis(1,cex=.25,tick=FALSE,line=-1)
        
        legend("topright",paste0(siteNames[j],': ',years[i],"\n n=",length(p)),bty='n')
      } 
      
    }
    mtext(paste0("Counts of seeds per damaged fruit (",years[i],")"), side = 1, outer = TRUE, line = .75)
    mtext("Density", side = 2, outer = TRUE, line = 1.2)
  }
}
dev.off()


