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
fullDataDirectory <- paste0(fullDataDirectory,list.files(fullDataDirectory))
data <- readRDS(fullDataDirectory[[grep("seedlingFruitingPlantCountsPermanentPlots.RDS",fullDataDirectory)]])

# - +Read in MCMC samples ----
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedlingSurvivalPosteriorSamples.RDS",mcmcSampleDirectory)]])

# - +Read in site names ----
siteAbiotic <- read.csv("~/Dropbox/clarkia-demography-projects/data/siteAbioticData.csv",header=TRUE) 

# create variable that arranges populations by easting
siteIndex <- order(siteAbiotic$easting,decreasing=FALSE)
siteNames = unique(siteAbiotic$site)[siteIndex]


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

n.iter = dim(chi2.obs)[1]
n.iter = length(subsample)


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
      tmp.chi2.obs=chi2.obs[,index]
      tmp.chi2.sim=chi2.sim[,index]
      
      tmp.obs=y.obs[index]/n.obs[index]
      test.obs=model.fun(tmp.obs,na.rm=TRUE)
      
      if(is.matrix(tmp.chi2.obs)){
        tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
        test.sim=apply(tmp.sim,1,model.fun,na.rm=TRUE)
      } else if(is.vector(tmp.chi2.obs)){
        tmp.sim=y.sim[,index]/n.obs[index]
        test.sim=tmp.sim
        #test.sim=model.fun(tmp.sim,na.rm=TRUE)
      }
      
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

# seedlings.min=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=min)
# seedlings.max=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=max)
seedlings.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)
# seedlings.sd=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=sd)

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
# seedlings.min.mat<-f.convert(seedlings.min)
# seedlings.max.mat<-f.convert(seedlings.max)
seedlings.mean.mat<-f.convert(seedlings.mean)
# seedlings.sd.mat<-f.convert(seedlings.sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(15)



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


pdf(file=paste0(outputDirectory,"pvals-seedlingSurvivalFruiting.pdf"),height=4,width=6)

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
# 
# par(mfrow=c(4,5),mar=c(1,1,1,1),oma=c(1,1,1,1))
# for(i in 1:20){
#   plot(1:15,seedlings.mean.mat[i,],type='b',ylim=c(0,1),
#        pch=ifelse(seedlings.mean.mat[i,]<.05|seedlings.mean.mat[i,]>.95,21,16))
# }
# 
# par(mfrow=c(3,5),mar=c(1,1,1,1),oma=c(1,1,1,1))
# for(i in 1:15){
#   plot(1:20,seedlings.mean.mat[,i][siteIndex],type='b',ylim=c(0,1),
#        pch=ifelse(seedlings.mean.mat[,i][siteIndex]<.05|seedlings.mean.mat[,i][siteIndex]>.95,21,16))
# }
# 
# var.mat = matrix(NA,nrow=20,ncol=15)
# for(i in 1:20){
#   for(j in 1:15){
#     index=data$site==i&data$year==j
#     n_sdlg=data$seedlingNumber[index]
#     n_frts=data$fruitplNumber[index]
#     #p=n_frts/n_sdlg
#     var.mat[i,j]=var(p,na.rm=TRUE)
#     var.mat[i,j]=sum(n_sdlg,na.rm=TRUE)
#   }
# }
# plot(var.mat,seedlings.mean.mat)
# 
# par(mfrow=c(4,5),mar=c(1,1,1,1),oma=c(1,1,1,1))
# for(i in 1:20){
#   plot(var.mat[i,],seedlings.mean.mat[i,],type='p',ylim=c(0,1),
#        pch=ifelse(seedlings.mean.mat[i,]<.05|seedlings.mean.mat[i,]>.95,21,16))
# }


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

par(mfrow = c(3,5),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(2,1,0))
  
for(j in 1){
    
  for(i in 1:length(years)){
    
    iter.ind = sample(1:n.sim,n.samples)
    tmpSite = data$site_observed[data$year_observed==i]
    
    j.tmp = siteIndex[j]
    if(j.tmp %in% tmpSite){
      
      index=data$site==j.tmp&data$year==i
      
      p.obs=y.obs[index]/n.obs[index]
     # p.obs=sweep(y.obs[index], 2, n.obs[index], FUN = '/')
      
  #    p.sim=as.matrix(y.sim[,index])
      p.sim=sweep(t(y.sim[,index]), 2, n.obs[index], FUN = '/')
      p.sim=as.matrix(t(p.sim)[iter.ind,])
      
      if(is.matrix(p.sim)&dim(p.sim)[2]>1){
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
                col=ifelse(p.chi.mat[j,i]>.95|p.chi.mat[j,i]<.05,
                             'orange','gray75'))
        }
        
        p = data$fruitplNumber[index]
        m.max=max(p,na.rm=TRUE)
        m.min=min(p,na.rm=TRUE)
        
        dens.x = density(p,from=m.min,to=m.max,na.rm=TRUE)
        
        lines(x=dens.x$x,y=dens.x$y)
        
        if(length(unique(p))==1){ abline(v=unique(p),lty='dotted',col='black')}
        
      } else if(is.matrix(p.sim)&dim(p.sim)[2]==1){
        # lines below are for the case where there is only 1 observation
        
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
        segments(x0=tb.x,y0=0,y1=tb/sum(tb),col='orange')
        
      }
      
      ifelse(i%in%c(1),axis(2L),NA)
      
      axis(2,cex=.25,tick=FALSE,line=-1)
      axis(1,cex=.25,tick=FALSE,line=-1)
      
      legend("topright"
             ,paste0(siteNames[j],": ",years[i],"\n n=",length(p)),
             bty='n')
    } else {
      plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
    }
  }
  
  mtext(paste0("Probability of seedling survival"), 
        side = 1, outer = TRUE, line = .75)
  mtext("Density", side = 2, outer = TRUE, line = 1.2)
}

dev.off()



# -------
# Chi-2 ----
# -------

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.obs"))[subsample,]
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sim"))[subsample,]
# calculations are rowwise
fit.obs=apply(chi2.obs,1,sum)
fit.sim=apply(chi2.sim,1,sum)
p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
mean(p.chi2.calc)

p.pop = matrix(NA,nrow=dim(chi2.obs)[1],ncol=20)
for(i in 1:20){
  tmp.chi2.obs=chi2.obs[,data$site==i]
  tmp.chi2.sim=chi2.sim[,data$site==i]
  fit.obs=apply(tmp.chi2.obs,1,sum)
  fit.sim=apply(tmp.chi2.sim,1,sum)
  p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
  p.pop[,i] = p.chi2.calc
}
apply(p.pop,2,mean)


pdf(file=paste0(outputDirectory,"population-seedlingSurvivalFruiting.pdf"),height=6,width=8)

par(mfrow=c(1,2))
pop.sample = 1:20
plot(y=(pop.sample),x=rev(apply(p.pop,2,mean)),
     xlim=c(0,1),pch=16,
     ylab="Population",xlab="p-Value (chi-squared)",
     main="Fruiting plant counts",    
     axes=FALSE,frame=FALSE,xaxt='n',yaxt='n')
abline(v=c(.1,.9),lty='dotted')
abline(v=c(.2,.8),lty='dotted',col='gray')

axis(2, (1:20),
     labels = rev(siteNames), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(1,  seq(0,1,by=.2), col.ticks = 1)

# population-wide, this seems to pass checks okay

# n.iter = dim(chi2.obs)[1]
# 
# p.chi.list = list()
# p.pop = matrix(NA,nrow=n.iter,ncol=length(years))
# for(j in 1:20){
#   for(i in 1:length(years)){
#     tmp.chi2.obs=chi2.obs[,data$site==j& data$year==i]
#     tmp.chi2.sim=chi2.sim[,data$site==j& data$year==i]
#     fit.obs=apply(tmp.chi2.obs,1,sum)
#     fit.sim=apply(tmp.chi2.sim,1,sum)
#     #   exclude=fit.obs!=fit.sim
#     #   p.chi2.calc=ifelse(fit.sim[exclude]-fit.obs[exclude]>=0,1,0)
#     p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
#     p.pop[,i] = p.chi2.calc
#   }
#   p.chi.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
# }
# p.chi.mat=do.call(rbind,p.chi.list)

time.sample = 1:length(years)+2005
plot(time.sample,p.chi.mat[1,],
     ylim=c(0,1),pch=16,xlim=c(2006,2019),
     xlab="Year",ylab="p-Value (chi-squared)",
     main="Fruiting plant counts",type='n')
for(i in 1:20){
  points(time.sample+rnorm(1,0,sd=.05),
         p.chi.mat[i,],pch=19,cex=.5,col='black')
}
abline(h=c(.1,.9),lty='dotted')
abline(h=c(.2,.8),lty='dotted',col='gray')
dev.off()

# per year*population shows some values that are not well modeled this way

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

time.sample = 1:length(years)+2005
for(i in 1:20){
  
  tot.seeds=c()
  for(j in 1:length(years)){
    tot.seeds[j]=sum(data$seedlingNumber[data$site==i&data$year==j],na.rm=TRUE)
  }
  
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(2006,2019),
       ylab='',xlab='',xaxt='n',yaxt='n')
  points(time.sample+rnorm(1,0,sd=.05),
         p.chi.mat[i,],pch=19,cex=1,
         col="black")
  
  abline(h=c(.1,.9),lty='dotted')
  abline(h=c(.2,.8),lty='dotted',col='gray')
  
  text(x=2005.5,y=.04,siteNames[i],pos=4)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
  ifelse(i%in%c(5), legend(x = 15, y = 1,
                           col = c('gray','orange'),
                           lty = c(1,1),
                           legend = c("Persistence only","Persistence & viability"),
                           cex=.55,
                           box.lty=0), NA)
}
mtext("Year", side = 1, outer = TRUE, line = 2.2)
mtext("Bayesian p-value (chi-2)", side = 2, outer = TRUE, line = 2.2)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1,
    mgp=c(3,0,0))

for(i in 1:20){
  
  tot.fruitpl=c()
  for(j in 1:length(years)){
    tot.fruitpl[j]=sum(data$fruitplNumber[data$site==i&data$year==j],na.rm=TRUE)
  }
  
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(0,max(tot.fruitpl)),
       ylab='',xlab='',yaxt='n',xaxt='n')
  points(tot.fruitpl,
         p.chi.mat[i,],pch=19,cex=1,
         col=ifelse(p.chi.mat[i,]>.9|p.chi.mat[i,]<.1,"purple","black"))
  
  abline(h=c(.1,.9),lty='dotted')
  
  text(x=5,y=.04,siteNames[i],pos=4)
  axis(1,cex.axis=.75,tick=FALSE)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Total number of fruiting plants", side = 1, outer = TRUE, line = 2.2)
mtext("Bayesian p-value (chi-2)", side = 2, outer = TRUE, line = 2.2)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1,
    mgp=c(3,0,0))

for(i in 1:20){
  
  var.fruitpl=c()
  for(j in 1:length(years)){
    var.fruitpl[j]=var(data$fruitplNumber[data$site==i&data$year==j],na.rm=TRUE)
  }
  
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(0,max(var.fruitpl,na.rm=TRUE)),
       ylab='',xlab='',yaxt='n',xaxt='n')
  points(var.fruitpl,
         p.chi.mat[i,],pch=19,cex=1,
         col=ifelse(p.chi.mat[i,]>.9|p.chi.mat[i,]<.1,"purple","black"))
  
  abline(h=c(.1,.9),lty='dotted')
  
  text(x=5,y=.04,siteNames[i],pos=4)
  axis(1,cex.axis=.75,tick=FALSE)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Variance in fruiting plant number", side = 1, outer = TRUE, line = 2.2)
mtext("Bayesian p-value (chi-2)", side = 2, outer = TRUE, line = 2.2)

# Years in which almost no plants survived, this model does not capture that


pdf(file=paste0(outputDirectory,"populationYear-seedlingSurvivalFruiting.pdf"),height=6,width=6)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1,
    mgp=c(3,0,0))

for(i in 1:20){
  
  tot.seeds=c()
  for(j in 1:length(years)){
    tot.seeds[j]=sum(data$seedlingNumber[data$site==i&data$year==j],na.rm=TRUE)
  }
  
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(0,max(tot.seeds,na.rm=TRUE)),
       ylab='',xlab='',yaxt='n',xaxt='n')
  points(tot.seeds,
         p.chi.mat[i,],pch=19,cex=1,
         col=ifelse(p.chi.mat[i,]>.9|p.chi.mat[i,]<.1,"purple","black"))
  
  abline(h=c(.1,.9),lty='dotted')
  
  text(x=5,y=.04,siteNames[i],pos=4)
  axis(1,cex.axis=.75,tick=FALSE)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Total number of seedlings", side = 1, outer = TRUE, line = 2.2)
mtext("Bayesian p-value (chi-2)", side = 2, outer = TRUE, line = 2.2)
#mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1,
    mgp=c(3,0,0))

for(i in 1:20){
  
  var.fruitpl=c()
  for(j in 1:length(years)){
    var.fruitpl[j]=sum(ifelse(data$seedlingNumber[data$site==i&data$year==j]>0,1,0))
  }
  
  plot(NA,NA,
       ylim=c(0,1),pch=16,xlim=c(0,max(var.fruitpl)),
       ylab='',xlab='',yaxt='n',xaxt='n')
  points(var.fruitpl,
         p.chi.mat[i,],pch=19,cex=1,
         col=ifelse(p.chi.mat[i,]>.9|p.chi.mat[i,]<.1,"purple","black"))
  
  abline(h=c(.1,.9),lty='dotted')
  
  text(x=5,y=.04,siteNames[i],pos=4)
  axis(1,cex.axis=.75,tick=FALSE)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}
mtext("Plots with nonzero numbers of seedlings", side = 1, outer = TRUE, line = 2.2)
mtext("Bayesian p-value (chi-2)", side = 2, outer = TRUE, line = 2.2)
#mtext("Population*year-level", side = 1, outer = TRUE, line = 2.5,adj=-.05,cex=1.25)



# years with low number of seedlings (trials) and fruiting plants (successes)
# pull the mean towards the population-level average
# due to pooling


dev.off()



# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Graphical Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------



