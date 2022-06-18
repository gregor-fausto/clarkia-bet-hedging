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
data <- readRDS(modelDataDirectory[[grep("seedlingSurvivalData.RDS",modelDataDirectory)]])

# - +Read in MCMC samples ----
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedlingSurvivalSamplesSimple.RDS",mcmcSampleDirectory)]])
mcmcSamplesPlots <- readRDS(mcmcSampleDirectory[[grep("seedlingSurvivalSamplesSimplePlots.RDS",mcmcSampleDirectory)]])

# - +Read in site names ----
siteAbiotic <- read.csv("data/siteAbiotic.csv",header=TRUE) 

siteNames <- siteAbiotic$site

# create variable that arranges populations by easting
siteIndex <- order(siteAbiotic$easting,decreasing=FALSE)
siteNames = unique(siteAbiotic$site)[siteIndex]

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Observations vs. model predictions: temporal only
# -------------------------------------------------------------------
# -------------------------------------------------------------------

mu0 <- MCMCchains(mcmcSamples,params="mu0")
mu0.summary <- apply(boot::inv.logit(mu0),2,quantile,c(.025,.5,.975))

mu_year <- MCMCchains(mcmcSamples,params="mu")
mu.year.summary <- apply(boot::inv.logit(mu_year),2,quantile,c(.025,.5,.975))

tmp.list <- list()
par(mfrow=c(1,4),mar=c(2,2,1,1))
for(i in 1:4){
  
  plot(NA,xlim=c(0,15),ylim=c(0,1))
  
  for(j in 1:15){
    index= data$site==i&data$year==j
    tmp=(data$fruitplNumber/data$seedlingNumber)[index]
    
    points(rep(j,length(tmp)),tmp,col='gray',pch=1,cex=1)
    points(j,mean(tmp),pch=4)
    if(length(mu.year.summary[2,unique(data$siteYearIndex[index])])>0){
      points(j,mu.year.summary[2,unique(data$siteYearIndex[index])],
             pch=16,cex=1,col='red')}
  }
}

par(mfrow=c(1,4),mar=c(2,2,1,1))
for(i in 1:4){
  
  i.tmp = siteIndex[i]
  
  plot(NA,xlim=c(0,15),ylim=c(0,1))
  abline(h=mu0.summary[2,i.tmp],lty='dotted')
  
  
  
  for(j in 1:15){
    index= data$site==i.tmp&data$year==j
    tmp=(data$fruitplNumber/data$seedlingNumber)[index]
    
    if(length(mu.year.summary[2,unique(data$siteYearIndex[index])])>0){
      points(j,mu.year.summary[2,unique(data$siteYearIndex[index])],
             pch=16,cex=1,col='black')
      segments(j,y0=mu.year.summary[1,unique(data$siteYearIndex[index])],
               y1=mu.year.summary[3,unique(data$siteYearIndex[index])])
    }
  }
}

par(mfrow=c(1,1))
sigma0 <- MCMCchains(mcmcSamples,params="sigma0")
sigma0.sum<-apply(sigma0,2,quantile,c(0.025,.5,.975))
plot(siteIndex,sigma0.sum[2,siteIndex],ylim=c(0,3))
segments(siteIndex,y0=sigma0.sum[1,siteIndex],y1=sigma0.sum[3,siteIndex])

plot(siteAbiotic$easting,sigma0.sum[2,],ylim=c(0,3))
segments(siteAbiotic$easting,y0=sigma0.sum[1,],y1=sigma0.sum[3,])


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Observations vs. model predictions: temporal+plot 
# -------------------------------------------------------------------
# -------------------------------------------------------------------

mu0 <- MCMCchains(mcmcSamplesPlots,params="mu0")
mu0.summary <- apply(boot::inv.logit(mu0),2,quantile,c(.025,.5,.975))

mu_year <- MCMCchains(mcmcSamplesPlots,params="mu")
mu.year.summary <- apply(boot::inv.logit(mu_year),2,quantile,c(.025,.5,.975))

tmp.list <- list()
par(mfrow=c(4,5),mar=c(2,2,1,1))
for(i in 1:20){
  
  plot(NA,xlim=c(0,15),ylim=c(0,1))
  
  for(j in 1:15){
    index= data$site==i&data$year==j
    tmp=(data$fruitplNumber/data$seedlingNumber)[index]
    
    points(rep(j,length(tmp)),tmp,col='gray',pch=1,cex=1)
    points(j,mean(tmp),pch=4)
    if(length(mu.year.summary[2,unique(data$siteYearIndex[index])])>0){
      points(j,mu.year.summary[2,unique(data$siteYearIndex[index])],
             pch=16,cex=1,col='red')}
  }
}

# observed vs predicted with individual level observations
par(mfrow=c(4,5),mar=c(2,2,1,1))
for(i in 1:20){ 
  
  plot(NA,xlim=c(0,1),ylim=c(0,1))
  abline(a=0,b=1)
  
  for(j in 1:15){
    index= data$site==i&data$year==j
    tmp=(data$fruitplNumber/data$seedlingNumber)[index]
    
    if(length(mu.year.summary[2,unique(data$siteYearIndex[index])])>0){
      
      est.tmp = mu.year.summary[2,unique(data$siteYearIndex[index])]
      points(x=rep(est.tmp,length(tmp)),y=tmp,col='gray',pch=16,cex=1)

      points(est.tmp,mean(tmp),
             pch=16,cex=1.25,col='black')
      }
  }
}

par(mfrow=c(4,5),mar=c(2,2,1,1))
for(i in 1:20){
  
  i.tmp = siteIndex[i]
  
  plot(NA,xlim=c(0,15),ylim=c(0,1))
  abline(h=mu0.summary[2,i.tmp],lty='dotted')
  
  
  
  for(j in 1:15){
    index= data$site==i.tmp&data$year==j
    tmp=(data$fruitplNumber/data$seedlingNumber)[index]
    
    if(length(mu.year.summary[2,unique(data$siteYearIndex[index])])>0){
      points(j,mu.year.summary[2,unique(data$siteYearIndex[index])],
             pch=16,cex=1,col='black')
      segments(j,y0=mu.year.summary[1,unique(data$siteYearIndex[index])],
               y1=mu.year.summary[3,unique(data$siteYearIndex[index])])
    }
  }
}

par(mfrow=c(1,1))
sigma0 <- MCMCchains(mcmcSamplesPlots,params="sigma0")
sigma0.sum<-apply(sigma0,2,quantile,c(0.025,.5,.975))
plot(siteIndex,sigma0.sum[2,siteIndex],ylim=c(0,3))
segments(siteIndex,y0=sigma0.sum[1,siteIndex],y1=sigma0.sum[3,siteIndex])

plot(siteAbiotic$easting,sigma0.sum[2,],ylim=c(0,3))
segments(siteAbiotic$easting,y0=sigma0.sum[1,],y1=sigma0.sum[3,])
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Estimates vs. MLEs
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# full model: temporal only
mu0 <- MCMCchains(mcmcSamples,params="mu0")
mu0 <- apply(boot::inv.logit(mu0),2,quantile,c(.025,.5,.975))

mu_year <- MCMCchains(mcmcSamples,params="mu")
mu.year.summary <- apply(boot::inv.logit(mu_year),2,quantile,c(.025,.5,.975))

# full model: temporal+spatial
mu0.p  <- MCMCchains(mcmcSamplesPlots,params="mu0")
mu0.p <- apply(boot::inv.logit(mu0.p),2,quantile,c(.025,.5,.975))

mu_year.p  <- MCMCchains(mcmcSamplesPlots,params="mu")
mu.year.p.summary <- apply(boot::inv.logit(mu_year.p),2,quantile,c(.025,.5,.975))

plot(mu0[2,],mu0.p[2,],ylim=c(0,1),xlim=c(0,1))

plot(mu.year.summary[2,],
     mu.year.p.summary[2,])
abline(a=0,b=1)

# full model

# estimates across all sites
y=data$fruitplNumber
n=data$seedlingNumber
site = data$site

year = data$year
siteYearIndex = data$siteYearIndex
d=data.frame(site=site,year=year,siteYearIndex=siteYearIndex,
             n=n,y=y,p=y/n)

d.summary = d %>% 
  dplyr::group_by(siteYearIndex) %>% 
  dplyr::summarise(p.yr=mean(p,na.rm=TRUE),n.yr=sum(n))

df<-d.summary

# observed vs. predicted
par(mfrow=c(1,2))
plot(mu.year.summary[2,],d.summary$p.yr,xlim=c(0,1),ylim=c(0,1),
     pch=NA,xlab="Predicted seedling survival (mean of posterior)",
     ylab="Observed seedling survival (mean of observations)")
abline(a=0,b=1)
points(d.summary$p.yr,mu.year.summary[2,],pch=1)

plot(mu.year.p.summary[2,],d.summary$p.yr,xlim=c(0,1),ylim=c(0,1),
     pch=NA,xlab="Predicted seedling survival (mean of posterior)",
     ylab="Observed seedling survival (mean of observations)")
abline(a=0,b=1)
points(d.summary$p.yr,mu.year.p.summary[2,],pch=1)


par(mfrow=c(1,4),mar=c(4,4,1,1))
for(i in 1:4){
  
  index = data$siteYearIndex[data$site==i]
  index = unique(index)
  
  plot(mu.year.summary[2,index],d.summary$p.yr[index],xlim=c(0,1),ylim=c(0,1),
       pch=NA,xlab="Mean (MLE)",ylab="Posterior median",
       main=siteNames[i])
  abline(a=0,b=1)
  points(mu.year.summary[2,index],d.summary$p.yr[index],pch=16)
  segments(y0=d.summary$p.yr[index],x0=mu.year.summary[1,index],x1=mu.year.summary[3,index])
  
}

mtext(paste0("Predicted seedling survival (mean of posterior)"), side = 1, outer = TRUE, line = 1.5)
mtext("Observed seedling survival (mean of observations)", side = 2, outer = TRUE, line = 2.2)

par(mfrow=c(4,5),mar=c(4,4,1,1))
for(i in 1:20){
  
  index = data$siteYearIndex[data$site==i]
  index = unique(index)
  
  plot(mu.year.p.summary[2,index],d.summary$p.yr[index],xlim=c(0,1),ylim=c(0,1),
       pch=NA,xlab="Mean (MLE)",ylab="Posterior median",
       main=siteNames[i])
  abline(a=0,b=1)
  points(mu.year.p.summary[2,index],d.summary$p.yr[index],pch=16)
  segments(y0=d.summary$p.yr[index],x0=mu.year.p.summary[1,index],x1=mu.year.p.summary[3,index])
  
}

mtext(paste0("Predicted seedling survival (mean of posterior)"), side = 1, outer = TRUE, line = 1.5)
mtext("Observed seedling survival (mean of observations)", side = 2, outer = TRUE, line = 2.2)


# sample size
par(mfrow=c(1,1),mar=c(4,4,4,4))
plot(d.summary$n.yr,mu.year.summary[2,]-d.summary$p.yr,ylim=c(-1,1),pch=1,
     xlab="Sample size of year*site combination",ylab="Posterior mean - mean of observations")
segments(x0=d.summary$n.yr,y0=mu.year.summary[1,]-d.summary$p.yr,y1=mu.year.summary[3,]-d.summary$p.yr)
abline(h=0)


par(mfrow=c(4,5),mar=c(2,2,1,1))
for(i in 1:20){
  
  index = data$siteYearIndex[data$site==i]
  index = unique(index)
  
  for(j in 1:length(index)){
    index.j = unique(data$siteYearIndex[data$siteYearIndex==index[j]])
    data$fruitplNumber[data$siteYearIndex %in% index.j]
    }
  
  plot(mu.year.summary[2,index],d.summary$p.yr[index],xlim=c(0,1),ylim=c(0,1),
       pch=NA,xlab="Mean (MLE)",ylab="Posterior median",
       main=siteNames[i])
  abline(a=0,b=1)
  points(mu.year.summary[2,index],d.summary$p.yr[index],pch=16)
  segments(y0=d.summary$p.yr[index],x0=mu.year.summary[1,index],x1=mu.year.summary[3,index])
  
}


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Graphical Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ---
# Graphical checks ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "fruitplNumber_sim")
y.obs=data$fruitplNumber
n.obs=data$seedlingNumber
n.iter=dim(y.sim)[1]
subsample = sample(n.iter,5000)
y.sim = y.sim[subsample,]
years=2006:2020

pdf(file=paste0(outputDirectory,"seedlingSurvivalFruitingSimple-ppc-population.pdf"),height=6,width=6)

for(i in 1:length(years)){
  
  par(mfrow = c(1,4),
      oma = c(4,5,0,0) + 0.1,
      mar = c(1,0,1,1) + 0.1)
  
  n.samples = 50
  iter.ind = sample(1:length(subsample),n.samples)
  
  tmpSite = data$site_observed[data$year_observed==i]
  
  for(j in 1:4){
    
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
          lines(x=dens.x$x,y=dens.x$y,lwd=0.25,col='orange')
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
        segments(x0=tb.x,y0=0,y1=tb/sum(tb),col='orange')
        
      }
      
      # ifelse(i%in%c(1),axis(2L),NA)
      
      axis(2,cex=.25,tick=FALSE,line=-1)
      axis(1,cex=.25,tick=FALSE,line=-1)
      
      legend("topright",paste0(siteNames[j],"\n n=",length(p)),bty='n')
    } else {
      plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
    }
  }
  
  mtext(paste0("Probability of seedling survival (",years[i],")"), side = 1, outer = TRUE, line = 1.5)
  mtext("Density", side = 2, outer = TRUE, line = 2.2)
}

dev.off()


# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ---
# Test statistics ----
# ---
chi2.obs=MCMCchains(mcmcSamplesPlots,params=c("chi2.obs"))

n.iter=dim(chi2.obs)[1]
subsample = sample(n.iter,5000)

chi2.obs=MCMCchains(mcmcSamplesPlots,params=c("chi2.obs"))[subsample,]
chi2.sim=MCMCchains(mcmcSamplesPlots,params=c("chi2.sim"))[subsample,]

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
        tmp.sim=y.sim[index]/n.obs[index]
        test.sim=model.fun(tmp.sim,na.rm=TRUE)
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

sims=MCMCchains(mcmcSamplesPlots,params=c("fruitplNumber_sim"))[subsample,]
df=data$fruitplNumber
df2=data$seedlingNumber

seedlings.min=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=min)
seedlings.max=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=max)
seedlings.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)
seedlings.sd=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=sd)

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
seedlings.min.mat<-f.convert(seedlings.min)
seedlings.max.mat<-f.convert(seedlings.max)
seedlings.mean.mat<-f.convert(seedlings.mean)
seedlings.sd.mat<-f.convert(seedlings.sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(15)

pdf(file=paste0(outputDirectory,"seedlingSurvivalFruitingSimplePlots-pvals.pdf"),height=6,width=6)

par(mfrow=c(1,1))
time.sample = 1:length(years)
plot(NA,NA,
     ylim=c(0,1),pch=16,xlim=c(-.5,24.5),
     xlab="",ylab="",
     main="Seedling survival to fruiting",type='n',
     xaxt='n',yaxt='n',frame=FALSE)

start = c(0,5,10,15,20)
offset = seq(0,4,length.out=length(years))
polygon(x=c(4.5,4.5,9.5,9.5),y=c(-1,2,2,-1),col='gray90',border=0)
polygon(x=c(14.5,14.5,19.5,19.5),y=c(-1,2,2,-1),col='gray90',border=0)
box(which="plot",bty="l",col='black')
abline(h=0.5,lty='dotted')

f=function(x,y,width=.2){
  d=boxplot(y,plot=FALSE)
  segments(x0=x,y0=d$stats[1],y1=d$stats[5],lty='solid')
  rect(xleft=x-width/2,xright=x+width/2,
       ybottom=d$stats[2],ytop=d$stats[4],col='white')
  segments(x0=x-width/2,x1=x+width/2,
           y0=d$stats[3],lwd=2)
  
}

for(j in 1:length(years)){
  
  f(x=start[1]+offset[j],y=seedlings.min.mat[,j],width=.2)
  
  f(x=start[2]+offset[j],y=seedlings.max.mat[,j],width=.2)
  
  f(x=start[3]+offset[j],y=seedlings.mean.mat[,j],width=.2)
  
  f(x=start[4]+offset[j],y=seedlings.sd.mat[,j],width=.2)
  
  f(x=start[5]+offset[j],y=p.chi.mat[,j],width=.2)
  
}


axis(1, c(2,7,12,17,22),
     labels = c("min","max","mean","sd","Chi-2"), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     seq(0,1,by=.2), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)

dev.off()

# -------
# Chi-2 ----
# -------

chi2.obs=MCMCchains(mcmcSamplesPlots,params=c("chi2.obs"))[subsample,]
chi2.sim=MCMCchains(mcmcSamplesPlots,params=c("chi2.sim"))[subsample,]
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


pdf(file=paste0(outputDirectory,"seedlingSurvivalFruitingSimple-population.pdf"),height=6,width=8)

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


pdf(file=paste0(outputDirectory,"seedlingSurvivalFruitingSimple-populationyear.pdf"),height=6,width=6)


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

