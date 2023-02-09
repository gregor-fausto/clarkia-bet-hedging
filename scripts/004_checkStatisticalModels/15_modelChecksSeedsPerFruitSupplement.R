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

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Graphical Posterior Predictive Checks
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# ---
# Graphical checks: Seeds per undamaged fruit ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "y_sd.sim")
y.obs=data$sdno
n.iter=dim(y.sim)[1]
years=2006:2020
n.samples = 25
y.sim = y.sim[sample(1:n.iter,n.samples),]

pdf(file=paste0(outputDirectory,"ppc-seedsPerUndamagedFruit.pdf"),height=6,width=6)

for(i in 1:length(years)){
  par(mfrow = c(4,5),
      oma = c(4,5,0,0) + 0.1,
      mar = c(1,0,1,1) + 0.1)
  
  tmpSite = data$site_und_observed[data$year_und_observed==i]
  
  for(j in 1:20){
    
    j.tmp = siteIndex[j]
    if(j.tmp %in% tmpSite){
      index=data$site3==j.tmp&data$year3==i
    
    p.obs=y.obs[index]
    p.sim=(y.sim[,index])

    if(is.matrix(p.sim)){      
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
        lines(x=dens.x$x,y=dens.x$y,lwd=0.25,col='orange')
      }
      
      p = data$sdno[index]
      
      lines(x=dens.obs$x,y=dens.obs$y)
      
      if(length(unique(p))==1){ abline(v=unique(p),lty='dotted',col='black')}
      
     } else if(is.vector(p.sim)){
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
        segments(x0=tb.x,y0=0,y1=tb/sum(tb),col='orange')
        
      }
      
    axis(2,cex=.25,tick=FALSE,line=-1)
    axis(1,cex=.25,tick=FALSE,line=-1)
      
      legend("topright",paste0(siteNames[j],"\n n=",length(p)),bty='n')
    } else {
        plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
      }
    }
  
  mtext(paste0("Counts of seeds per undamaged frut (",years[i],")"), side = 1, outer = TRUE, line = 1.5)
  mtext("Density", side = 2, outer = TRUE, line = 2.2)
} 
dev.off()

# ---
# Graphical checks: Seeds per damaged fruit ----
# ---

y_dam.sim=MCMCchains(mcmcSamples, params = "y_sd_dam.sim")
y_dam.obs=data$sdno_dam
n.iter=dim(y_dam.sim)[1]
years=2013:2020
n.samples = 25
y_dam.sim = y_dam.sim[sample(1:n.iter,n.samples),]

pdf(file=paste0(outputDirectory,"ppc-seedsPerDamagedFruit.pdf"),height=6,width=6)

for(i in 1:length(years)){
  par(mfrow = c(4,5),
      oma = c(4,5,0,0) + 0.1,
      mar = c(1,0,1,1) + 0.1)
  
  tmpSite = data$site_dam_observed[data$year_dam_observed==i]
  
  for(j in 1:20){
    
    j.tmp = siteIndex[j]
    if(j.tmp %in% tmpSite){
      index=data$site4==j.tmp&data$year4==i
      
      p.obs=y_dam.obs[index]
      p.sim=(y_dam.sim[,index])

      if(is.matrix(p.sim)){      
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
          lines(x=dens.x$x,y=dens.x$y,lwd=0.25,col='orange')
        }
        
        p = data$sdno_dam[index]
        
        lines(x=dens.obs$x,y=dens.obs$y)
        
        if(length(unique(p))==1){ abline(v=unique(p),lty='dotted',col='black')}
        
      } else if(is.vector(p.sim)){
        # lines below are for the case where there is only 1 observation
        
        p.sim = p.sim
        all.max.x= max(density(p.sim,na.rm=TRUE)$x)
        all.max.y= max(density(p.sim,na.rm=TRUE)$y)
        
        plot(NA,NA,
             ylim=c(0,all.max.y),xlim=c(0,all.max.x),
             xaxt='n',xlab='',ylab='',yaxt='n')
        
        # plot observed value
        p = data$sdno_dam[index]
        abline(v=p,lty='dotted',col='black')
        
        # overlay simulated values
        tb = table(p.sim)
        tb.x = as.numeric(names(tb))
        segments(x0=tb.x,y0=0,y1=tb/sum(tb),col='orange')
        
      }
      
      axis(2,cex=.25,tick=FALSE,line=-1)
      axis(1,cex=.25,tick=FALSE,line=-1)
      
      legend("topright",paste0(siteNames[j],"\n n=",length(p)),bty='n')
    } else {
      plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE)
    }
  }
  
  mtext(paste0("Counts of seeds per damaged frut (",years[i],")"), side = 1, outer = TRUE, line = 1.5)
  mtext("Density", side = 2, outer = TRUE, line = 2.2)
}
dev.off()


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
sample.index = sample(1:n.iter,1000)
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

sdno.min=f(y.sim=sims,y.obs=df,model.fun=min)
sdno.max=f(y.sim=sims,y.obs=df,model.fun=max)
sdno.mean=f(y.sim=sims,y.obs=df,model.fun=mean)
sdno.sd=f(y.sim=sims,y.obs=df,model.fun=sd)

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
p.chi.mat<-f.convert(p.chi.list)
sdno.min.mat<-f.convert(sdno.min)
sdno.max.mat<-f.convert(sdno.max)
sdno.mean.mat<-f.convert(sdno.mean)
sdno.sd.mat<-f.convert(sdno.sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(length(years))

pdf(file=paste0(outputDirectory,"pvals-seedsPerUndamagedFruit.pdf"),height=6,width=6)

par(mfrow=c(1,1))
time.sample = 1:length(years)
plot(NA,NA,
     ylim=c(0,1),pch=16,xlim=c(-.5,24.5),
     xlab="",ylab="",
     main="Seeds per undamaged fruit",type='n',
     xaxt='n',yaxt='n',frame=FALSE)

start = c(0,5,10,15,20)
offset = seq(0,4,length.out=length(years))#c(-.125,.125)
polygon(x=c(4.5,4.5,9.5,9.5),y=c(-1,2,2,-1),col='gray90',border=0)
polygon(x=c(14.5,14.5,19.5,19.5),y=c(-1,2,2,-1),col='gray90',border=0)
box(which="plot",bty="l",col='black')
abline(h=0.5,lty='dotted')

f.box=function(x,y,width=.2){
  d=boxplot(y,plot=FALSE)
  segments(x0=x,y0=d$stats[1],y1=d$stats[5],lty='solid')
  rect(xleft=x-width/2,xright=x+width/2,
       ybottom=d$stats[2],ytop=d$stats[4],col='white')
  segments(x0=x-width/2,x1=x+width/2,
           y0=d$stats[3],lwd=2)
  
}

for(j in 1:length(years)){
  
  f.box(x=start[1]+offset[j],y=sdno.min.mat[,j],width=.2)
  
  f.box(x=start[2]+offset[j],y=sdno.max.mat[,j],width=.2)
  
  f.box(x=start[3]+offset[j],y=sdno.mean.mat[,j],width=.2)
  
  f.box(x=start[4]+offset[j],y=sdno.sd.mat[,j],width=.2)
  
  f.box(x=start[5]+offset[j],y=p.chi.mat[,j],width=.2)
  
}

axis(1, c(2,7,12,17,22),
     labels = c("min","max","mean","sd","Chi-2"), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     seq(0,1,by=.2), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)

dev.off()


# ---
# Test statistics: Seeds per damaged fruit ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "y_sd_dam.sim")
y.obs=data$sdno_dam
n.iter=dim(y.sim)[1]
years=2013:2020
sample.index = sample(1:n.iter,1000)
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

sdno.min=f(y.sim=sims,y.obs=df,model.fun=min)
sdno.max=f(y.sim=sims,y.obs=df,model.fun=max)
sdno.mean=f(y.sim=sims,y.obs=df,model.fun=mean)
sdno.sd=f(y.sim=sims,y.obs=df,model.fun=sd)

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
p.chi.mat<-f.convert(p.chi.list)
sdno.min.mat<-f.convert(sdno.min)
sdno.max.mat<-f.convert(sdno.max)
sdno.mean.mat<-f.convert(sdno.mean)
sdno.sd.mat<-f.convert(sdno.sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(length(years))

pdf(file=paste0(outputDirectory,"pvals-seedsPerDamagedFruitl.pdf"),height=6,width=6)

par(mfrow=c(1,1))
time.sample = 1:length(years)
plot(NA,NA,
     ylim=c(0,1),pch=16,xlim=c(-.5,24.5),
     xlab="",ylab="",
     main="Seeds per damaged fruit",type='n',
     xaxt='n',yaxt='n',frame=FALSE)

start = c(0,5,10,15,20)
offset = seq(0,4,length.out=length(years))#c(-.125,.125)
polygon(x=c(4.5,4.5,9.5,9.5),y=c(-1,2,2,-1),col='gray90',border=0)
polygon(x=c(14.5,14.5,19.5,19.5),y=c(-1,2,2,-1),col='gray90',border=0)
box(which="plot",bty="l",col='black')
abline(h=0.5,lty='dotted')

f.box=function(x,y,width=.2){
  d=boxplot(y,plot=FALSE)
  segments(x0=x,y0=d$stats[1],y1=d$stats[5],lty='solid')
  rect(xleft=x-width/2,xright=x+width/2,
       ybottom=d$stats[2],ytop=d$stats[4],col='white')
  segments(x0=x-width/2,x1=x+width/2,
           y0=d$stats[3],lwd=2)
  
}

for(j in 1:length(years)){
  
  f.box(x=start[1]+offset[j],y=sdno.min.mat[,j],width=.2)
  
  f.box(x=start[2]+offset[j],y=sdno.max.mat[,j],width=.2)
  
  f.box(x=start[3]+offset[j],y=sdno.mean.mat[,j],width=.2)
  
  f.box(x=start[4]+offset[j],y=sdno.sd.mat[,j],width=.2)
  
  f.box(x=start[5]+offset[j],y=p.chi.mat[,j],width=.2)
  
}

axis(1, c(2,7,12,17,22),
     labels = c("min","max","mean","sd","Chi-2"), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     seq(0,1,by=.2), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)

dev.off()

