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
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("viabilityTrialSamples.RDS",mcmcSampleDirectory)]])

# - +Read in site names ----
siteAbiotic <- read.csv("data/siteAbiotic.csv",header=TRUE)
siteNames <- siteAbiotic$site

# create variable that arranges populations by easting
siteIndex <- order(siteAbiotic$easting,decreasing=FALSE)
siteNames = unique(siteAbiotic$site)[siteIndex]

# ---
# Graphical checks ----
# ---

# ---
# *Germination trials ----
# ---

par(mfrow=c(1,3))

y.sim=MCMCchains(mcmcSamples, params = "germCount_sim")
n.iter=dim(y.sim)[1]

f = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density(x, from= x.min, to = x.max)
  return(dres)
}

pdf(file=paste0(outputDirectory,"viability-trials-germination.pdf"),height=6,width=6)

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
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]

  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
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
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]
  
  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
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
  
  index=data$siteViab==tmp.i&data$indexViab==3
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$germStart[index], FUN = '/')
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]
  
  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
  }
  
  dres = f(data$germCount[index]/data$germStart[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  mtext(siteNames[i], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination in lab trial\n (3 years after seed production, round 1)", 
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
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]
  
  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
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
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]
  
  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
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
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]
  
  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
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
# *Viability trials ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "viabStain_sim")
n.iter=dim(y.sim)[1]

pdf(file=paste0(outputDirectory,"viability-trials-viability.pdf"),height=6,width=6)

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
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]
  
  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
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
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]
  
  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
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
  
  index=data$siteViab_v==tmp.i&data$indexViab_v==3
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$viabStart[index], FUN = '/')
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]
  
  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
  }
  
  dres = f(data$viabStain[index]/data$viabStart[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  mtext(siteNames[i], side = 3, adj = 0.05, 
        line = -1.3,cex=.75)
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of viability, conditional on not germinating, in lab trial\n (3 years after seed production, round 1)", 
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
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]
  
  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
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
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]
  
  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
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
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]
  
  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
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


# ---
# Posterior predictive checks ----
# ---

# ---
# *Germination ----
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

germ.min=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=min)
germ.max=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=max)
germ.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)
germ.sd=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=sd)

pdf(file=paste0(outputDirectory,"viabilityTrial-germination-ppc.pdf"),height=6,width=6)

par(mfrow=c(1,1))

## Make plot of p-values
trial.sample = 1:6
plot(trial.sample,germ.min[1,],
     ylim=c(0,1),pch=16,xlim=c(.5,13),
     xlab="",ylab="p-Value",
     main="Germinant counts, lab trials",type='n',
     xaxt='n',yaxt='n',frame=FALSE)

box(which="plot",bty="l",col='black')
abline(h=0.5,lty='dotted')

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(6)

start = c(1,3.5,6,8.5,11)
offset = seq(0,1.5,length.out=6)
  for(j in 1:6){
  points(start[1]+offset[j]+rnorm(20,0,sd=.025),
         germ.min[,j],pch=21,cex=.5,bg=col.vec[j])

    points(start[2]+offset[j]+rnorm(20,0,sd=.025),
           germ.max[,j],pch=21,cex=.5,bg=col.vec[j])

    points(start[3]+offset[j]+rnorm(20,0,sd=.025),
           germ.mean[,j],pch=21,cex=.5,bg=col.vec[j])

    points(start[4]+offset[j]+rnorm(20,0,sd=.025),
           germ.sd[,j],pch=21,cex=.5,bg=col.vec[j])

    points(start[5]+offset[j]+rnorm(20,0,sd=.025),
           p.chi.mat[,j],pch=21,cex=.5,bg=col.vec[j])
  }
abline(v=c(3,5.5,8,10.5))
axis(1, (start+.75),
     labels = c("min","max","mean","sd","Chi-2"), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     seq(0,1,by=.2), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)


legend("topright",legend=c("1/1","1/2","1/3","2/1","2/2","3/1"),
       title="Seed age/experiment round",
       pch=21, cex =.5, pt.bg=col.vec)


# Simulate values of the minimum 
f2=function(y.sim=chains,y.obs=data,n.obs=data2,class=1,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  test.list = list()
  for(j in 1:20){
      index=data$siteViab==j&data$indexViab==class
      tmp.obs=y.obs[index]/n.obs[index]
      tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
      test.obs=model.fun(tmp.obs)
      test.sim=apply(tmp.sim,1,model.fun)
      test.list[[j]] = list(test.obs,test.sim)
  }
  return(test.list)
}

sims=MCMCchains(mcmcSamples,params=c("germCount_sim"))

germ.min=f2(y.sim=sims,y.obs=df,n.obs=df2,class=6,model.fun=min)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(1,0,1,1) + 0.1)

for(i in 1:20){
  tmp = germ.min[[i]]
  hist( tmp[[2]], breaks=100,main='',freq=FALSE,ylab='',yaxt='n')
  abline(v=tmp[[1]],col='red',lwd=2,lty='dotted')
  y.max = max(hist( tmp[[2]], breaks=100,plot=FALSE)$density)*.75
  x.max = max(hist( tmp[[2]], breaks=100,plot=FALSE)$mids)*.5
  
  mtext(siteNames[i], side = 3, adj = 0.1, 
        line = -1.3,cex=.75)
}

mtext("Minimum probability of germination (age 3 seeds, round 1)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

dev.off()

# ---
# *Viable, staining seeds count ----
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
  
 # classes = list(c(1,2,3),c(4, 5),c(6))
  
  for(i in 1:6){
  #  vIndex=classes[[i]]
    index = data$siteViab_v==j&data$indexViab_v==i#%in%vIndex
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
    
   # classes = list(c(1,2,3),c(4, 5),c(6))
    
    for(i in 1:6){
    #  vIndex=classes[[i]]
      index = data$siteViab_v==j&data$indexViab_v==i#%in%vIndex
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

seeds.min=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=min)
seeds.max=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=max)
seeds.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)
seeds.sd=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(6)

pdf(file=paste0(outputDirectory,"viabilityTrial-viability-ppc.pdf"),height=6,width=6)

par(mfrow=c(1,1))
time.sample = 1:6
plot(NA,NA,
     ylim=c(0,1),pch=16,xlim=c(.5,13),
     xlab="",ylab="p-Value",
     main="Counts of stained seeds (test for viability)",type='n',
     xaxt='n',yaxt='n',frame=FALSE)

polygon(x=c(3,3,5.5,5.5),y=c(-1,2,2,-1),col='gray90',border=0)
polygon(x=c(8,8,10.5,10.5),y=c(-1,2,2,-1),col='gray90',border=0)
box(which="plot",bty="l",col='black')
abline(h=0.5,lty='dotted')

start = c(1,3.5,6,8.5,11)
offset = seq(0,1.5,length.out=6)
for(j in 1:6){
  points(start[1]+offset[j]+rnorm(20,0,sd=.025),
         seeds.min[,j],pch=21,cex=.5,bg=col.vec[j])

  points(start[2]+offset[j]+rnorm(20,0,sd=.025),
         seeds.max[,j],pch=21,cex=.5,bg=col.vec[j])
  
  points(start[3]+offset[j]+rnorm(20,0,sd=.025),
         seeds.mean[,j],pch=21,cex=.5,bg=col.vec[j])
  
  points(start[4]+offset[j]+rnorm(20,0,sd=.025),
         seeds.sd[,j],pch=21,cex=.5,bg=col.vec[j])
  
  points(start[5]+offset[j]+rnorm(20,0,sd=.025),
         p.chi.mat[,j],pch=21,cex=.5,bg=col.vec[j])
}
abline(v=c(3,5.5,8,10.5))
axis(1, (start+.75),
     labels = c("min","max","mean","sd","Chi-2"), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     seq(0,1,by=.2), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
legend("topright",legend=c("1/1","1/2","1/3","2/1","2/2","3/1"),
       title="Seed age/experiment round",
       pch=21, cex =.5, pt.bg=col.vec)

f2=function(y.sim=chains,y.obs=data,n.obs=data2,class=4,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  test.list = list()
  for(j in 1:20){
    
    index = data$siteViab_v==j&data$indexViab_v==class
    
    tmp.obs=y.obs[index]/n.obs[index]
    tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
    test.obs=model.fun(tmp.obs)
    test.sim=apply(tmp.sim,1,model.fun)
    test.list[[j]] = list(test.obs,test.sim)

  }
  return(test.list)
}

seeds.min=f2(y.sim=sims,y.obs=df,n.obs=df2,class=6,model.fun=min)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(1,0,1,1) + 0.1)

for(i in 1:20){
  tmp = seeds.min[[i]]
  hist( tmp[[2]], breaks=100,main='',freq=FALSE,ylab='',yaxt='n')
  abline(v=tmp[[1]],col='red',lwd=2,lty='dotted')
  y.max = max(hist( tmp[[2]], breaks=100,plot=FALSE)$density)*.75
  x.max = min(hist( tmp[[2]], breaks=100,plot=FALSE)$mids)*1
  
  text(x=x.max,y=y.max,
       siteNames[i],pos=4)
}

mtext("Minimum probability of seed staining in viability trial (age 3 seeds, round 1)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)
dev.off()

