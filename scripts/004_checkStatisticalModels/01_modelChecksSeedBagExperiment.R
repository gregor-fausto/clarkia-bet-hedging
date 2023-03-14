####
####
# Script for model checks for seed bag experiment
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
modelDataDirectory = "outputs/002_fitStatisticalModels/data/"
outputDirectory = "outputs/004_checkStatisticalModels/"

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

# - Function to plot Bayesian p-values ----
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
           y=diagnostic.test.mat[,j],
           pch=21,cex=.8,
           bg=ifelse(diagnostic.test.mat[,j]>0.95|
                       diagnostic.test.mat[,j]<0.05,'orange','gray95'))
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

# - +Function to obtain posterior predictive densities----
f.dens.plot = function(x){
  x.max <- max(x)
  x.min <- min(x)
  dres <- density(x, from= x.min, to = x.max)
  return(dres)
}

# - Germination ----
# *Posterior predictive checks ----

# - ++Get chi-squared calculated in JAGS ----
chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sim"))
n.iter=dim(chi2.obs)[1]

# - ++Calculate Bayesian p-value per population/germination index ----
# chi-squared
p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=6)
for(j in 1:20){
  for(i in 1:6){
    tmp.chi2.obs=chi2.obs[,data$siteGermination==j&data$germinationIndex %in% i]
    tmp.chi2.sim=chi2.sim[,data$siteGermination==j&data$germinationIndex %in% i]
    fit.obs=apply(tmp.chi2.obs,1,sum)
    fit.sim=apply(tmp.chi2.sim,1,sum)
    p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
    p.pop[,i] = p.chi2.calc
  }
  p.chi.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
}
p.chi.mat=do.call(rbind,p.chi.list)

# - +++function to calculate Bayesian p-value for specific test statistic ----
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

# - ++get simulated observations and observed data ----
sims=MCMCchains(mcmcSamples,params=c("seedlingJan_sim"))
df=data$seedlingJan
df2=data$totalJan

# - ++Calculate Bayesian p-value per population/germination index ----
# mean
germ.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)

# - +++function to convert list to matrix ----
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

# - +Plot Bayesian p-vals for germination ----
pdf(file=paste0(outputDirectory,"pvals-germinationSeedBagExperiment.pdf"),height=3,width=6)

par(mfrow = c(1,2),
    oma = c(2,2.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

# germination g1 in 2006, 2007, and 2008 
# correspond to germinationIndex 1, 4, 6 
# see lines 102-107 in scripts/001_prepareDataForModels/01_prepDataForSeedModel.R
f.plot(diagnostic.test.mat = germ.mean.mat[,c(1,4,6)],years=2006:2008)
mtext("Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat[,c(1,4,6)],years=2006:2008)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)

dev.off()

# *Posterior predictive distribution ----
# - ++get simulated observations----
y.sim=MCMCchains(mcmcSamples, params = "seedlingJan_sim")
n.iter=dim(y.sim)[1]

# - +Plot posterior predictive distribution for germination ----
# all three years in the seed bag burial experiment

pdf(file=paste0(outputDirectory,"ppc-germinationSeedBagExperiment.pdf"),height=6,width=6)

# - ++ Plot for 2006 ----
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

# for each population
for(i in 1:20){
  
  # order populations geographically
  i.tmp = siteIndex[i]
  
  # for each site get index for germination of age 1 seeds in 2006
  index=data$siteGermination==i.tmp&data$germinationIndex%in%c(1)
  
  # calculate the probability of germination for simulated data
  tmp=sweep(y.sim[,index], 2, data$totalJan[index], FUN = '/')
  # sample 25 datasets
  sample.index = sample(1:n.iter,25)
  tmp=tmp[sample.index,]
  # calculate the density
  list.dens=apply(tmp,1,density,na.rm=TRUE)
  # get max values for plot dimensions
  all.max.y=max(unlist(lapply(list.dens, "[", "y")))
  all.max.x=max(unlist(lapply(list.dens, "[", "x")))
  
  plot(NA,NA,
       ylim=c(0,all.max.y),xlim=c(0,1),
       xaxt='n',xlab='',ylab='',yaxt='n')
  
  # for each dataset, plot the posterior predictive distribution
  # coloring it according to whether or not the p-value is extreme
  # use the p-values for the mean test statistic
  for(j in 1:25){
    dres=f.dens.plot(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,
          col=ifelse(germ.mean.mat[i,1]>.95|germ.mean.mat[i,1]<.05,
                     'orange','gray75'))
  }
  
  # get the distribution for the observed data
  # plot it as a thick black line overlaying the ppd
  dres = f.dens.plot(data$seedlingJan[index]/data$totalJan[index])
  lines(dres$x,dres$y,lwd=1,col='black')
  
  # add a legend with the site name
  # and number of bags from which the data comes
  legend("topright",paste0(siteNames[i],"\n n=",length(data$totalJan[index])),bty='n')
  
  # add axes on the outer plots
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

# axis labels
mtext("Probability of germination in 2006 (year t+1 after seed production)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

# - ++ Plot for 2007 ----
# code follows the same structure as the plot for 2006
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

for(i in 1:20){
  
  i.tmp = siteIndex[i]
  
  # for each site get index for germination of age 1 seeds in 2007
  index=data$siteGermination==i.tmp&data$germinationIndex%in%c(4)
  
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
    dres=f.dens.plot(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,
          col=ifelse(germ.mean.mat[i,4]>.95|germ.mean.mat[i,4]<.05,
                     'orange','gray75'))
  }
  
  dres = f.dens.plot(data$seedlingJan[index]/data$totalJan[index])
  lines(dres$x,dres$y,lwd=1,col='black')
  
  legend("topright",paste0(siteNames[i],"\n n=",length(data$totalJan[index])),bty='n')
  
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination in 2007 (year t+1 after seed production)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

# - ++ Plot for 2008 ----
# code follows the same structure as the plot for 2006
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

for(i in 1:20){
  
  i.tmp = siteIndex[i]
  
  # for each site get index for germination of age 1 seeds in 2008
  index=data$siteGermination==i.tmp&data$germinationIndex%in%c(6)
  
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
    dres=f.dens.plot(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,
          col=ifelse(germ.mean.mat[i,6]>.95|germ.mean.mat[i,6]<.05,
                     'orange','gray75'))
  }
  
  dres = f.dens.plot(data$seedlingJan[index]/data$totalJan[index])
  lines(dres$x,dres$y,lwd=1,col='black')
  
  legend("topright",paste0(siteNames[i],"\n n=",length(data$totalJan[index])),bty='n')
  
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination in 2008 (year t+1 after seed production)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

dev.off()

# - +Plot observation-level posterior predictive distribution for one populations ----
# this shows the simulated observations for each seed bag for OKRW
# across all three years of the seed bag burial experiment

pdf(file=paste0(outputDirectory,"ppc-germinationSeedBagExperiment-populationLevel.pdf"),height=4,width=6)

par(mfrow = c(1,3),
    oma = c(2,2,0,0) + 0.1,
    mar = c(1,0,1,1) + 0.1)

# for OKRW 
for(i in 9){
  
  # for all three years
  for(k in c(1:3)){
    
    # get the right site
    i.tmp = siteIndex[i]
    
    # get index for age 1 germination in 2006, 2007 and 2008
    index=data$siteGermination==i.tmp&data$germinationIndex%in%c(1,4,6)[k]
    
    # get the observed and simulated values
    p.obs=data$seedlingJan[index]
    p.sim=as.matrix(y.sim[,index])
    
    # the plot is organized so that observations (indexed by 1, 2, ...) are on the x-axis
    # and the counts of seedlings are on the y-axis
    plot(NA,NA,
         xlim=c(1,length(p.obs)+1),
         ylim=c(0,max(p.sim,p.obs)),
         xaxt='n',xlab='',ylab='',yaxt='n')
    
    # create a gray background on alternating observations
    rect(xleft=seq(1,50,by=2),
         xright=seq(1,50,by=2)+1,
         ybottom=-1000,ytop=1000,col='gray99',border='gray99')
    
    # for each observation
    # draw 25 simulated values
    # and plot their frequency as segments
    for(h in 1:length(p.obs)){
      index.rand = sample(1:45000,25)
      tmp<-table(p.sim[index.rand,h])
      segments(y0=as.numeric(names(tmp)),y1=as.numeric(names(tmp)),
               x0=h, x1=h+(tmp/max(tmp))*.8)
    }
    
    # add the observed data as points
    points(1:length(p.obs),p.obs,pch=16)
    
    # make the plot look nice and add labels
    box()
    axis(2,cex=.25,tick=FALSE,line=-1)
    axis(1,cex=.25,tick=FALSE,line=-1)
    
    # add the population name, year, the number of bags from which data came
    # and the total number of seeds across those bags
    legend("topright",paste0((2006:2008)[k],"\n",siteNames[i],
                             "\n bags=",length(p.obs),
                             "\n total seeds=",sum(p.obs)),bty='n')
  } 
}

# axis label
mtext(paste0("Germinant counts"), side = 2, outer = TRUE, line = 1)
mtext("Observation (index)", side = 1, outer = TRUE, line = 0.2)

dev.off()

# - Intact seed counts ----
# *Posterior predictive checks ----
# structure of the code here largely follows that for checks for the germination model
# I only note differences

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.yobs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.ysim"))
n.iter = dim(chi2.obs)[1]

p.chi.list = list()
# there are 12 time points with observations for the intact seed dataset
p.pop = matrix(NA,nrow=n.iter,ncol=12)
for(j in 1:20){
  for(i in 1:12){
    # observations are indexed by compIndex from 1-12
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

# - +++function to calculate Bayesian p-value for specific test statistic ----
# this is different because of the indexing
f=function(y.sim=chains,y.obs=data,n.obs=data2,model.fun=mean){
  n.iter=dim(y.sim)[1]
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=12)
  for(j in 1:20){
    for(i in 1:12){
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

# - ++get simulated observations and observed data ----
sims=MCMCchains(mcmcSamples,params=c("y_sim"))
df=data$y
df2=data$seedStart

# - ++Calculate Bayesian p-value per population/survival index ----
# mean
seeds.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)

# - +++function to convert list to matrix ----
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

# - +Plot Bayesian p-vals for intact seeds/seed survival ----
# three plots - one each for seed survival in January and October of year t+1
# and one for seed survival in January of year t+2
pdf(file=paste0(outputDirectory,"pvals-seedSurvivalSeedBagExperiment.pdf"),height=3,width=6)

par(mfrow = c(1,2),
    oma = c(2,3.5,0,0) + 0.1,
    mar = c(0,.5,1,1) + 0.1,
    mgp=c(3,.5,0))

# seed survival s1 in 2006, 2007, and 2008 
# correspond to compIndex  c(1,7,11)
f.plot(diagnostic.test.mat = seeds.mean.mat[,c(1,7,11)],years=2006:2008)
mtext("Intact seeds in January (year t+1) \n Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("A. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat[,c(1,7,11)],years=2006:2008)
box(which="plot",bty="l",col='black')
title("B. Chi-squared",adj=0)

# seed survival s2 in 2006, 2007, and 2008 
# correspond to compIndex  c(2,8,12)
f.plot(diagnostic.test.mat = seeds.mean.mat[,c(2,8,12)],years=2006:2008)
mtext("Intact seeds in October (year t+1) \n Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("C. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat[,c(2,8,12)],years=2006:2008)
box(which="plot",bty="l",col='black')
title("D. Chi-squared",adj=0)

# seed survival s3 in 2007, 2008
# correspond to compIndex  c(3,9)
f.plot(diagnostic.test.mat = seeds.mean.mat[,c(3,9)],years=2007:2008)
mtext("Intact seeds in January (year t+2) \n Bayesian p-value", side=2, line=2, cex.lab=1,las=0, col="black")
title("E. Mean", adj=0)

f.plot(diagnostic.test.mat = p.chi.mat[,c(3,9)],years=2007:2008)
box(which="plot",bty="l",col='black')
title("F. Chi-squared",adj=0)

dev.off()

# *Posterior predictive distribution ----
# - ++get simulated observations----
y.sim=MCMCchains(mcmcSamples, params = "y_sim")
n.iter=dim(y.sim)[1]

# - +Plot posterior predictive distribution for intact seeds ----
# these plots are organized as follows
# for an individual population, I plot 3 panels
# 1st panel is for intact seeds in January t+1
# 2nd panel is for intact seeds in October t+1
# 3rd panel is for intact seeds in January t+2
# this shows the probability of seeds surviving after 4, 12 and 16 months
# and highlights how issues with model with are related across the experiment
# because the model jointly fits all observations

pdf(file=paste0(outputDirectory,"ppc-seedSurvivalSeedBagExperiment.pdf"),height=6,width=6)

par(mfrow=c(2,3),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

# names of time
months.names = c("Four months", "Twelve months", "Sixteen months")

# get populations that have extreme test statistics here
x1=(1:20)[seeds.mean.mat[,2]>.95&!is.na(seeds.mean.mat[,2])]
# and then get a population with a normal test statistic
add.index <- c(sample((1:20)[!(1:20%in%x1)],1))
# create a vector of 1 population with extreme p-value and one that is not extreme
group.tmp=c(sample(x1,1),add.index)

# for both of these populations
for(j in 1:2){
  
  j.tmp=group.tmp[j]
  
  # sample 25 simulated datasets
  n.samples = 25
  iter.ind = sample(1:n.iter,n.samples)

  # for each of the 3 timepoints  
  for(i in 1:3){
    
    # get the index corresponding to the site and times
    index = data$siteSurvival==j.tmp&data$compIndex==i
    
    # get the simulated values
    y_sim=y.sim[,index]
    # calculate the probability of seeds being intact from simulations
    tmp=sweep(y_sim, 2, data$seedStart[index], FUN = '/')
    tmp=tmp[iter.ind,]
    # get density info for plots
    list.dens=apply(tmp,1,density)
    all.max=max(unlist(lapply(list.dens, "[", "y")))
    
    # plots are vertical so that the probability is on the y-axis and 
    # density is on the x-axis
    plot(NA,NA,
         xlim=c(0,all.max),ylim=c(0,1),
         main=ifelse(j==1,months.names[i],''),
         xaxt='n',xlab='',ylab='',yaxt='n')
    
    # for each dataset, plot the posterior predictive distribution
    # coloring it according to whether or not the p-value is extreme
    # use the p-values for the mean test statistic
    # use both orange and purple to note direction of test statistic
    for(h in 1:n.samples){
      m.max=max(tmp[h,],na.rm=TRUE)
      m.min=min(tmp[h,],na.rm=TRUE)
      dens.x = density(tmp[h,],from=m.min,to=m.max)
      lines(y=dens.x$x,x=dens.x$y,lwd=0.25,
            col=if(seeds.mean.mat[j.tmp,i]>.95){'orange'}
            else if(seeds.mean.mat[j.tmp,i]<.05){'purple'}
            else{ 'gray75'})
    }
  
    # get the distribution for the observed data
    # plot it as a black line overlaying the ppd  
    p = data$y[index]/data$seedStart[index]
    m.max=max(p,na.rm=TRUE)
    m.min=min(p,na.rm=TRUE)
    dens.x = density(p,from=m.min,to=m.max, na.rm=TRUE)
    lines(y=dens.x$x,x=dens.x$y)
    
    # add labels and make plot look nice
    ifelse(i%in%c(1),axis(2L),NA)
    ifelse(j%in%c(3),axis(1L),NA)
    if(i==1){
      legend("topright",
             paste0(siteAbiotic$site[group.tmp[j]],"\n n=",
                    length(data$seedStart[index])),bty='n') 
    } else {
      legend("topright",
             paste0("n=",
                    length(data$seedStart[index])),bty='n') 
    }
    axis(1,cex=.25,tick=FALSE,line=-1)
  }
}
mtext("Probability of a seed being intact (first round seed cohort)", side = 2, outer = TRUE, line = 2.2)
mtext("Density", side = 1, outer = TRUE, line = 2.2)

dev.off()

# - Seedlings emerging in plots/s0 ----
# *Posterior predictive checks ----
# structure here largely follows that of the earlier checks

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.plot.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.plot.sim"))
n.iter = dim(chi2.obs)[1]

# - ++Calculate Bayesian p-value per population ----
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

# - +++function to calculate Bayesian p-value for specific test statistic ----
f=function(y.sim=chains,y.obs=data,model.fun=mean){
  n.iter=dim(y.sim)[1]
  p.test.list = list()
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

# - ++Get data ----
sims=MCMCchains(mcmcSamples,params=c("y_sim"))
df=data$plotSeedlings

# - ++Calculate Bayesian p-value per population/survival index ----
seedlings.mean=f(y.sim=sims,y.obs=df,model.fun=mean)

# convert to matrix 
seedlings.mean.mat<-seedlings.mean
p.chi.mat<-do.call(rbind,p.chi.list)

# - +Plot posterior predictive distribution for seedlings emerging in plots ----
pdf(file=paste0(outputDirectory,"pvals-seedlingsEmergingPlots.pdf"),height=3,width=6)

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

# *Posterior predictive distribution ----

y.sim=MCMCchains(mcmcSamples, params = "plotSeedlings_sim")

pdf(file=paste0(outputDirectory,"ppc-seedlingsEmergingPlots.pdf"),height=6,width=6)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

for(i in 1:20){
  
  # get indexes
  n.samples = 25
  iter.ind = sample(1:n.iter,n.samples)
  
  # get site index
  i.tmp = siteIndex[i]
  index=data$sitePlot==i.tmp
  
  # draw simulated data and get density
  y_sim=y.sim[,index]
  tmp=y_sim
  tmp=tmp[iter.ind,]
  list.dens=apply(tmp,1,density)
  
  # get density of observed
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
  
  # plot each simulated dataset
  for(h in 1:n.samples){
    m.max=max(tmp[h,])
    m.min=min(tmp[h,])
    dens.x = density(tmp[h,],from=m.min,to=m.max)
    lines(x=dens.x$x,y=dens.x$y,lwd=0.25,
          col=if(seedlings.mean.mat[i.tmp,1]>.95){'orange'}
          else if(seedlings.mean.mat[i.tmp,1]<.05){'purple'}
          else{ 'gray75'})
  }
  
  # add observed data
  p = data$plotSeedlings[index]
  m.max=max(p)
  m.min=min(p)
  dens.x = density(p,from=m.min,to=m.max)
  lines(x=dens.x$x,y=dens.x$y)
  
  # make plot look nice
  legend("topright",paste0(siteNames[i],"\n n=",sum(data$plotSeedlings[index])),bty='n')
  axis(1,cex=.25,tick=FALSE,line=-1)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Seeds emerging in permanent plots in January 2008, \n from seeds produced in year t-1 and t-2", side = 1, outer = TRUE, line = 3.5)
mtext("Density", side = 2, outer = TRUE, line = 2.2)
dev.off()

