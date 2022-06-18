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
data <- readRDS(modelDataDirectory[[grep("seedData.RDS",modelDataDirectory)]])

# - +Read in MCMC samples ----
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedBagExperimentSamplesNonparametric.RDS",mcmcSampleDirectory)]])

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
# *Germination ----
# ---

par(mfrow=c(1,3))

y.sim=MCMCchains(mcmcSamples, params = "seedlingJan_sim")
n.iter=dim(y.sim)[1]

f = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density(x, from= x.min, to = x.max)
  return(dres)
}

pdf(file=paste0(outputDirectory,"germination-population.pdf"),height=6,width=6)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  
  i.tmp = siteIndex[i]
  
  # for each site get index for year 1 germination
  index=data$siteGermination==i.tmp&data$germinationIndex%in%c(1,4,6)
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$totalJan[index], FUN = '/')
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]

  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
  }

  dres = f(data$seedlingJan[index]/data$totalJan[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  legend("topright",paste0(siteNames[i],"\n n=",length(data$totalJan[index])),bty='n')
  
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination (year t+1 after seed production)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)


par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  
  i.tmp = siteIndex[i]
  
  # for each site get index for year 1 germination
  index=data$siteGermination==i.tmp&data$germinationIndex%in%c(2,5)
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$totalJan[index], FUN = '/')
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]
  
  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
  }
  
  dres = f(data$seedlingJan[index]/data$totalJan[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  legend("topright",paste0(siteNames[i],"\n n=",length(data$totalJan[index])),bty='n')
  
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination (year t+2 after seed production)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)



par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i in 1:20){
  
  i.tmp = siteIndex[i]
  
  # for each site get index for year 1 germination
  index=data$siteGermination==i.tmp&data$germinationIndex%in%c(3)
  
  plot(NA,NA,xlim=c(0,1),ylim=c(0,10), main='',
       ylab='',xlab='',xaxt='n',yaxt='n')
  tmp=sweep(y.sim[,index], 2, data$totalJan[index], FUN = '/')
  
  sample.index = sample(1:n.iter,20)
  tmp=tmp[sample.index,]
  
  for(j in 1:20){
    dres=f(tmp[j,])
    lines(dres$x,dres$y,lwd=0.25,col='gray')
  }
  
  dres = f(data$seedlingJan[index]/data$totalJan[index])
  
  lines(dres$x,dres$y,lwd=1,col='black')
  
  legend("topright",paste0(siteNames[i],"\n n=",length(data$totalJan[index])),bty='n')
  
  ifelse(i%in%c(16:20),axis(1L),NA)
  ifelse(i%in%c(1,6,11,16),axis(2L),NA)
}

mtext("Probability of germination (year t+3 after seed production)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

dev.off()

# ---
# *Intact seed counts ----
# ---

# par(mfrow=c(2,3))
# 
 y.sim=MCMCchains(mcmcSamples, params = "y_sim")
# 
# for(i in c(1,7,11)){
#   hist(data$y[data$siteSurvival==1&data$compIndex==i], breaks = 10, 
#        freq = FALSE, main = "Simulated and real data for germination", 
#        xlab = expression(paste("germinant count")), cex.lab = 1.2) 
#   y_sim=y.sim[,data$siteSurvival==1&data$compIndex==i]
#   hist(y_sim, col = "red",add=TRUE,freq=FALSE)
# }
# 
# dev.off()



pdf(file=paste0(outputDirectory,"decay-population.pdf"),height=8,width=6)

par(mfrow=c(5,3),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

group.index = list(1:5,6:10,11:15,16:20)
group = list(siteIndex[1:5],siteIndex[6:10],siteIndex[11:15],siteIndex[16:20])


time.sample = 1:3
months.names = c("Four months", "Twelve months", "Sixteen months")


for(g in 1:4){
  
  group.tmp=group[[g]]
  
  for(j in 1:5){
    
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
           main=ifelse(j==1,months.names[i],''),xaxt='n',xlab='',ylab='',yaxt='n')
      for(h in 1:n.samples){
        m.max=max(tmp[h,],na.rm=TRUE)
        m.min=min(tmp[h,],na.rm=TRUE)
        
        dens.x = density(tmp[h,],from=m.min,to=m.max)
        
        lines(y=dens.x$x,x=dens.x$y,lwd=0.25,col='orange')
      }
      
      p = data$y[index]/data$seedStart[index]
      m.max=max(p,na.rm=TRUE)
      m.min=min(p,na.rm=TRUE)
      
      dens.x = density(p,from=m.min,to=m.max, na.rm=TRUE)
      
      lines(y=dens.x$x,x=dens.x$y)
      
      ifelse(i%in%c(1),axis(2L),NA)
      #ifelse(i%in%c(1,6,11,16),axis(2L),NA)
      ifelse(j%in%c(5),axis(1L),NA)
      if(i==1)  mtext(siteNames[group.index[[g]][j]], side = 3, adj = 0.05, 
                       line = -1.3,cex=.75) else NA
      
    }
    
    
    mtext("Probability of a seed being intact", side = 2, outer = TRUE, line = 2.2)
    mtext("Density", side = 1, outer = TRUE, line = 2.2)
  }
}
dev.off()



# ---
# *s0 ----
# ---

y.sim=MCMCchains(mcmcSamples, params = "plotSeedlings_sim")

pdf(file=paste0(outputDirectory,"s0-population.pdf"),height=6,width=6)

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

mtext("Seeds emerging in permanent plots in January, from seeds produced in year t-1 and t-2", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)
dev.off()


# ---
# Posterior predictive checks ----
# ---

# ---
# *Germination ----
# ---

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.obs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.sim"))
# calculations are rowwise
fit.obs=apply(chi2.obs,1,sum)
fit.sim=apply(chi2.sim,1,sum)
p.chi2.calc=ifelse(fit.sim-fit.obs>=0,1,0)
mean(p.chi2.calc)

gIndex = list(c(1,4,6),c(2,5),c(3))
p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=3)
for(j in 1:20){
  for(i in 1:3){
    tmp.chi2.obs=chi2.obs[,data$siteGermination==j&data$germinationIndex %in% gIndex[[i]]]
    tmp.chi2.sim=chi2.sim[,data$siteGermination==j&data$germinationIndex %in% gIndex[[i]]]
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

# age.sample = 1:3
# plot(age.sample,p.chi.mat[1,],
#      ylim=c(0,1),pch=16,xlim=c(1,3),
#      xlab="Year",ylab="p-Value",
#      main="Germination",type='n')
# for(i in 1:20){
#   points(age.sample+rnorm(1,0,sd=.05),
#          p.chi.mat[i,],pch=19,cex=.5,col='black')
# }
# abline(h=c(.1,.9),lty='dotted')
# abline(h=c(.2,.8),lty='dotted',col='gray')
# 
# 
# 
# y.sim=MCMCchains(mcmcSamples,params=c("seedlingJan_sim"))
# y.obs=data$seedlingJan
# test.sim=apply(y.sim,1,mean)
# test.obs=mean(y.obs)
# p.test = ifelse(test.sim-test.obs>=0,1,0)
# mean(p.test)
# 
# n.iter=dim(y.sim)[1]
# 
# p.test.list = list()
# p.pop = matrix(NA,nrow=n.iter,ncol=3)
# for(j in 1:20){
#   for(i in 1:3){
#     tmp.obs=y.obs[data$siteGermination==j&data$gIndex==i]
#     tmp.sim=y.sim[,data$siteGermination==j&data$gIndex==i]
#     test.obs=mean(tmp.obs)
#     test.sim=apply(tmp.sim,1,mean)
#     p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
#     p.pop[,i] = p.test.calc
#   }
#   p.test.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
# }
# p.test.mat=do.call(rbind,p.test.list)
# 
# 
# par(mfrow=c(1,1))
# age.sample = 1:3
# plot(age.sample,p.test.mat[1,],
#      ylim=c(0,1),pch=16,xlim=c(0,4),
#      xlab="Months",ylab="p-Value",
#      main="Germinant counts",type='n')
# for(i in 1:20){
#   points(age.sample+rnorm(1,0,sd=.05),
#          p.test.mat[i,],pch=1,cex=.5)
# }
# abline(h=c(.1,.9),lty='dotted')



f=function(y.sim=chains,y.obs=data,n.obs=data2,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  p.test.list = list()
  p.pop = matrix(NA,nrow=n.iter,ncol=3)
  for(j in 1:20){
    for(i in 1:3){
      index=data$siteGermination==j&data$germinationIndex %in% gIndex[[i]]
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

germ.min=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=min)
germ.max=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=max)
germ.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)
germ.sd=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=sd)

pdf(file=paste0(outputDirectory,"germination-ppc.pdf"),height=6,width=6)

par(mfrow=c(1,1))

## Make plot of p-values
age.sample = 1:3
plot(age.sample,germ.sd[1,],
     ylim=c(0,1),pch=16,xlim=c(.5,5.5),
     xlab="Months",ylab="p-Value",
     main="Germinant counts",type='n',
     xaxt='n',yaxt='n',frame=FALSE)

polygon(x=c(1.5,1.5,2.5,2.5),y=c(-1,2,2,-1),col='gray90',border=0)
polygon(x=c(3.5,3.5,4.5,4.5),y=c(-1,2,2,-1),col='gray90',border=0)
box(which="plot",bty="l",col='black')
abline(h=0.5,lty='dotted')

col.vec = c("white","gray","black")
offset = c(-.25,0,.25)
  for(j in 1:3){
  points(1+offset[j]+rnorm(20,0,sd=.025),
         germ.min[,j],pch=21,cex=.5,bg=col.vec[j])
    
    points(2+offset[j]+rnorm(20,0,sd=.025),
           germ.max[,j],pch=21,cex=.5,bg=col.vec[j])
    
    points(3+offset[j]+rnorm(20,0,sd=.025),
           germ.mean[,j],pch=21,cex=.5,bg=col.vec[j])
    
    points(4+offset[j]+rnorm(20,0,sd=.025),
           germ.sd[,j],pch=21,cex=.5,bg=col.vec[j])
    
    points(5+offset[j]+rnorm(20,0,sd=.025),
           p.chi.mat[,j],pch=21,cex=.5,bg=col.vec[j])
  }
abline(v=c(1.5,2.5,3.5,4.5))
axis(1, (1:5),
     labels = c("min","max","mean","sd","Chi-2"), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     seq(0,1,by=.2), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
legend(x=.5,y=.1,legend=c("year t+1", "year t+2", "year t+3"),
       pch =c(21,21), cex =.5, pt.bg=col.vec)


# Simulate values of the minimum 
f2=function(y.sim=chains,y.obs=data,n.obs=data2,age=1,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  test.list = list()
  for(j in 1:20){
      index=data$siteGermination==j&data$germinationIndex==gIndex[[age]]
      tmp.obs=y.obs[index]/n.obs[index]
      tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
      test.obs=model.fun(tmp.obs)
      test.sim=apply(tmp.sim,1,model.fun)
      test.list[[j]] = list(test.obs,test.sim)
  }
  return(test.list)
}

sims=MCMCchains(mcmcSamples,params=c("seedlingJan_sim"))
df=data$seedlingJan
df2=data$totalJan

germ.min=f2(y.sim=sims,y.obs=df,n.obs=df2,age=1,model.fun=min)

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

mtext("Minimum probability of germination (year t+1 after seed production)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)

dev.off()

# ---
# *Intact seed counts ----
# ---

chi2.obs=MCMCchains(mcmcSamples,params=c("chi2.yobs"))
chi2.sim=MCMCchains(mcmcSamples,params=c("chi2.ysim"))

n.iter = dim(chi2.obs)[1]

p.chi.list = list()
p.pop = matrix(NA,nrow=n.iter,ncol=6)
for(j in 1:20){
  
  classes = list(c(1,7,11),c(2,8,12),c(3,9),c(4,10),c(5),c(6))
  
  for(i in 1:3){
    cIndex=classes[[i]]
    index = data$siteSurvival==j&data$compIndex%in%cIndex
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
  p.pop = matrix(NA,nrow=n.iter,ncol=6)
  for(j in 1:20){
    
    classes = list(c(1,7,11),c(2,8,12),c(3,9),c(4,10),c(5),c(6))
    
    for(i in 1:3){
      cIndex=classes[[i]]
      index = data$siteSurvival==j&data$compIndex%in%cIndex
      
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

seeds.min=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=min)
seeds.max=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=max)
seeds.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)
seeds.sd=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(3)

pdf(file=paste0(outputDirectory,"decay-ppc.pdf"),height=6,width=6)

par(mfrow=c(1,1))
time.sample = 1:3
plot(NA,NA,
     ylim=c(0,1),pch=16,xlim=c(.5,13),
     xlab="",ylab="p-Value",
     main="Intact seed counts",type='n',
     xaxt='n',yaxt='n',frame=FALSE)

polygon(x=c(3,3,5.5,5.5),y=c(-1,2,2,-1),col='gray90',border=0)
polygon(x=c(8,8,10.5,10.5),y=c(-1,2,2,-1),col='gray90',border=0)
box(which="plot",bty="l",col='black')
abline(h=0.5,lty='dotted')

start = c(1,3.5,6,8.5,11)
offset = seq(0,1.5,length.out=3)
for(j in 1:3){
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
legend(x=11.5,y=1,legend=paste0(c(4,12,16)),
       title="Time (months)",
       pch=21, cex =.5, pt.bg=col.vec)

f2=function(y.sim=chains,y.obs=data,n.obs=data2,class=4,model.fun=mean){
  
  n.iter=dim(y.sim)[1]
  
  test.list = list()
  for(j in 1:20){
    
    index = data$siteSurvival==j&data$compIndex==class
    
    tmp.obs=y.obs[index]/n.obs[index]
    tmp.sim=sweep(y.sim[,index], 2, n.obs[index], FUN = '/')
    test.obs=model.fun(tmp.obs)
    test.sim=apply(tmp.sim,1,model.fun)
    test.list[[j]] = list(test.obs,test.sim)

  }
  return(test.list)
}

sims=MCMCchains(mcmcSamples,params=c("y_sim"))
df=data$y
df2=data$seedStart

seeds.min=f2(y.sim=sims,y.obs=df,n.obs=df2,class=1,model.fun=max)

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

mtext("Maximum intact seed count (January in year t+1 after seed production)", side = 1, outer = TRUE, line = 2.2)
mtext("Density", side = 2, outer = TRUE, line = 2.2)
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

# 
# f=function(y.sim=chains,y.obs=data,n.obs=data2,model.fun=mean){
#   
#   n.iter=dim(y.sim)[1]
#   
#   p.test.list = list()
#   p.pop = matrix(NA,nrow=n.iter,ncol=1)
#   for(j in 1:20){
#     
#       index = data$sitePlot==j
#       tmp.obs=y.obs[index]/n.obs[index]
#       tmp.sim=sweep(as.matrix(y.sim[,index]), 2, n.obs[index], FUN = '/')
#       test.obs=model.fun(tmp.obs)
#       test.sim=apply(tmp.sim,1,model.fun)
#       p.test.calc=ifelse(test.sim-test.obs>=0,1,0)
#       p.pop[,1] = p.test.calc
#     
#     p.test.list[[j]] = apply(p.pop,2,mean,na.rm=TRUE)
#   }
#   p.test.mat=do.call(rbind,p.test.list)
#   return(p.test.mat)
# }
# 
# sims=MCMCchains(mcmcSamples,params=c("plotSeedlings_sim"))
# df=data$plotSeedlings
# df2=data$fec
# 
# seeds.min=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=min)
# seeds.max=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=max)
# seeds.mean=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=mean)
# seeds.sd=f(y.sim=sims,y.obs=df,n.obs=df2,model.fun=sd)

colfunc <- colorRampPalette(c("white", "black"))
col.vec=colfunc(2)

pdf(file=paste0(outputDirectory,"s0-ppc.pdf"),height=6,width=6)

par(mfrow=c(1,1))
time.sample = 1:2
plot(NA,NA,
     ylim=c(0,1),pch=16,xlim=c(0,1),
     xlab="",ylab="",
     main="Emerging seedlings",type='n',
     xaxt='n',yaxt='n',frame=FALSE)

box(which="plot",bty="l",col='black')
abline(h=0.5,lty='dotted')

start = c(.5)
offset = c(0)
for(j in 1){
  
  points(start[1]+offset[j]+rnorm(20,0,sd=.05),
         p.chi.mat[,j],pch=21,cex=.5,bg=col.vec[j])
}
axis(1, (.5),
     labels = c("Chi-2"), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
axis(2, seq(0,1,by=.2),
     seq(0,1,by=.2), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)
legend("topleft",legend=paste0("Seedlings emerging in 2008"),
       title="Year",
       pch =21, cex =.5, pt.bg=col.vec)
dev.off()



