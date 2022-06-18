# -------------------------------------------------------------------
# Density-independent model of germination
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(HDInterval)
library(bayesplot)

# -------------------------------------------------------------------
# Functions for use when analyzing data
# -------------------------------------------------------------------
temporal_variance <- function(x,fun=var){
  apply(x,1,fun)
}

cols_fun <- function(x,fun=var){
  apply(x,2,fun)
}

# geometric var
gsd <- function(x){
  y <- exp(sd(log(x)))
  return(y)
}


f<-function(x="parm",chain){
  chain<-MCMCchains(chain=belowground,params = x)
  p<-boot::inv.logit(chain)
  BCI <- t(apply(p,2,FUN = function(x) quantile(x, c(.025, .5, .975))))
  return(BCI)
}

# geometric mean from 
# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm_mean = function(x, na.rm=TRUE){
exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x), na.rm=na.rm) / length(x))
}

# -------------------------------------------------------------------

# read in samples from posterior distributions
s0 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/s0-pop.RDS")
s1 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/s1-pop.RDS")
g1 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/g1-pop.RDS")
s2 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/s2-pop.RDS")
s3 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/s3-pop.RDS")
rs <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/rsPosterior.RDS")

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) 
siteNames = unique(position$site)

siteIndex <- data.frame(site=siteNames,siteIndex=1:20)
yearIndex <- data.frame(year=2006:2020,yearIndex=1:15)
index=expand.grid(1:20,1:15)

# -------------------------------------------------------------------

fitness <- function(g=g1,s0=s0,s1=s1,s2=s2,s3=s3,rs=rs){
  p1 = g*rs*s1*s0
  p2 = (1-g)*(s2*s3)
  return(as.numeric(p1+p2))
}

posterior.median <- function(x){return(apply(x,2,median))}

posterior.mode = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]
  return(modeParam)
}

modeEst <- function(x){return(apply(x,2,posterior.mode))}

# use mode as Bayesian estimator
s0.hat = modeEst(s0)
g1.hat = modeEst(g1)
s1.hat = modeEst(s1)
s2.hat = modeEst(s2)
s3.hat = modeEst(s3)
rs.hat = modeEst(rs)

# -------------------------------------------------------------------
# calculate autocorrelation of per-capita reproductive success (mode)
# no significant autocorrelation, include in appendix
# -------------------------------------------------------------------

allPopulations = order((position %>% dplyr::select(site,easting))$easting)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for( k in allPopulations){
  pop.index=index[,1]==k
  rs.tmp = rs.hat[pop.index]
  yt = rs.tmp
  acf(yt)
  text(5,.9,(position %>% dplyr::select(site,easting))$site[k])
}

# -------------------------------------------------------------------
# calculate harmonic mean for each population
# -------------------------------------------------------------------

f = function(x){1/mean(1/x)}

pop.harmonicMean =c ()
for( k in 1:20){
  pop.index=index[,1]==k
  rs.tmp = c(rs.hat[pop.index])
  hm = f(rs.tmp)
  pop.harmonicMean[k] = hm
}

par(mfrow=c(1,1))
plot(pop.harmonicMean,s2.hat*s3.hat,type='n');
abline(a=0,b=1);
text(pop.harmonicMean,s2.hat*s3.hat,siteNames)

# -------------------------------------------------------------------
# Simulation to calculate optimal value of germination
# -------------------------------------------------------------------

g = seq(0,1,by=.01)
g.seq = g
fit = c()
g.mat = matrix(NA,nrow=10000,ncol=length(g))
g.sites = list()
yt = c()
n.iter=dim(rs)[1]

# -------------------------------------------------------------------
# for each population k (use mode)
# get the median of per-capita reproductive success in each year
# sample a sequence of 1000 years to use in the analysis for each value of g
# for each g, calculate the growth rate for 1000 years
# -------------------------------------------------------------------

mcmc.draw = sample(1:length(s0[,1]),50)



vec=c()
vec.list = list()
for(k in 1){

  # - ++get the population index  ----
  pop.index=index[,1]==k
  y_t_unc = rs[,pop.index]
  s0_unc = s0[,k]
  s1_unc = s1[,k]
  s2_unc = s2[,k]
  s3_unc = s3[,k]
  
  for(j in 1:50){

  # - ++draw the 15 years of reproductive success estimates  ----

  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = y_t_unc[mcmc.draw[j],]
  s0_est = s0_unc[mcmc.draw[j]]
  s1_est = s1_unc[mcmc.draw[j]]
  s2_est = s2_unc[mcmc.draw[j]]
  s3_est = s3_unc[mcmc.draw[j]]
  
  
  # - ++draw 10000 samples for reproductive success with replacement  ----
  y_t.resample = sample(y_t,10000,replace=TRUE)

  fr <- function(x){
    x1 = x[1]
    fit <- fitness(g=x1,s0=s0_est,s1=s1_est,s2=s2_est,s3=s3_est,rs=y_t.resample)
    gm_mean(fit)
  }
  
  tmp = optimize(fr,c(0,1),maximum=TRUE)
  vec[j] = tmp$maximum
  
  }
vec.list[[k]]  = vec
}

vec.list[[k]]

par(mfrow=c(1,1))
#plot(g.seq,apply(g.mat,2,gm_mean))

apply(g.mat,2,prod)
# -------------------------------------------------------------------
# Plots of geometric mean fitness
# note that without reproductive failure fitness
# remains high at high values of g
# -------------------------------------------------------------------

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(k in 1:20){
  plot(g.seq,apply(g.sites[[k]],2,gm_mean),type='l')
  mtext(siteNames[k],side=3,adjust=-1,line=-2)
}

# -------------------------------------------------------------------
# Plots of arithmetic mean fitness (log)
# -------------------------------------------------------------------

for(k in 1:20){
  plot(g.seq,log(apply(g.sites[[k]],2,mean)),type='l')
}

# -------------------------------------------------------------------
# plots of standard deviation of pop growth rate
# -------------------------------------------------------------------

for(k in 1:20){
  plot(g.seq,apply(g.sites[[k]],2,sd),type='l')
}

# -------------------------------------------------------------------
# Calculate optimal germination fraction
# -------------------------------------------------------------------
gmean <- function(x){apply(x,2,gm_mean)}
site.optima<-lapply(g.sites,gmean)

maxfun <- function(x){g[which(x %in% max(x))]}
optima<-unlist(lapply(site.optima,maxfun))

dev.off()

# -------------------------------------------------------------------
# Make plots
# -------------------------------------------------------------------

# # - ++save manuscript data  ----
# 
# optimalGerminationFractions = data.frame(site=siteNames,optima=optima,g1.hat=g1.hat)
# 
# outputDirectory <- "~/Dropbox/clarkiaSeedBanks/analysis/outputs/007_hypothesisTesting/objects/"
# saveRDS(optimalGerminationFractions,paste0(outputDirectory,"optimalGerminationFractions.RDS"))
# 
# 
# pt12 = 1
# pt10 = 10/12
# pt9 = 9/12
# pt8 = 8/12
# pt7 = 7/12
# pt6 = 6/12
# pt5 = 5/12
# 
# # tiff(filename=paste0("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/products/figures/optimalGerminationFraction.tif"),
# #      height=3.75,width=3.5,units="in",res=300,compression="lzw",pointsize=12)
# 
# par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(2,2.2,.7,1)+.1,mgp=c(3,.45,0))
# 
# plot(NA,NA,xlim=c(.7,1.025),ylim=c(0,.3),
#       #cex.lab = 1.25, cex.axis = 1,    
#      xlab = "",
#      ylab = "",bty='n',
#      xaxt= "n", yaxt="n",)
# 
# points(x=optima,y=g1.hat,pch=21,col='black',bg='white',cex=2.5)
# 
# # Now, define a custom axis
# d.plot=data.frame(optima,
#                   g1.hat,site=siteNames)
# d.plot[c(3:6,10,13,11),1]=1.025
# d.plot[c(3,10,11,13),2]=c(.12,.095,.089,.08)
# d.plot[12,c(1:2)]=d.plot[12,c(1:2)]*c(1,1.05)
# 
# text(d.plot[,1:2],siteNames,cex=4/12)
# 
# box()
# 
# 
# abline(a=0,b=1,lty='dotted')
# 
# axis(1, seq(0,1,by=.1), padj = -.5,
#      labels = seq(0,1,by=.1), line = 0,
#      col = NA, col.ticks = 1, cex.axis = pt8)
# axis(1, seq(.25,1,by=.1),labels=FALSE)
# axis(2, seq(0,1,by=.1),
#      labels = seq(0,1,by=.1), las = 1, line = 0, hadj= 1.2,
#      col = NA, col.ticks = 1, cex.axis = pt8)
# axis(2, seq(.05,1,by=.1),labels=FALSE)
# 
# mtext("Observed germination fraction",
#       side=2,line=1.5,adj=.5,col='black',cex=pt10)
# mtext("Predicted germination fraction",
#       side=1,line=1,adj=.5,col='black',cex=pt10)
# 
# mtext("A.", adj = 0, cex=pt10)
# 
# #dev.off()
# 
# 
# 
# plot(position$easting,optima)
# abline(h=1,col='gray90')
# points(position$easting,mode.lambda.s.sb,pch=16)
# segments(x0=position$easting,y0=HPDI.lambda.s.sb[1,],y1=HPDI.lambda.s.sb[2,])
# 
# 
# # -------------------------------------------------------------------
# # Repeat with low fitness years set to 0
# # -------------------------------------------------------------------
# 
# siteIndex <- data.frame(site=siteNames,siteIndex=1:20)
# yearIndex <- data.frame(year=2006:2020,yearIndex=1:15)
# 
# lowFitnessYears <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/007_hypothesisTesting/objects/lowFitnessYears.RDS")
# 
# lowFitnessYears=lowFitnessYears %>% 
#   dplyr::left_join(yearIndex) %>% 
#   dplyr::left_join(siteIndex)
# 
# # set years without any plants at the population-level to zero fitness
# zeroYearIndex=colnames(rs)%in%paste0("mu[",lowFitnessYears$siteIndex,",",lowFitnessYears$yearIndex,"]")
# rs[,zeroYearIndex] = 0
# 
# index=expand.grid(1:20,1:15)
# 
# # - +Estimate posterior mode for seed survival, germination, and reproductive success  ----
# s0.hat = modeEst(s0)
# g1.hat = modeEst(g1)
# s1.hat = modeEst(s1)
# s2.hat = modeEst(s2)
# s3.hat = modeEst(s3)
# rs.hat = modeEst(rs)
# 
# # -------------------------------------------------------------------
# # Simulation to calculate optimal value of germination
# # -------------------------------------------------------------------
# 
# g = seq(0,1,by=.01)
# g.seq = g
# fit = c()
# g.mat = matrix(NA,nrow=10000,ncol=length(g))
# g.sites = list()
# yt = c()
# n.iter=dim(rs)[1]
# 
# # -------------------------------------------------------------------
# # for each population k (use mode)
# # get the median of per-capita reproductive success in each year
# # sample a sequence of 1000 years to use in the analysis for each value of g
# # for each g, calculate the growth rate for 1000 years
# # -------------------------------------------------------------------
# 
# for( k in 1:20){
#   # - ++get the population index  ----
#   pop.index=index[,1]==k
#   
#   # - ++draw the 15 years of reproductive success estimates  ----
#   y_t = rs.hat[pop.index]
#   
#   # - ++draw 10000 samples for reproductive success with replacement  ----
#   y_t.resample = sample(y_t,10000,replace=TRUE)
#   
#   for( i in 1:length(g)){
#     fit<-fitness(g=g[i],s0.hat[k],s1.hat[k],s2.hat[k],s3.hat[k],rs=y_t.resample)
#     #logfit<-log(fit)
#     g.mat[,i] <- fit
#     #g.mat[,i] <- logfit
#   }
#   g.sites[[k]] <- g.mat
# }
# 
# 
# par(mfrow=c(1,1))
# #plot(g.seq,apply(g.mat,2,gm_mean))
# 
# apply(g.mat,2,prod)
# # -------------------------------------------------------------------
# # Plots of geometric mean fitness
# # note that without reproductive failure fitness
# # remains high at high values of g
# # -------------------------------------------------------------------
# 
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# for(k in 1:20){
#   plot(g.seq,apply(g.sites[[k]],2,gm_mean),type='l')
#   mtext(siteNames[k],side=3,adjust=-1,line=-2)
# }
# 
# # -------------------------------------------------------------------
# # Plots of arithmetic mean fitness (log)
# # -------------------------------------------------------------------
# 
# for(k in 1:20){
#   plot(g.seq,log(apply(g.sites[[k]],2,mean)),type='l')
# }
# 
# # -------------------------------------------------------------------
# # plots of standard deviation of pop growth rate
# # -------------------------------------------------------------------
# 
# for(k in 1:20){
#   plot(g.seq,apply(g.sites[[k]],2,sd),type='l')
# }
# 
# # -------------------------------------------------------------------
# # Calculate optimal germination fraction
# # -------------------------------------------------------------------
# gmean <- function(x){apply(x,2,gm_mean)}
# site.optima<-lapply(g.sites,gmean)
# 
# maxfun <- function(x){g[which(x %in% max(x))]}
# optima<-unlist(lapply(site.optima,maxfun))
# 
# #dev.off()
# 
# # -------------------------------------------------------------------
# # Make plots
# # -------------------------------------------------------------------
# 
# # - ++save manuscript data  ----
# 
# optimalGerminationFractionsLowFitness = data.frame(site=siteNames,optima=optima,g1.hat=g1.hat)
# 
# outputDirectory <- "~/Dropbox/clarkiaSeedBanks/analysis/outputs/007_hypothesisTesting/objects/"
# saveRDS(optimalGerminationFractionsLowFitness,paste0(outputDirectory,"optimalGerminationFractions-lowFitness.RDS"))
# 
# # tiff(filename=paste0("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/products/figures/optimalGerminationFraction-lowFitness.tif"),
# #      height=3.75,width=3.5,units="in",res=300,compression="lzw",pointsize=12)
# 
# par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(2,2.2,.7,1)+.1,mgp=c(3,.45,0))
# 
# plot(NA,NA,xlim=c(.6,1.025),ylim=c(0,.3),
#      #cex.lab = 1.25, cex.axis = 1,    
#      xlab = "",
#      ylab = "",bty='n',
#      xaxt= "n", yaxt="n",)
# 
# points(x=optima,y=g1.hat,pch=21,col='black',bg='white',cex=2.5)
# 
# # Now, define a custom axis
# d.plot=data.frame(optima,
#                   g1.hat,site=siteNames)
# d.plot[c(3:6,10,13),1]=1.03
# d.plot[c(3,10),2]=c(.12,.096)
# #d.plot[12,c(1:2)]=d.plot[12,c(1:2)]*c(1,1.05)
# 
# text(d.plot[,1:2],siteNames,cex=4/12)
# 
# box()
# 
# 
# abline(a=0,b=1,lty='dotted')
# 
# axis(1, seq(0,1,by=.1), padj = -.5,
#      labels = seq(0,1,by=.1), line = 0,
#      col = NA, col.ticks = 1, cex.axis = pt8)
# axis(1, seq(.25,1,by=.1),labels=FALSE)
# axis(2, seq(0,1,by=.1),
#      labels = seq(0,1,by=.1), las = 1, line = 0, hadj= 1.2,
#      col = NA, col.ticks = 1, cex.axis = pt8)
# axis(2, seq(.05,1,by=.1),labels=FALSE)
# 
# mtext("Observed germination fraction",
#       side=2,line=1.5,adj=.5,col='black',cex=pt10)
# mtext("Predicted germination fraction",
#       side=1,line=1,adj=.5,col='black',cex=pt10)
# 
# mtext("B.", adj = 0, cex=pt10)
# 
# #dev.off()
# 
# plot(position$easting,optima)
