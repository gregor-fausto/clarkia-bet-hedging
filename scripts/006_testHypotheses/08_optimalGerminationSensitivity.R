# -------------------------------------------------------------------
# Sensitivity of density-independent model of germination
# to seed survival in the soil seed bank 
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
s0 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s0-population-level.RDS")
g1 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/g1-population-level.RDS")
s1 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s1-population-level.RDS")
s2 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s2-population-level.RDS")
s3 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s3-population-level.RDS")
rs <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/reproductiveSuccess-population-year-level-mat.RDS")

siteAbiotic <- read.csv("data/siteAbiotic.csv",header=TRUE)

position<-siteAbiotic %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteNames = unique(position$site)

siteIndex <- data.frame(site=siteNames,siteIndex=1:20)
yearIndex <- data.frame(year=2006:2020,yearIndex=1:15)
index=expand.grid(1:20,1:15)

# -------------------------------------------------------------------

fitness <- function(g=g1,s0=s0,s1=s1,s2=s2,s3=s3,rs=rs){
  p1 = g*rs*s0*s1
  p2 = (1-g)*(s2*s3)
  return(as.numeric(p1+p2))
}

posterior.mode = function(x){
  if(!is.na(x[1])){ x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]}else if(is.na(x[1])){
    modeParam <- NA
  }
  return(modeParam)
}

modeEst <- function(x){return(apply(x,2,posterior.mode))}

# use mode as Bayesian estimator
s0.hat  <- apply(s0,2,posterior.mode)
s1.hat  <- apply(s1,2,posterior.mode)
g1.hat  <- apply(g1,2,posterior.mode)
s2.hat  <- apply(s2,2,posterior.mode)
s3.hat  <- apply(s3,2,posterior.mode)
rs.hat <- lapply(rs,apply,2,posterior.mode)


# try optimizer

fr <- function(x,s0){
  x1 = x[1]
  fit <- fitness(g=x1,s0,s1.hat[k],s2.hat[k],s3.hat[k],rs=y_t.resample)
  gm_mean(fit)
}

vec = c()
vec.list = list()
s.seq = seq(0,1,by=0.001)

index=expand.grid(1:20,1:15)

for(k in 1:20){
  
  # - ++get the population index  ----
  pop.index=k
  
  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = rs.hat[[pop.index]]
  y_t = y_t[!is.na(y_t)]
  
  # - ++draw 10000 samples for reproductive success with replacement  ----
  y_t.resample = sample(y_t,1000,replace=TRUE)
  
  for(i in 1:length(s.seq)){
    
    tmp = optimize(fr,c(0,1),s0=s.seq[i],maximum=TRUE)
    vec[i] = tmp$maximum

  }
  vec.list[[k]] = vec
}


geo=order(position$easting)

HPDI.s0 <- apply(s0,2,FUN = function(x) hdi(x, .68))


# add data

# pdf("~/Dropbox/clarkiaSeedBanks/analysis/products/figures/optimalGerminationSensitivity-s0.pdf",
#     height = 7, width = 7)
pdf("~/Desktop/figures/optimalGerminationSensitivity-s0.pdf",
    height = 7, width = 7)

par(mfrow=c(4,5),mar=c(0,.25,.25,0),
    oma=c(4,4,1,1))
for(i in 1:20){
  k = geo[i]
  plot(s.seq,vec.list[[k]],ylim=c(0,1),type='l',xaxt='n',yaxt='n')
  # get closest value in vector
  # https://stackoverflow.com/questions/43472234/fastest-way-to-find-nearest-value-in-vector  
  index = which.min(abs(s.seq - s0.hat[k]))
  
  x.lo = which.min(abs(s.seq - HPDI.s0[1,k]))
  x.hi = which.min(abs(s.seq - HPDI.s0[2,k]))
  
  polygon(c(s.seq[x.lo:x.hi],rev(s.seq[x.lo:x.hi])),
          c(vec.list[[k]][x.lo:x.hi],rep(0,length(x.lo:x.hi))),
          col='gray99',border=FALSE)
  
  segments(x0=s0.hat[k],y0=-1,y1=vec.list[[k]][index],lty='dotted')
  
  
  lines(s.seq[x.lo:x.hi],vec.list[[k]][x.lo:x.hi],lwd=3,col='red')
  
  # add estimates
  points(x=s0.hat[k],y=0,pch=16,cex=1.2)
  segments(x0=HPDI.s0[1,k],x1=HPDI.s0[2,k],y0=0)
  
  
  legend("bottomright",siteNames[k],bty='n',cex=.9)
  
  ifelse(i%in%c(1,6,11,16),axis(2L,las=1),NA)
  ifelse(i%in%c(16:20),axis(1L),NA)
}

mtext("Predicted optimal germination fraction", side = 2, outer = TRUE, line = 2)
mtext(expression('Observed parameter s'[0]), side = 1, outer = TRUE, line = 2.3)

dev.off()

# now compare with s2

fr2 <- function(x,s2){
  x1 = x[1]
  fit <- fitness(g=x1,s0.hat[k],s1.hat[k],s2,s3.hat[k],rs=y_t.resample)
  gm_mean(fit)
}

vec = c()
vec.list = list()
s.seq = seq(0,1,by=0.001)

index=expand.grid(1:20,1:15)

for(k in 1:20){
  
  # - ++get the population index  ----
  pop.index=k
  
  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = rs.hat[[pop.index]]
  y_t = y_t[!is.na(y_t)]
  
  # - ++draw 10000 samples for reproductive success with replacement  ----
  y_t.resample = sample(y_t,1000,replace=TRUE)
  
  for(i in 1:length(s.seq)){
    
    tmp = optimize(fr2,c(0,1),s2=s.seq[i],maximum=TRUE)
    vec[i] = tmp$maximum
    
  }
  vec.list[[k]] = vec
}

  
par(mfrow=c(1,1))
plot(s.seq,vec,type='n')
for(k in 1:20){
  lines(s.seq,vec.list[[k]])
}

HPDI.s2 <- apply(s2,2,FUN = function(x) hdi(x, .68))


# pdf("~/Dropbox/clarkiaSeedBanks/analysis/products/figures/optimalGerminationSensitivity-s2.pdf",
#     height = 7, width = 7)

pdf("~/Desktop/figures/optimalGerminationSensitivity-s2.pdf",
    height = 7, width = 7)

par(mfrow=c(4,5),mar=c(0,.25,.25,0),
    oma=c(4,4,1,1))
for(i in 1:20){
  k = geo[i]
  plot(s.seq,vec.list[[k]],ylim=c(0,1),type='l',xaxt='n',yaxt='n')
  # get closest value in vector
  # https://stackoverflow.com/questions/43472234/fastest-way-to-find-nearest-value-in-vector  
  index = which.min(abs(s.seq - s2.hat[k]))
  
  x.lo = which.min(abs(s.seq - HPDI.s2[1,k]))
  x.hi = which.min(abs(s.seq - HPDI.s2[2,k]))
  
  polygon(c(s.seq[x.lo:x.hi],rev(s.seq[x.lo:x.hi])),
          c(vec.list[[k]][x.lo:x.hi],rep(.5,length(x.lo:x.hi))),
          col='gray99',border=FALSE)
  
  segments(x0=s2.hat[k],y0=-1,y1=vec.list[[k]][index],lty='dotted')
  
  
  lines(s.seq[x.lo:x.hi],vec.list[[k]][x.lo:x.hi],lwd=3,col='red')
  
  # add estimates
  points(x=s2.hat[k],y=0,pch=16,cex=1.2)
  segments(x0=HPDI.s2[1,k],x1=HPDI.s2[2,k],y0=0)
  
  
  legend("bottomleft",siteNames[k],bty='n',cex=.9)
  
  ifelse(i%in%c(1,6,11,16),axis(2L,las=1),NA)
  ifelse(i%in%c(16:20),axis(1L),NA)
}

mtext("Predicted optimal germination fraction", side = 2, outer = TRUE, line = 2)
mtext(expression('Observed parameter s'[2]), side = 1, outer = TRUE, line = 2.3)

dev.off()



# try optimizer

fr3 <- function(x,s1){
  x1 = x[1]
  fit <- fitness(g=x1,s0.hat[k],s1,s2.hat[k],s3.hat[k],rs=y_t.resample)
  gm_mean(fit)
}

vec = c()
vec.list = list()
s.seq = seq(0,1,by=0.001)

index=expand.grid(1:20,1:15)

for(k in 1:20){
  
  # - ++get the population index  ----
  pop.index=k
  
  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = rs.hat[[pop.index]]
  y_t = y_t[!is.na(y_t)]
  
  # - ++draw 10000 samples for reproductive success with replacement  ----
  y_t.resample = sample(y_t,1000,replace=TRUE)
  
  for(i in 1:length(s.seq)){
    
    tmp = optimize(fr3,c(0,1),s1=s.seq[i],maximum=TRUE)
    vec[i] = tmp$maximum
    
  }
  vec.list[[k]] = vec
}


HPDI.s1 <- apply(s1,2,FUN = function(x) hdi(x, .68))

# add data

# pdf("~/Dropbox/clarkiaSeedBanks/analysis/products/figures/optimalGerminationSensitivity-s1.pdf",
#     height = 7, width = 7)
pdf("~/Desktop/figures/optimalGerminationSensitivity-s1.pdf",
    height = 7, width = 7)

par(mfrow=c(4,5),mar=c(0,.25,.25,0),
    oma=c(4,4,1,1))
for(i in 1:20){
  k = geo[i]
  plot(s.seq,vec.list[[k]],ylim=c(0,1),type='l',xaxt='n',yaxt='n')
  # get closest value in vector
  # https://stackoverflow.com/questions/43472234/fastest-way-to-find-nearest-value-in-vector  
  index = which.min(abs(s.seq - s1.hat[k]))
  
  x.lo = which.min(abs(s.seq - HPDI.s1[1,k]))
  x.hi = which.min(abs(s.seq - HPDI.s1[2,k]))
  
  polygon(c(s.seq[x.lo:x.hi],rev(s.seq[x.lo:x.hi])),
          c(vec.list[[k]][x.lo:x.hi],rep(0,length(x.lo:x.hi))),
          col='gray99',border=FALSE)
  
  segments(x0=s1.hat[k],y0=-1,y1=vec.list[[k]][index],lty='dotted')
  
  
  lines(s.seq[x.lo:x.hi],vec.list[[k]][x.lo:x.hi],lwd=3,col='red')
  
  # add estimates
  points(x=s1.hat[k],y=0,pch=16,cex=1.2)
  segments(x0=HPDI.s1[1,k],x1=HPDI.s1[2,k],y0=0)
  
  
  legend("bottomright",siteNames[k],bty='n',cex=.9)
  
  ifelse(i%in%c(1,6,11,16),axis(2L,las=1),NA)
  ifelse(i%in%c(16:20),axis(1L),NA)
}

mtext("Predicted optimal germination fraction", side = 2, outer = TRUE, line = 2)
mtext(expression('Observed parameter s'[1]), side = 1, outer = TRUE, line = 2.3)

dev.off()

# now compare with s3

fr4 <- function(x,s3){
  x1 = x[1]
  fit <- fitness(g=x1,s0.hat[k],s1.hat[k],s2.hat[k],s3,rs=y_t.resample)
  gm_mean(fit)
}

vec = c()
vec.list = list()
s.seq = seq(0,1,by=0.001)

index=expand.grid(1:20,1:15)

for(k in 1:20){
  
  # - ++get the population index  ----
  pop.index=k
  
  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = rs.hat[[pop.index]]
  y_t = y_t[!is.na(y_t)]
  
  # - ++draw 10000 samples for reproductive success with replacement  ----
  y_t.resample = sample(y_t,1000,replace=TRUE)
  
  for(i in 1:length(s.seq)){
    
    tmp = optimize(fr4,c(0,1),s3=s.seq[i],maximum=TRUE)
    vec[i] = tmp$maximum
    
  }
  vec.list[[k]] = vec
}


par(mfrow=c(1,1))
plot(s.seq,vec,type='n')
for(k in 1:20){
  lines(s.seq,vec.list[[k]])
}

HPDI.s3 <- apply(s3,2,FUN = function(x) hdi(x, .68))


# pdf("~/Dropbox/clarkiaSeedBanks/analysis/products/figures/optimalGerminationSensitivity-s3.pdf",
#     height = 7, width = 7)

pdf("~/Desktop/figures/optimalGerminationSensitivity-s3.pdf",
    height = 7, width = 7)

par(mfrow=c(4,5),mar=c(0,.25,.25,0),
    oma=c(4,4,1,1))
for(i in 1:20){
  k = geo[i]
  plot(s.seq,vec.list[[k]],ylim=c(0,1),type='l',xaxt='n',yaxt='n')
  # get closest value in vector
  # https://stackoverflow.com/questions/43472234/fastest-way-to-find-nearest-value-in-vector  
  index = which.min(abs(s.seq - s3.hat[k]))
  
  x.lo = which.min(abs(s.seq - HPDI.s3[1,k]))
  x.hi = which.min(abs(s.seq - HPDI.s3[2,k]))
  
  polygon(c(s.seq[x.lo:x.hi],rev(s.seq[x.lo:x.hi])),
          c(vec.list[[k]][x.lo:x.hi],rep(0,length(x.lo:x.hi))),
          col='gray99',border=FALSE)
  
  segments(x0=s3.hat[k],y0=-1,y1=vec.list[[k]][index],lty='dotted')
  
  
  lines(s.seq[x.lo:x.hi],vec.list[[k]][x.lo:x.hi],lwd=3,col='red')
  
  # add estimates
  points(x=s3.hat[k],y=0,pch=16,cex=1.2)
  segments(x0=HPDI.s3[1,k],x1=HPDI.s3[2,k],y0=0)
  
  
  legend("bottomleft",siteNames[k],bty='n',cex=.9)
  
  ifelse(i%in%c(1,6,11,16),axis(2L,las=1),NA)
  ifelse(i%in%c(16:20),axis(1L),NA)
}

mtext("Predicted optimal germination fraction", side = 2, outer = TRUE, line = 2)
mtext(expression('Observed parameter s'[3]), side = 1, outer = TRUE, line = 2.3)

dev.off()