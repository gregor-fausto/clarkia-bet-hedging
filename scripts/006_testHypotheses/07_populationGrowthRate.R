#################################################################################
################################################################################
################################################################################
# Simulation
#
# Scripts by Gregor Siegmund
# fausto.siegmund@gmail.com
# last updated 06-01-2021
#################################################################################
#################################################################################
#################################################################################

library(tidyverse)
library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)
# -------------------------------------------------------------------
# Functions for use when analyzing data
# -------------------------------------------------------------------
posterior.mode = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]
  return(modeParam)
}

modeEst <- function(x){return(apply(x,2,posterior.mode))}

# ---
# - Read in germination and survival estimates ----
# ---

s0 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/s0-pop.RDS")
s1 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/s1-pop.RDS")
g1 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/g1-pop.RDS")
s2 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/s2-pop.RDS")
s3 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/s3-pop.RDS")
rs <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/rsPosterior.RDS")

# ---
# - Site names ----
# ---
position<-read.csv(file="data/siteAbioticData.csv",header=TRUE) %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteNames <- unique(position$site)

# ---
# - Create data frame to match sites and years  ----
# ---
siteIndex <- data.frame(site=siteNames,siteIndex=1:20)
yearIndex <- data.frame(year=2006:2020,yearIndex=1:15)
index=expand.grid(1:20,1:15)

# ---
# - Test of bet hedging with equal resampling; pooled estimates  ----
# ---

# - +Estimate posterior mode for seed survival, germination, and reproductive success  ----
s0.hat = modeEst(s0)
g1.hat = modeEst(g1)
s1.hat = modeEst(s1)
s2.hat = modeEst(s2)
s3.hat = modeEst(s3)
rs.hat = modeEst(rs)

# - +Create matrix to hold objects  ----
fit.sb.point = matrix(NA,ncol=20,nrow=10000)
fit.nosb.point = matrix(NA,ncol=20,nrow=10000)

# - +Create vector to hold objects  ----
lambda.sb.point = lambda.nosb.point = c()

# - Quantities with posterior modes ----
# - +With posterior mode  ----
for( k in 1:20){
  # - ++get the population index  ----
  pop.index=index[,1]==k
  
  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = rs.hat[pop.index]
  
  # - ++draw 10000 samples for reproductive success with replacement  ----
  y_t.resample = sample(y_t,10000,replace=TRUE)
  
  # - ++for each resample year, calculate the population growth rate with and without the seed bank  ----
  fit.sb.point[,k] = g1.hat[k]*y_t.resample*s0.hat[k]*s1.hat[k]+(1-g1.hat[k])*s2.hat[k]*s3.hat[k]
  fit.nosb.point[,k] = y_t.resample*s0.hat[k]*s1.hat[k]
  
  # - ++calculate the stochastic population growth rate by approximation, with and without the seed bank  ----
  lambda.sb.point[k] = exp(sum(log(fit.sb.point[,k]))/10000)
  lambda.nosb.point[k] = exp(sum(log(fit.nosb.point[,k]))/10000)
}

# - ++calculate the variance in population growth rates, with and without the seed bank  ----
var.lambda.point =apply(fit.sb.point,2,var)
var.lambda.nosb.point = apply(fit.nosb.point,2,var)


# - Calculate quantities with parameter uncertainty  ----

# - +Draw 1000 samples from posterior  ----
sample.index = 1:4500
sample.length = length(sample.index)

# - +Create matrix to hold objects  ----
lambda.s.sb = array(NA,dim=c(sample.length,20))
lambda.s.nosb = array(NA,dim=c(sample.length,20))

var.sb = array(NA,dim=c(sample.length,20))
var.nosb = array(NA,dim=c(sample.length,20))

# - +With parameter uncertainty  ----
for( k in 1:20){
  # - ++get the population index  ----
  pop.index=index[,1]==k
  
  # - ++draw the 15 years of reproductive success estimates  ----
  y_t.tmp = rs[sample.index,pop.index]
  g1.tmp = g1[sample.index,k]
  s0.tmp = s0[sample.index,k]
  s1.tmp = s1[sample.index,k]
  s2.tmp = s2[sample.index,k]
  s3.tmp = s3[sample.index,k]
  
  for(j in 1:length(sample.index)){
    
    # - ++draw 10000 samples for reproductive success with replacement  ----
    y_t.resample = sample(y_t.tmp[j,],10000,replace=TRUE)
    
    # - ++for each resample year, calculate the population growth rate with and without the seed bank  ----
    # fit[,j,k] = g1.tmp[j]*y_t.resample*s0.tmp[j]*s1.tmp[j]+(1-g1.tmp[j])*s2.tmp[j]*s3.tmp[j]
    #fit.nosb[,j,k] = y_t.resample*s0.tmp[j]*s1.tmp[j]
    fit = g1.tmp[j]*y_t.resample*s0.tmp[j]*s1.tmp[j]+(1-g1.tmp[j])*s2.tmp[j]*s3.tmp[j]
    fit.nosb = y_t.resample*s0.tmp[j]*s1.tmp[j]
    
    # - ++calculate the stochastic population growth rate by approximation, with and without the seed bank  ----
    lambda.s.sb[j,k] = exp(sum(log(fit))/10000)
    lambda.s.nosb[j,k] = exp(sum(log(fit.nosb))/10000)
    
    var.sb[j,k] = var(fit)
    var.nosb[j,k] =  var(fit.nosb)
    
  }
}

# - +Summary statistics on lambda_s with sb ----
# - ++Calculate the 95% highest posterior density interval ----
HPDI.lambda.s.sb <- apply(lambda.s.sb,2,FUN = function(x) hdi(x, .68))
# - ++Calculate the posterior mode ----
mode.lambda.s.sb <- apply(lambda.s.sb,2, FUN = posterior.mode)

# - +Summary statistics on lambda_s without sb----
# - ++Calculate the 95% highest posterior density interval ----
HPDI.lambda.s.nosb <- apply(lambda.s.nosb,2,FUN = function(x) hdi(x, .68))
# - ++Calculate the posterior mode ----
mode.lambda.s.nosb <- apply(lambda.s.nosb,2, FUN = posterior.mode)

# - +Summary statistics on var(lambda) with sb ----
# - ++Calculate the 95% highest posterior density interval ----
HPDI.var.s.sb <- apply(var.sb,2,FUN = function(x) hdi(x, .68))
# - ++Calculate the posterior mode ----
mode.var.s.sb <- apply(var.sb,2, FUN = posterior.mode)

# - +Summary statistics on var(lambda) without sb----
# - ++Calculate the 95% highest posterior density interval ----
HPDI.var.s.nosb <- apply(var.nosb,2,FUN = function(x) hdi(x, .68))
# - ++Calculate the posterior mode ----
mode.var.s.nosb <- apply(var.nosb,2, FUN = posterior.mode)


dev.off()
par(mfrow=c(1,1))

plot(position$easting,lambda.sb.point,ylim=c(0,3))
abline(h=1,col='gray90')
points(position$easting,mode.lambda.s.sb,pch=16)
segments(x0=position$easting,y0=HPDI.lambda.s.sb[1,],y1=HPDI.lambda.s.sb[2,])

plot(position$easting,lambda.nosb.point,ylim=c(0,18))
abline(h=1,col='gray90')
points(position$easting,mode.lambda.s.nosb,pch=16)
segments(x0=position$easting,y0=HPDI.lambda.s.nosb[1,],y1=HPDI.lambda.s.nosb[2,])

plot(position$easting,var.lambda.point,ylim=c(0,80))
abline(h=1,col='gray90')
points(position$easting,mode.var.s.sb,pch=16)
segments(x0=position$easting,y0=HPDI.var.s.sb[1,],y1=HPDI.var.s.sb[2,])

plot(position$easting,var.lambda.nosb.point)
abline(h=1,col='gray90')
points(position$easting,mode.var.s.nosb,pch=16)
segments(x0=position$easting,y0=HPDI.var.s.nosb[1,],y1=HPDI.var.s.nosb[2,])
