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

# ---
# - Read in germination and survival estimates ----
# ---

s0 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s0-population-level.RDS")
g1 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/g1-population-level.RDS")
s1 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s1-population-level.RDS")
s2 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s2-population-level.RDS")
s3 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s3-population-level.RDS")
rs <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/reproductiveSuccess-population-year-level-mat.RDS")

# ---
# - Site names ----
# ---
siteAbiotic <- read.csv("data/siteAbioticData.csv",header=TRUE)
siteNames <- siteAbiotic$site

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
g1.hat  <- apply(g1,2,posterior.mode)
s1.hat  <- apply(s1,2,posterior.mode)
s2.hat  <- apply(s2,2,posterior.mode)
s3.hat  <- apply(s3,2,posterior.mode)
s0.hat  <- apply(s0,2,posterior.mode)
rs.hat <- lapply(rs,apply,2,posterior.mode)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

for( k in 1:20){
  pop.index==k
  rs.tmp = rs.hat[[pop.index]]
  yt = rs.tmp
  acf(yt,na.action=na.pass);text(1,.75,siteNames[k])
}

# - +Create matrix to hold objects  ----
fit.point = matrix(NA,ncol=20,nrow=10000)
fit.nosb.point = matrix(NA,ncol=20,nrow=10000)

# - +Create vector to hold objects  ----
lambda.point = lambda.nosb.point = c()
lambda.a.point = lambda.a.nosb.point = c()

# - Quantities with posterior modes ----
# - +With posterior mode  ----
for( k in 1:20){
  # - ++get the population index  ----
  pop.index=k
  
  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = rs.hat[[pop.index]]
  y_t = y_t[!is.na(y_t)]
  
  lambda.a.point[k] = mean(g1.hat[k]*y_t*s0.hat[k]*s1.hat[k]+(1-g1.hat[k])*s2.hat[k]*s3.hat[k])
  lambda.a.nosb.point[k] = mean(y_t*s0.hat[k]*s1.hat[k])
  
  # - ++draw 10000 samples for reproductive success with replacement  ----
  y_t.resample = sample(y_t,1000,replace=TRUE)
  
  # - ++for each resample year, calculate the population growth rate with and without the seed bank  ----
  fit.point[,k] = g1.hat[k]*y_t.resample*s0.hat[k]*s1.hat[k]+(1-g1.hat[k])*s2.hat[k]*s3.hat[k]
  fit.nosb.point[,k] = y_t.resample*s0.hat[k]*s1.hat[k]
  
  # - ++calculate the stochastic population growth rate by approximation, with and without the seed bank  ----
  lambda.point[k] = exp(sum(log(fit.point[,k]))/1000)
  lambda.nosb.point[k] = exp(sum(log(fit.nosb.point[,k]))/1000)
}

# - ++calculate the variance in population growth rates, with and without the seed bank  ----
var.lambda.point =apply(fit.point,2,var)
var.lambda.nosb.point = apply(fit.nosb.point,2,var)


# - Calculate quantities with parameter uncertainty  ----

# - +Draw 1000 samples from posterior  ----
sample.index = sample(1:45000,1000,replace=TRUE)

# - +Create vector to hold objects  ----
lambda.a = lambda.a.nosb = matrix(NA,ncol=20,nrow=length(sample.index))


# - Average population growth rate ----
# - +For each population  ----
for( k in 1:20){
  # - ++get the population index  ----
  pop.index=k
  
  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = rs[[pop.index]][sample.index,]
  y_t = y_t[,!is.na(y_t[1,])]
  
  # - ++calculate the arithmetic mean population growth rate with and without the seed bank  ----
  lambda.a.mat=g1[sample.index,k]*y_t*s0[sample.index,k]*s1[sample.index,k]+(1-g1[sample.index,k])*s2[sample.index,k]*s3[sample.index,k]
  lambda.a[,k] = apply(lambda.a.mat,1,mean)
  
  lambda.a.nosb.mat = y_t*s0[sample.index,k]*s1[sample.index,k]
  lambda.a.nosb[,k] = apply(lambda.a.nosb.mat,1,mean)
  
}


# - +Summary statistics on lambda_a with sb ----
# - ++Calculate the 95% highest posterior density interval ----
HPDI.lambda.a.sb <- apply(lambda.a,2,FUN = function(x) hdi(x, .68))
# - ++Calculate the posterior mode ----
mode.lambda.a.sb <- apply(lambda.a,2, FUN = posterior.mode)

# - +Summary statistics on lambda_a without sb----
# - ++Calculate the 95% highest posterior density interval ----
HPDI.lambda.a.nosb <- apply(lambda.a.nosb,2,FUN = function(x) hdi(x, .68))
# - ++Calculate the posterior mode ----
mode.lambda.a.nosb <- apply(lambda.a.nosb,2, FUN = posterior.mode)

# - +Create matrix to hold objects  ----
lambda.s.sb = array(NA,dim=c(length(sample.index),20))
lambda.s.nosb = array(NA,dim=c(length(sample.index),20))

var.sb = array(NA,dim=c(length(sample.index),20))
var.nosb = array(NA,dim=c(length(sample.index),20))

# - +With parameter uncertainty  ----
for( k in 1:20){
  # - ++get the population index  ----
  pop.index=k
  
  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = rs[[pop.index]][sample.index,]
  y_t = y_t[,!is.na(y_t[1,])]
  g1.tmp = g1[sample.index,k]
  s0.tmp = s0[sample.index,k]
  s1.tmp = s1[sample.index,k]
  s2.tmp = s2[sample.index,k]
  s3.tmp = s3[sample.index,k]
  
  for(j in 1:length(sample.index)){
    
    # - ++draw 10000 samples for reproductive success with replacement  ----
    y_t.resample = sample(y_t[j,],1000,replace=TRUE)
    
    # - ++for each resample year, calculate the population growth rate with and without the seed bank  ----
    fit = g1.tmp[j]*y_t.resample*s0.tmp[j]*s1.tmp[j]+(1-g1.tmp[j])*s2.tmp[j]*s3.tmp[j]
    fit.nosb = y_t.resample*s0.tmp[j]*s1.tmp[j]
    
    # - ++calculate the stochastic population growth rate by approximation, with and without the seed bank  ----
    lambda.s.sb[j,k] = exp(sum(log(fit))/1000)
    lambda.s.nosb[j,k] = exp(sum(log(fit.nosb))/1000)
    
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

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,.5,.5) + 0.1)

for(i in 1:20){
  hist(lambda.s.sb[,i],breaks=100,main='')
  abline(v=posterior.mode(lambda.s.sb[,i]),col='red',lty='dotted')
  abline(v=lambda.point[i],col='red',lwd=2);abline(v=1,col='purple')
  mtext(siteNames[i],side=3,line=-3,adj=.75,cex=.8)
}

for(i in 1:20){
  hist(lambda.s.nosb[,i],breaks=100,main='')
  abline(v=posterior.mode(lambda.s.nosb[,i]),col='red',lty='dotted')
  abline(v=lambda.nosb.point[i],col='red',lwd=2);abline(v=1,col='purple')
mtext(siteNames[i],side=3,line=-3,adj=.75,cex=.8)
}

# - Manuscript figure ----

# - +set font sizes ----
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12
pt5 = 5/12


# tiff(filename=paste0("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/products/figures/betHedgingTest-robust.tif"),
#     height=2.5,width=7,units="in",res=300,compression="lzw",pointsize=12)
tiff(filename=paste0("~/Desktop/figures/betHedgingTest-robust.tif"),
     height=2.5,width=7,units="in",res=300,compression="lzw",pointsize=12)

par(mfrow=c(1,3),mar=c(0,3,0,0)+.1,oma=c(2.5,0,.8,0)+.1,mgp=c(3,.45,0))
plot(mode.lambda.a.sb,mode.lambda.a.nosb,
     pch=21, cex = .75, bg='white', type='n',
     xlim=c(0,7),ylim=c(0,40),
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt8)

segments(x0=HPDI.lambda.a.sb[1,],x1=HPDI.lambda.a.sb[2,],
         y0=mode.lambda.a.nosb, lwd=1)
segments(x0=mode.lambda.a.sb,
         y0=HPDI.lambda.a.nosb[1,], y1=HPDI.lambda.a.nosb[2,],
         lwd=1)
points(mode.lambda.a.sb,mode.lambda.a.nosb,
       pch=21,col='black',bg='white',cex=.75)

axis(1, seq(0,10,by=2), padj = -.5,
     labels = seq(0,10,by=2), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(1, seq(1,9,by=2),labels=FALSE)
axis(2, seq(0,40,by=10),
     labels = seq(0,40,by=10), las = 1, line = 0, hadj= 1.2,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(2, seq(5,35,10),labels=FALSE)

mtext(expression(lambda[a] ~ 'without seedbank'),
      side=2,line=1.5,adj=.5,col='black',cex=pt8)
mtext(expression(lambda[a] ~ 'with seedbank'),
      side=1,line=1.5,col='black',cex=pt8,adj=.5)

abline(a=0,b=1,lty='dotted')
mtext("A.", adj = 0, cex=pt10)

#mtext("A. Arithmetic population growth rate\nis greater without a seedbank",adj=0,cex=.75)

plot(NA,type='n',
     pch=21, cex = .75, bg='white',
     xlim=c(0,80),ylim=c(0,3500),
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt8)

segments(x0=HPDI.var.s.sb[1,],x1=HPDI.var.s.sb[2,],
         y0=mode.var.s.nosb, lwd=1)
segments(x0=mode.var.s.sb,
         y0=HPDI.var.s.nosb[1,], y1=HPDI.var.s.nosb[2,],
         lwd=1)
points(mode.var.s.sb,mode.var.s.nosb,
       pch=21,col='black',bg='white',cex=.75)

axis(1, seq(0,80,by=10), padj = -.5,
     labels = seq(0,80,by=10), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(1, seq(5,75,by=10),labels=FALSE)
axis(2, seq(0,3500,by=500),
     labels = seq(0,3500,by=500), las = 1, line = 0, hadj= 1.2,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(2, seq(250,3250,500),labels=FALSE)

mtext(expression(Var(lambda) ~ 'without seedbank'),
      side=2,line=1.75,adj=.5,col='black',cex=pt8)
mtext(expression(Var(lambda) ~ 'with seedbank'),
      side=1,line=1.5,col='black',cex=pt8,adj=.5)

abline(a=0,b=1,lty='dotted')
mtext("B.", adj = 0, cex=pt10)

#mtext("B. Variation in population growth\nis greater without a seedbank",adj=0,cex=.75)


plot(mode.lambda.s.sb,mode.lambda.s.nosb,
     type='n',
     xlim=c(0,3.5),ylim=c(0,9),
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt8)
abline(a=0,b=1,lty='dotted')


d.plot=data.frame(lambda=mode.lambda.s.sb,
                  lambda.nosb=mode.lambda.s.nosb)

segments(x0=HPDI.lambda.s.sb[1,],x1=HPDI.lambda.s.sb[2,],
         y0=mode.lambda.s.nosb, lwd=1)
segments(x0=mode.lambda.s.sb,
         y0=HPDI.lambda.s.nosb[1,], y1=HPDI.lambda.s.nosb[2,],
         lwd=1)
points(mode.lambda.s.sb,mode.lambda.s.nosb,
       pch=21,col='black',bg='white',cex=1)

axis(1, seq(0,3,by=1), padj = -.5,
     labels = seq(0,3,by=1), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(1, seq(.5,2.5,by=1),labels=FALSE)
axis(2, seq(0,20,by=2),
     labels = seq(0,20,by=2), las = 1, line = 0, hadj= 1.2,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(2, seq(1,19,2),labels=FALSE)

mtext(expression(lambda[s] ~ 'without seedbank'),
      side=2,line=1.5,adj=.5,col='black',cex=pt8)
mtext(expression(lambda[s] ~ 'with seedbank'),
      side=1,line=1.5,col='black',cex=pt8,adj=.5)

mtext("C.", adj = 0, cex=pt10)


#mtext("C. Stochastic population growth is\ngreater without a seedbank (not expected)",adj=0,cex=.75)


dev.off()



