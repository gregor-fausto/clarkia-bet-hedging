#################################################################################
# Script to conduct demographic bet hedging test
#################################################################################

# - Environment ----
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)

# - Source functions for analysis ----
source("scripts/006_testHypotheses/00_utilityFunctions.R")

# - Libraries ----
library(rjags) # jags interface
library(MCMCvis)
library(tidyverse)
library(reshape2)
library(HDInterval)
library(bayesplot)

# - Site names ----
position<-read.csv(file="data/siteAbioticData.csv",header=TRUE) %>%
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteNames <- unique(position$site)

# - Create data frame to match sites and years  ----
siteIndex <- data.frame(site=siteNames,siteIndex=1:20)
yearIndex <- data.frame(year=2006:2020,yearIndex=1:15)
index=expand.grid(1:20,1:15)

# - Read in germination and survival estimates ----
s0 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s0-population-level.RDS")
g1 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/g1-population-level.RDS")
s1 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s1-population-level.RDS")
s2 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s2-population-level.RDS")
s3 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s3-population-level.RDS")
sigma <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/sigma-population-year-level-mat.RDS")
fec <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/combinedF-population-year-level-mat.RDS")
phi <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/phi-population-year-level-mat.RDS")
rs <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/reproductiveSuccess-population-year-level-mat.RDS")

# get mode from full posterior and calculate var in RS based on modes
g1.hat  <- apply(g1,2,posterior.mode)
s1.hat  <- apply(s1,2,posterior.mode)
s2.hat  <- apply(s2,2,posterior.mode)
s3.hat  <- apply(s3,2,posterior.mode)
s0.hat  <- apply(s0,2,posterior.mode)
sigma.mode <- lapply(sigma,apply,2,posterior.mode)
fec.mode <- lapply(fec,apply,2,posterior.mode)
phi.mode <- lapply(phi,apply,2,posterior.mode)
rs.mode <- lapply(rs,apply,2,posterior.mode)

# ---
# - Test of bet hedging with equal resampling; pooled estimates  ----
# ---

# - +Estimate posterior mode for seed survival, germination, and reproductive success  ----

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for( k in 1:20){
 # pop.index=index[,1]==k
  rs.tmp = rs.mode[[k]]
  yt = rs.tmp
  acf(yt,na.action=na.pass);text(1,.75,siteNames[k])
}

# - +Create matrix to hold objects  ----
fit = matrix(NA,ncol=20,nrow=1000)
fit.nosb = matrix(NA,ncol=20,nrow=1000)

# - +Create vector to hold objects  ----
lambda = lambda.nosb = c()
lambda.a = lambda.a.nosb = c()

# - +For each population  ----
for( k in 1:20){
  # - ++get the population index  ----
 # pop.index=index[,1]==k
  rs.hat = rs.mode
  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = rs.hat[[k]]
  y_t = y_t[!is.na(y_t)]

  # - ++calculate the arithmetic mean population growth rate with and without the seed bank  ----
  lambda.a[k] = mean(g1.hat[k]*y_t*s0.hat[k]*s1.hat[k]+(1-g1.hat[k])*s2.hat[k]*s3.hat[k])
  lambda.a.nosb[k] = mean(y_t*s0.hat[k]*s1.hat[k])

  # - ++draw 10000 samples for reproductive success with replacement  ----
  y_t.resample = sample(y_t,1000,replace=TRUE)

  # - ++for each resample year, calculate the population growth rate with and without the seed bank  ----
  fit[,k] = g1.hat[k]*y_t.resample*s0.hat[k]*s1.hat[k]+(1-g1.hat[k])*s2.hat[k]*s3.hat[k]
  fit.nosb[,k] = y_t.resample*s0.hat[k]*s1.hat[k]

  # - ++calculate the stochastic population growth rate by approximation, with and without the seed bank  ----
  lambda[k] = exp(sum(log(fit[,k]))/1000)
  lambda.nosb[k] = exp(sum(log(fit.nosb[,k]))/1000)

}

# - ++calculate the variance in population growth rates, with and without the seed bank  ----
var.lambda =apply(fit,2,var)
var.lambda.nosb = apply(fit.nosb,2,var)

# - Manuscript figure ----

# - +set font sizes ----
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12
pt5 = 5/12

# - ++save manuscript data  ----

demographicBetHedgingTestData = data.frame(site=siteNames,
                                           lambda.a=lambda.a,lambda.a.nosb=lambda.a.nosb,
                                           var.lambda=var.lambda,var.lambda.nosb=var.lambda.nosb,
                                           lambda=lambda,lambda.nosb=lambda.nosb)

outputDirectory <- "outputs/006_hypothesisTesting/"
saveRDS(demographicBetHedgingTestData,paste0(outputDirectory,"demographicBetHedgingTest.RDS"))

# - PLOT FIGURE  ----

tiff(filename=paste0("products/figures/betHedgingTest.tif"),
     height=2.5,width=6.5,units="in",res=800,compression="lzw",pointsize=12)

par(mfrow=c(1,3),mar=c(0,3,0,0)+.1,oma=c(2.5,0,.8,0)+.1,mgp=c(3,.45,0))

# - +PANEL A  ----
plot(lambda.a,lambda.a.nosb,
     pch=21, cex = .75, bg='white',
     xlim=c(0,8),ylim=c(0,40),
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt8,lwd=.75)

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

# - +PANEL B  ----
plot(var.lambda,var.lambda.nosb,
     pch=21, cex = .75, bg='white',
     xlim=c(0,60),ylim=c(0,3100),
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt8,lwd=.75)

axis(1, seq(0,60,by=10), padj = -.5,
     labels = seq(0,60,by=10), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(1, seq(5,55,by=10),labels=FALSE)
axis(2, seq(0,3000,by=500),
     labels = seq(0,3000,by=500), las = 1, line = 0, hadj= 1.2,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(2, seq(250,2750,500),labels=FALSE)

mtext(expression(Var(lambda) ~ 'without seedbank'),
      side=2,line=1.75,adj=.5,col='black',cex=pt8)
mtext(expression(Var(lambda) ~ 'with seedbank'),
      side=1,line=1.5,col='black',cex=pt8,adj=.5)

abline(a=0,b=1,lty='dotted')
mtext("B.", adj = 0, cex=pt10)

# - +PANEL C  ----
plot(lambda,lambda.nosb,
     type='n',
     xlim=c(0,3),ylim=c(0,5),
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt8)
abline(a=0,b=1,lty='dotted')


d.plot=data.frame(lambda=lambda,
                  lambda.nosb=lambda.nosb,
                  site=siteNames)

points(lambda,lambda.nosb,
       pch=21,col='black',bg='white',cex=3,lwd=.75)
points(lambda[c(5,17)],lambda.nosb[c(5,17)],
       pch=21,col='black',bg='white',cex=3,lwd=.75)

axis(1, seq(0,3,by=1), padj = -.5,
     labels = seq(0,3,by=1), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(1, seq(.5,2.5,by=1),labels=FALSE)
axis(2, seq(0,14,by=2),
     labels = seq(0,14,by=2), las = 1, line = 0, hadj= 1.2,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(2, seq(1,13,2),labels=FALSE)

mtext(expression(lambda[s] ~ 'without seedbank'),
      side=2,line=1.5,adj=.5,col='black',cex=pt8)
mtext(expression(lambda[s] ~ 'with seedbank'),
      side=1,line=1.5,col='black',cex=pt8,adj=.5)

mtext("C.", adj = 0, cex=pt10)

dev.off()
