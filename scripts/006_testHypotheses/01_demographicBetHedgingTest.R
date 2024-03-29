#################################################################################
# Script to conduct demographic bet hedging test
# Produces Figure 3A-C
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
s0 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s0-population-level.RDS")
g1 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/g1-population-level.RDS")
s1 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s1-population-level.RDS")
s2 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s2-population-level.RDS")
s3 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s3-population-level.RDS")
sigma <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/sigma-population-year-level-mat.RDS")
fec <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/combinedF-population-year-level-mat.RDS")
phi <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/phi-population-year-level-mat.RDS")
perCapitaRS <- readRDS("outputs/005_calculatePopulationModelParameters/04_reproductiveSuccess/reproductiveSuccessWithCorrectionForMissingness-populationYear-mat.RDS")

# - +Estimate posterior mode for seed survival, germination, and reproductive success  ----
g1.hat  <- apply(g1,2,posterior.mode)
s1.hat  <- apply(s1,2,posterior.mode)
s2.hat  <- apply(s2,2,posterior.mode)
s3.hat  <- apply(s3,2,posterior.mode)
s0.hat  <- apply(s0,2,posterior.mode)
sigma.mode <- lapply(sigma,apply,2,posterior.mode)
fec.mode <- lapply(fec,apply,2,posterior.mode)
phi.mode <- lapply(phi,apply,2,posterior.mode)
perCapitaRS.mode <- lapply(perCapitaRS,apply,2,posterior.mode)

# ---
# - Test of bet hedging with equal resampling of years  ----
# ---

# check for autocorrelation
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for( k in 1:20){
 # pop.index=index[,1]==k
  rs.tmp = perCapitaRS.mode[[k]]
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
  rs.hat = perCapitaRS.mode
  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = rs.hat[[k]]
  y_t = y_t[!is.na(y_t)]

  # - ++calculate the arithmetic mean population growth rate with and without the seed bank  ----
  lambda.a[k] = mean(g1.hat[k]*y_t*s0.hat[k]*s1.hat[k]+(1-g1.hat[k])*s2.hat[k]*s3.hat[k])
  lambda.a.nosb[k] = mean(y_t*s0.hat[k]*s1.hat[k])

  # - ++draw 1000 samples for reproductive success with replacement  ----
  y_t.resample = sample(y_t,1000,replace=TRUE)

  # - ++for each resample year, calculate the population growth rate with and without the seed bank  ----
  fit[,k] = g1.hat[k]*y_t.resample*s0.hat[k]*s1.hat[k]+(1-g1.hat[k])*s2.hat[k]*s3.hat[k]
  fit.nosb[,k] = y_t.resample*s0.hat[k]*s1.hat[k]

  # - ++calculate the stochastic population growth rate by approximation, with and without the seed bank  ----
  lambda[k] = exp(sum(log(fit[,k]))/1000)
  lambda.nosb[k] = exp(sum(log(fit.nosb[,k]))/1000)

}

var.lambda =apply(fit,2,var)
var.lambda.nosb = apply(fit.nosb,2,var)

# - +For each population  ----
lambda.mat = lambda.nosb.mat = matrix(NA,nrow = 1000,ncol=20)
var.lambda.mat = var.lambda.nosb.mat = matrix(NA,nrow = 1000,ncol=20)
for(j in 1:1000){
  for( k in 1:20){
    # - ++get the population index  ----
    # pop.index=index[,1]==k
    rs.hat = perCapitaRS.mode
    # - ++draw the 15 years of reproductive success estimates  ----
    y_t = rs.hat[[k]]
    y_t = y_t[!is.na(y_t)]
    
    # - ++draw 1000 samples for reproductive success with replacement  ----
    y_t.resample = sample(y_t,1000,replace=TRUE)
    
    # - ++for each resample year, calculate the population growth rate with and without the seed bank  ----
    fit[,k] = g1.hat[k]*y_t.resample*s0.hat[k]*s1.hat[k]+(1-g1.hat[k])*s2.hat[k]*s3.hat[k]
    fit.nosb[,k] = y_t.resample*s0.hat[k]*s1.hat[k]
    
    # - ++calculate the stochastic population growth rate by approximation, with and without the seed bank  ----
    lambda[k] = exp(sum(log(fit[,k]))/1000)
    lambda.nosb[k] = exp(sum(log(fit.nosb[,k]))/1000)
    
  }
  
  lambda.mat[j,] = lambda
  lambda.nosb.mat[j,] = lambda.nosb
  var.lambda.mat[j,] = apply(fit,2,var)
  var.lambda.nosb.mat[j,] = apply(fit.nosb,2,var)
}

# - ++calculate the mean of population growth rates after 1000 simulations  ----
lambda.mean = apply(lambda.mat,2,mean)
lambda.nosb.mean = apply(lambda.nosb.mat,2,mean)

# - ++calculate mean of the variance in population growth rates, with and without the seed bank  ----
var.lambda.mean = apply(var.lambda.mat,2,mean)
var.lambda.nosb.mean = apply(var.lambda.nosb.mat,2,mean)

# - Manuscript figure  ----

# - +set font sizes ----
pt14 = 14/12
pt12 = 1
pt11 = 11/12
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12
pt5 = 5/12

# - ++save manuscript data  ----

demographicBetHedgingTestData = data.frame(site=siteNames,
                                           lambda.a=lambda.a,lambda.a.nosb=lambda.a.nosb,
                                           var.lambda=var.lambda.mean,var.lambda.nosb=var.lambda.nosb.mean,
                                           lambda=lambda.mean,lambda.nosb=lambda.nosb.mean)

outputDirectory <- "outputs/006_hypothesisTesting/"
saveRDS(demographicBetHedgingTestData,paste0(outputDirectory,"demographicBetHedgingTest.RDS"))

# uncomment to replot the figure 
demographicBetHedgingTestData<-readRDS(paste0(outputDirectory,"demographicBetHedgingTest.RDS"))
siteNames=demographicBetHedgingTestData$site;
lambda.a=demographicBetHedgingTestData$lambda.a;lambda.a.nosb=demographicBetHedgingTestData$lambda.a.nosb;
var.lambda.mean=demographicBetHedgingTestData$var.lambda;var.lambda.nosb.mean=demographicBetHedgingTestData$var.lambda.nosb;
lambda.mean=demographicBetHedgingTestData$lambda;lambda.nosb.mean=demographicBetHedgingTestData$lambda.nosb

# - Plot Figure 3  ----

tiff(filename=paste0("products/figures/betHedgingTest.tif"),
     height=6.5,width=2.5,units="in",res=800,compression="lzw",pointsize=12,
     family = "ArialMT")

par(mfrow=c(3,1),mar=c(3,3.75,0,0)+.1,oma=c(0,0,1,0)+.1,mgp=c(3,.45,0))

# - +PANEL A  ----
plot(lambda.a,lambda.a.nosb,
     pch=21, cex = .75, bg='white',
     xlim=c(0,8),ylim=c(0,44),
     xlab = "",type='n',
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt8,lwd=.75)

points(lambda.a,lambda.a.nosb,
       pch=21,cex=1,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 1),lwd=0)
points(lambda.a,lambda.a.nosb,lwd = 0.5,
       pch=21,col='black',cex=1,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 0.5))

axis(1, seq(0,10,by=2), padj = 0, tcl = -0.4,
     labels = seq(0,10,by=2), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt12)
axis(1, seq(1,9,by=2),labels=FALSE,tck=-0.02)
axis(2, seq(0,40,by=10), tcl = -0.4,
     labels = seq(0,40,by=10), las = 1, line = 0, hadj= 1.1,
     col = NA, col.ticks = 1, cex.axis = pt12)
axis(2, seq(5,35,10),labels=FALSE,tck=-0.02)

mtext(expression(lambda[a] ~ 'without seedbank'),
      side=2,line=2,adj=.5,col='black',cex=pt11)
mtext(expression(lambda[a] ~ 'with seedbank'),
      side=1,line=1.9,col='black',cex=pt11,adj=.5)

abline(a=0,b=1,lty='dotted')
mtext("A.", adj = 0, cex=pt12)

# - +PANEL B  ----
plot(var.lambda.mean,var.lambda.nosb.mean,
     pch=21, cex = .75, bg='white',
 #    xlim=c(0,60),ylim=c(0,3100),
     xlab = "", type='n',
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt8,lwd=.75)

points(var.lambda.mean,var.lambda.nosb.mean,
       pch=21,cex=1,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 1),lwd=0)
points(var.lambda.mean,var.lambda.nosb.mean,
       pch=21,col='black',cex=1,lwd = 0.5,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 0.5))

axis(1, seq(0,60,by=10), padj = 0, tcl = -0.4,
     labels = seq(0,60,by=10), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt12)
axis(1, seq(5,55,by=10),labels=FALSE,tck=-0.02)
axis(2, seq(0,3000,by=500), tcl = -0.4,
     labels = seq(0,3000,by=500), las = 1, line = 0, hadj= 1.1,
     col = NA, col.ticks = 1, cex.axis = pt12)
axis(2, seq(250,2750,500),labels=FALSE,tck=-0.02)

mtext(expression(Var(lambda) ~ 'without seedbank'),
      side=2,line=2.5,adj=.5,col='black',cex=pt11)
mtext(expression(Var(lambda) ~ 'with seedbank'),
      side=1,line=1.9,col='black',cex=pt11,adj=.5)

abline(a=0,b=1,lty='dotted')
mtext("B.", adj = 0, cex=pt12)

# - +PANEL C  ----
plot(lambda.mean,lambda.nosb.mean,
     type='n',
     xlim=c(0,4),ylim=c(0,20),
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt12)
abline(a=0,b=1,lty='dotted')


d.plot=data.frame(lambda=lambda.mean,
                  lambda.nosb=lambda.nosb.mean,
                  site=siteNames)

points(lambda.mean,lambda.nosb.mean,
       pch=21,cex=1,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 1),lwd=0)
points(lambda.mean,lambda.nosb.mean,
       pch=21,col='black',cex=1, lwd = 0.5,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 0.5))

axis(1, seq(0,4,by=1), padj = 0, tcl = -0.4,
     labels = seq(0,4,by=1), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt12)
axis(1, seq(.5,3.5,by=1),labels=FALSE,tck=-0.02)
axis(2, seq(0,20,by=4),tcl = -0.4,
     labels = seq(0,20,by=4), las = 1, line = 0, hadj= 1.1,
     col = NA, col.ticks = 1, cex.axis = pt12)
axis(2, seq(2,20,4),labels=FALSE,tck = -.02)

mtext(expression(lambda[s] ~ 'without seedbank'),
      side=2,line=2,adj=.5,col='black',cex=pt11)
mtext(expression(lambda[s] ~ 'with seedbank'),
      side=1,line=1.9,col='black',cex=pt11,adj=.5)

mtext("C.", adj = 0, cex=pt12)

dev.off()

# - Revision/response to reviewers  ----
# a reviewer asked us about the range of points observed in panel B
# to respond to this question, we produced a figure that zoomed into 
# a section of panel B and better reveals the range of the data
# we include this figure in our response to reviewers

tiff(filename=paste0("products/figures/rev-betHedgingTestB.tif"),
     height=2.5,width=4.5,units="in",res=800,compression="lzw",pointsize=12,
     family = "ArialMT")

par(mfrow=c(1,2),mar=c(0,3,0,0)+.1,oma=c(2.5,0,.8,0)+.1,mgp=c(3,.45,0))
# - +original PANEL B plus gray highlight ----
plot(var.lambda.mean,var.lambda.nosb.mean,
     pch=21, cex = .75, bg='white',
     #xlim=c(0,10),ylim=c(0,250),
     xlab = "", type='n',
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt8,lwd=.75)

rect(-10,-1000,10,600,col='gray90',border=0)

points(var.lambda.mean,var.lambda.nosb.mean,
       pch=21,cex=1,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 1),lwd=0)
points(var.lambda.mean,var.lambda.nosb.mean,
       pch=21,col='black',cex=1,lwd = 0.5,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 0.5))

axis(1, seq(0,60,by=10), padj = -.5,
     labels = seq(0,60,by=10), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(1, seq(5,55,by=10),labels=FALSE,tck = -.02)
axis(2, seq(0,3000,by=500),
     labels = seq(0,3000,by=500), las = 1, line = 0, hadj= 1.2,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(2, seq(250,2750,500),labels=FALSE,tck = -.02)

mtext(expression(Var(lambda) ~ 'without seedbank'),
      side=2,line=1.75,adj=.5,col='black',cex=pt8)
mtext(expression(Var(lambda) ~ 'with seedbank'),
      side=1,line=1.5,col='black',cex=pt8,adj=.5)

box()

abline(a=0,b=1,lty='dotted')
mtext("A. Figure 3B: full plot", adj = 0, cex=pt9)


# - +zoomed PANEL B plus gray highlight ----
plot(var.lambda.mean,var.lambda.nosb.mean,
     pch=21, cex = .75, bg='white',
     xlim=c(0,10),ylim=c(0,600),
     xlab = "", type='n',
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt8,lwd=.75)

rect(-10,-1000,20,700,col='gray90',border=0)

points(var.lambda.mean,var.lambda.nosb.mean,
       pch=21,cex=1,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 1),lwd=0)
points(var.lambda.mean,var.lambda.nosb.mean,
       pch=21,col='black',cex=1,lwd = 0.5,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 0.5))

axis(1, seq(0,10,by=2), padj = -.5,
     labels = seq(0,10,by=2), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(1, seq(1,9,by=2),labels=FALSE,tck = -.02)
axis(2, seq(0,600,by=200),
     labels = seq(0,600,by=200), las = 1, line = 0, hadj= 1.2,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(2, seq(100,500,100),labels=FALSE,tck = -.02)

mtext(expression(Var(lambda) ~ 'without seedbank'),
      side=2,line=1.75,adj=.5,col='black',cex=pt8)
mtext(expression(Var(lambda) ~ 'with seedbank'),
      side=1,line=1.5,col='black',cex=pt8,adj=.5)

box()

abline(a=0,b=1,lty='dotted')
mtext("B. Figure 3B: zoomed in plot", adj = 0, cex=pt9)

dev.off()
