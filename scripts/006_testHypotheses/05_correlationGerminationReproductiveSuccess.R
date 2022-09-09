# -------------------------------------------------------------------
# Analysis of correlation between germination and RS
# Produces Figure 5B
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)

# -------------------------------------------------------------------
# Load required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)
library(tidyverse)
library(reshape2)
library(HDInterval)
library(bayesplot)
library(rethinking)

# - +source scripts ----
source("scripts/006_testHypotheses/00_utilityFunctions.R")

# - Site names by position ----
siteAbiotic <- read.csv("data/siteAbioticData.csv",header=TRUE)

position<-siteAbiotic %>%
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteIndex <- order(position$easting,decreasing=FALSE)
siteNames = unique(position$site)[siteIndex]

# -------------------------------------------------------------------
# Read in samples from posterior distributions
# -------------------------------------------------------------------

g1 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/g1-population-level.RDS")
sigma <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/sigma-population-year-level-mat.RDS")
fec <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/combinedF-population-year-level-mat.RDS")
phi <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/phi-population-year-level-mat.RDS")
perCapitaRS <- readRDS("outputs/005_calculatePopulationModelParameters/04_reproductiveSuccess/reproductiveSuccessWithCorrectionForMissingness-populationYear-mat.RDS")

# get mode from full posterior and calculate var in RS based on modes
sigma.mode <- lapply(sigma,apply,2,posterior.mode)
fec.mode <- lapply(fec,apply,2,posterior.mode)
phi.mode <- lapply(phi,apply,2,posterior.mode)
perCapitaRS.mode <- lapply(perCapitaRS,apply,2,posterior.mode)

# -------------------------------------------------------------------
# Compare calculation of geometric SD in reproductive success
# -------------------------------------------------------------------
par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

plot(NA,NA,xlim=c(0,10),ylim=c(0,10),pch=19);

abline(a=0,b=1,col='black',lty='dotted')

for(i in 1:20){
  text(gsd.am(sigma.mode[[i]]*fec.mode[[i]]*phi.mode[[i]]),gsd.am(perCapitaRS.mode[[i]]),i);
}

mtext(side=1,"GSD individual estimates",line=1.5)
mtext(side=2,"GSD reproductive success",line=1.5)

# -------------------------------------------------------------------
# Calculate geometric SD in reproductive success
# -------------------------------------------------------------------

df.list = list()
for(i in 1:20){
  obj = perCapitaRS[[i]]
  tmp.df=apply(obj,1,gsd.am)
  df.list[[i]] = tmp.df
}

gsdSummary=do.call(cbind,df.list)

# -------------------------------------------------------------------
# Analyze correlation of germination and GSD per-capita RS
# -------------------------------------------------------------------

probability.g1 = g1

# - +Summary statistics on germination ----
# - ++Calculate the 95% credible interval ----
CI.g1 <- apply(probability.g1,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# - ++Calculate the 68% highest posterior density interval ----
HPDI.g1 <- apply(probability.g1,2,FUN = function(x) hdi(x, .68))
# - ++Calculate the posterior mode ----
mode.g1 <- apply(probability.g1,2, FUN = posterior.mode)

# - ++Construct data frame ----
# use highest posterior density interval and mode
g1PosteriorSummary <- data.frame(cbind(t(HPDI.g1),mode.g1))
names(g1PosteriorSummary) <- c("lo.g1","hi.g1","mode.g1")

# - +Summary statistics on variability in RS ----
# - ++Calculate the 95% credible interval ----
CI.rs <- apply(gsdSummary,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# - ++Calculate the 68% highest posterior density interval ----
HPDI.rs <- apply(gsdSummary,2,FUN = function(x) hdi(x, .68))
# - ++Calculate the posterior mode ----
mode.rs <- apply(gsdSummary,2, FUN = posterior.mode)

# - ++Construct data frame ----
rsPosteriorSummary<-data.frame(cbind(t(HPDI.rs),mode.rs))
names(rsPosteriorSummary) <- c("lo.rs","hi.rs","mode.rs")

# create empty vector for the correlation
posterior.correlation<-c()

n.iter = dim(gsdSummary)[1]
# calculate correlation for each draw from the posterior
for(i in 1:n.iter){
  posterior.correlation[i]<-cor(probability.g1[i,],gsdSummary[i,])
}

# - +Summary statistics on correlation ----
# - ++Calculate the 95% credible interval ----
CI.correlation <- quantile(posterior.correlation, c(.025, .5, .975))
# - ++Calculate the 95% highest posterior density interval ----
HPDI.correlation <- hdi(posterior.correlation, .95)
# - ++Calculate the posterior mode ----
mode.correlation <- posterior.mode(posterior.correlation)

# - ++Construct data frame ----
correlationPosteriorSummary<-data.frame(cbind(t(HPDI.correlation),mode.correlation))
names(correlationPosteriorSummary) <- c("lo.corr","hi.corr","mode.corr")

signif(correlationPosteriorSummary,3)

# - Manuscript figure ----

# - +set font sizes ----
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12
pt5 = 5/12

dev.off()

tiff(filename=paste0("products/figures/correlationGerminationVariabilityRS.tif"),
     height=3.2,width=3.2,units="in",res=800,compression="lzw",pointsize=12)

#par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(2,2,.7,0)+.1,mgp=c(3,.45,0))
par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(2,2.2,.7,0)+.1,mgp=c(3,.45,0))
plot(x = NA,
     y = NA,
     xlim=c(0,10),ylim=c(0,.425),
     pch=16, cex = 0.5,
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt8)

# for plotting site names
# d.plot=data.frame(s=rsPosteriorSummary$mode.rs,
#                   g=g1PosteriorSummary$mode.g1,
#                   site=siteNames)

segments(x0=rsPosteriorSummary$lo.rs,x1=rsPosteriorSummary$hi.rs,
         y0=g1PosteriorSummary$mode.g1, lwd=1)
segments(x0=rsPosteriorSummary$mode.rs,
         y0=g1PosteriorSummary$lo.g1, y1=g1PosteriorSummary$hi.g1,
         lwd=1)
points(rsPosteriorSummary$mode.rs,g1PosteriorSummary$mode.g1,
       pch=21,cex=1.5,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 1),lwd=0)
points(rsPosteriorSummary$mode.rs,g1PosteriorSummary$mode.g1,
       pch=21,col='black',cex=1.5,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 0.5))

axis(1, seq(0,10,by=2), padj = -.5,
     labels = seq(0,10,by=2), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(1, seq(0,10,by=1),labels=FALSE)
axis(2, seq(0,1,by=.1),
     labels = seq(0,1,by=.1), las = 1, line = 0, hadj= 1.2,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(2, seq(.05,1,by=.1),labels=FALSE)

mtext("Germination probability",
      side=2,line=1.5,adj=.5,col='black',cex=pt10)
mtext("Geometric SD of per-capita reproductive success",
      side=1,line=1,adj=-.025,col='black',cex=pt10,at=-2)

box()
legend("topleft",bty='n',inset=c(-.05,0),
       paste0("Pearson's r=",round(correlationPosteriorSummary[3],3),
              " (",round(correlationPosteriorSummary[1],3),", ",
              round(correlationPosteriorSummary[2],3),")"),
       cex=pt7)

mtext("B.", adj = 0, cex=pt10)

dev.off()
