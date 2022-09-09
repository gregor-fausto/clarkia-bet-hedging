# -------------------------------------------------------------------
# Density-independent model of germination
# Figure output IS NOT USED in the manuscript
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)

# - Source functions for analysis ----
source("scripts/006_testHypotheses/00_utilityFunctions.R")

# - Function to calculate one-step population growth rate based on vital rates ----
fitness <- function(g=g1,s0=s0,s1=s1,s2=s2,s3=s3,rs=rs){
  p1 = g*rs*s0*s1
  p2 = (1-g)*(s2*s3)
  return(as.numeric(p1+p2))
}

# - Libraries ----
library(rjags) # jags interface
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(HDInterval)
library(bayesplot)

# - Read in site abiotic data ----

siteAbiotic <- read.csv("data/siteAbioticData.csv",header=TRUE)

position<-siteAbiotic %>%
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteNames = unique(position$site)

# - Create index of site numbers and years ----

siteIndex <- data.frame(site=siteNames,siteIndex=1:20)
yearIndex <- data.frame(year=2006:2020,yearIndex=1:15)
index=expand.grid(1:20,1:15)

# - Read in samples from posterior distribution ----

s0 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s0-population-level.RDS")
g1 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/g1-population-level.RDS")
s1 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s1-population-level.RDS")
s2 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s2-population-level.RDS")
s3 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s3-population-level.RDS")
perCapitaRS <- readRDS("outputs/005_calculatePopulationModelParameters/04_reproductiveSuccess/reproductiveSuccessWithCorrectionForMissingness-populationYear-mat.RDS")

# get mode from full posterior and calculate var in RS based on modes
g1.hat  <- apply(g1,2,posterior.mode)
s1.hat  <- apply(s1,2,posterior.mode)
s2.hat  <- apply(s2,2,posterior.mode)
s3.hat  <- apply(s3,2,posterior.mode)
s0.hat  <- apply(s0,2,posterior.mode) 
perCapitaRS.hat <- lapply(perCapitaRS,apply,2,posterior.mode)

# -------------------------------------------------------------------
# calculate autocorrelation of per-capita reproductive success (mode)
# CP3 is the only population showing significant autocorrelation at lag>1
# which is within expected number of observations by chance
# -------------------------------------------------------------------

allPopulations = order((position %>% dplyr::select(site,easting))$easting)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for( k in allPopulations){
  rs.tmp = perCapitaRS.hat[[k]]
  yt = rs.tmp
  acf(yt,na.action=na.pass)
  text(5,.9,(position %>% dplyr::select(site,easting))$site[k])
}

# -------------------------------------------------------------------
# calculate harmonic mean for each population
# -------------------------------------------------------------------

# function for harmonic mean
harmonicMean = function(x){1/mean(1/x,na.rm=TRUE)}

pop.harmonicMean =c ()
for( k in 1:20){
  # pop.index=index[,1]==k
  rs.tmp = c(perCapitaRS.hat[[k]])
  hm = harmonicMean(rs.tmp)
  pop.harmonicMean[k] = hm
}

# plot harmonic mean of per-capita reproductive success vs. seed survival
par(mfrow=c(1,1))
plot(pop.harmonicMean,s2.hat*s3.hat,type='n',xlim=c(0,20),ylim=c(-.1,1));
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
n.iter=dim(perCapitaRS)[1]

# -------------------------------------------------------------------
# for each population k (use mode)
# get the median of per-capita reproductive success in each year
# sample a sequence of 1000 years to use in the analysis for each value of g
# for each g, calculate the growth rate for 1000 years
# -------------------------------------------------------------------

for( k in 1:20){
  # - ++get the population index  ----
  # pop.index=index[,1]==k
  pop.index = k

  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = perCapitaRS.hat[[pop.index]]
  y_t = y_t[!is.na(y_t)]

  # - ++draw 10000 samples for reproductive success with replacement  ----
  y_t.resample = sample(y_t,10000,replace=TRUE)
  
  for( i in 1:length(g)){
    fit<-fitness(g=g[i],s0.hat[k],s1.hat[k],s2.hat[k],s3.hat[k],rs=y_t.resample)
    #logfit<-log(fit)
    g.mat[,i] <- fit
    #g.mat[,i] <- logfit
  }
  g.sites[[k]] <- g.mat
}

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
  mtext(siteNames[k],side=3,line=-2)
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

# - Save manuscript data  ----

optimalGerminationFractions = data.frame(site=siteNames,optima=optima,g1.hat=g1.hat)

outputDirectory <- "outputs/006_hypothesisTesting/"
write.csv(optimalGerminationFractions,paste0(outputDirectory,"optimalGerminationFractions.csv"),row.names=FALSE)

# - PLOT  ----

pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12
pt5 = 5/12

tiff(filename=paste0("products/figures/optimalGerminationFraction.tif"),
     height=3.75,width=3.5,units="in",res=300,compression="lzw",pointsize=12)

par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(2,2.2,.7,1)+.1,mgp=c(3,.45,0))

plot(NA,NA,xlim=c(.4,1.04),ylim=c(0,.4),
     xlab = "",
     ylab = "",bty='n',
     xaxt= "n", yaxt="n",)

#points(x=optima,y=g1.hat,pch=21,col='black',bg='white',cex=2.5)

points(x=optima,y=g1.hat,
       pch=21,cex=2.5,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 1),lwd=0)
points(x=optima,y=g1.hat,lwd = 0.5,
       pch=21,col='black',cex=2.5,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 0.5))


# Now, define a custom axis
d.plot=data.frame(optima,
                  g1.hat,site=siteNames)


# # Now, define a custom axis

d.plot[c(4,5,6),1]=1.04
d.plot[c(3,12,13),2]=c(.1665,.14,0.12)
d.plot[18,c(1:2)]=d.plot[18,c(1:2)]*c(.995,1.01)

text(d.plot[,1:2],siteNames,cex=4/12)

box()

#abline(a=0,b=1,lty='dotted')

axis(1, seq(0,1,by=.1), padj = -.5,
     labels = seq(0,1,by=.1), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(1, seq(.25,1,by=.1),labels=FALSE)
axis(2, seq(0,1,by=.1),
     labels = seq(0,1,by=.1), las = 1, line = 0, hadj= 1.2,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(2, seq(.05,1,by=.1),labels=FALSE)

mtext("Observed germination fraction",
      side=2,line=1.5,adj=.5,col='black',cex=pt10)
mtext("Predicted germination fraction",
      side=1,line=1,adj=.5,col='black',cex=pt10)

# uncomment if more than one panel
# mtext("A.", adj = 0, cex=pt10)

dev.off()
