# -------------------------------------------------------------------
# Density-independent model of germination
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
library(bayesplot

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

s0 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s0-population-level.RDS")
g1 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/g1-population-level.RDS")
s1 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s1-population-level.RDS")
s2 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s2-population-level.RDS")
s3 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s3-population-level.RDS")
rs <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/reproductiveSuccess-population-year-level-mat.RDS")


# get mode from full posterior and calculate var in RS based on modes
g1.hat  <- apply(g1,2,posterior.mode)
s1.hat  <- apply(s1,2,posterior.mode)
s2.hat  <- apply(s2,2,posterior.mode)
s3.hat  <- apply(s3,2,posterior.mode)
s0.hat  <- apply(s0,2,posterior.mode)
rs.hat <- lapply(rs,apply,2,posterior.mode)

# -------------------------------------------------------------------
# calculate autocorrelation of per-capita reproductive success (mode)
# no significant autocorrelation, include in appendix
# -------------------------------------------------------------------

allPopulations = order((position %>% dplyr::select(site,easting))$easting)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for( k in allPopulations){
  rs.tmp = rs.hat[[k]]
  yt = rs.tmp
  acf(yt,na.action=na.pass)
  text(5,.9,(position %>% dplyr::select(site,easting))$site[k])
}

# -------------------------------------------------------------------
# calculate harmonic mean for each population
# -------------------------------------------------------------------

# function for harmonic mean
f = function(x){1/mean(1/x,na.rm=TRUE)}

pop.harmonicMean =c ()
for( k in 1:20){
 # pop.index=index[,1]==k
  rs.tmp = c(rs.hat[[k]])
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

for( k in 1:20){
  # - ++get the population index  ----
 # pop.index=index[,1]==k
  pop.index = k

  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = rs.hat[[pop.index]]
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

# -------------------------------------------------------------------
# Make plots
# -------------------------------------------------------------------

# - ++save manuscript data  ----

optimalGerminationFractions = data.frame(site=siteNames,optima=optima,g1.hat=g1.hat)

outputDirectory <- "outputs/006_hypothesisTesting/"
saveRDS(optimalGerminationFractions,paste0(outputDirectory,"optimalGerminationFractions.RDS"))

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

plot(NA,NA,xlim=c(0,1),ylim=c(0,1),
     xlab = "",
     ylab = "",bty='n',
     xaxt= "n", yaxt="n",)

points(x=optima,y=g1.hat,pch=21,col='black',bg='white',cex=2.5)

# Now, define a custom axis
d.plot=data.frame(optima,
                  g1.hat,site=siteNames)

text(d.plot[,1:2],siteNames,cex=4/12)

box()


abline(a=0,b=1,lty='dotted')

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

mtext("A.", adj = 0, cex=pt10)

dev.off()
