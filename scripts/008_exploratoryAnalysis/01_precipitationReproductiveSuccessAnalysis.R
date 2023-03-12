# -------------------------------------------------------------------
# Script to conduct analyze relationship between 
# spring precipitation and per-capita reproductive success
# analysis follows that in Venable 2007
# https://doi.org/10.1890/06-1495
# -------------------------------------------------------------------

rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)

# - Libraries ----
library(rjags) # jags interface
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(HDInterval)
library(bayesplot)

# - Source functions for analysis ----
source("scripts/006_testHypotheses/00_utilityFunctions.R")

# - Site names ----
position<-read.csv(file="~/Dropbox/clarkia-demography-projects/data/siteAbioticData.csv",header=TRUE) %>%
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteNames <- unique(position$site)

# - Read in posterior distributions of per-capita reproductive success ----

perCapitaRS <- readRDS("outputs/005_calculatePopulationModelParameters/04_reproductiveSuccess/reproductiveSuccessWithCorrectionForMissingness-populationYear-mat.RDS")

# - Read in climate data ----
climate <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData-2021.RDS")
climate <- climate[climate$site%in%siteNames,] %>%
  # specifically the data on spring precipitation
  dplyr::filter(season=="spring"&variable=='p') %>%
  dplyr::select(site,year,value) %>%
  dplyr::rename(springPrecipitation=value)
write.csv(climate,file="~/Downloads/springPrecipitationData.csv")
# - Create climate+RS data frames ----

# empty list
df.list=list()

# for each population
for(i in 1:20){
  # get the climate data for site i
  tmp.clim <- climate[climate$site==siteNames[i],] %>%
    # specifically the data on spring precipitation
    dplyr::filter(season=="spring"&variable=='p')
  # calculate the posterior mode of per-capita RS for that population
  tmp <- apply(perCapitaRS[[i]],2,posterior.mode)
  # bind the data on spring precipitation and per-capita RS to create a new data frame
  df <- data.frame(tmp.clim,rs=tmp)
  # save data frames in a list
  df.list[[i]] <- df
}

# bind the data frames into a single data frame
df <- do.call(rbind,df.list)

# - Fit regression of RS on log of spring precipitation ----

# empty vector for slope parameters
b.vec <- c()

# for each population
for(i in 1:20){
  # get the data frame
  tmp <- df.list[[i]]
  # take the log of spring precipitation
  tmp$value.log <- log(tmp$value+1)
  # take the log of per-capita reproductive success
  tmp$rs.log <- log(tmp$rs+.5)
  # fit a regression of log per capita RS on log spring precipitation
  m1 <- lm((tmp$rs.log)~(tmp$value.log))
  # extract the slope parameter
  b <- coef(m1)[2];
  # put the slope parameter into a vector
  b.vec[i] <- b
}

# - Plot figure  ----

pdf("~/Dropbox/clarkia-bet-hedging/products/figures/rs-climate-sensitivity.pdf", height = 8, width = 8)

index=order(position$easting)
par(mfrow=c(4,5),mar=c(0,.5,.5,0),
    oma=c(4,4,1,1))

p.val <- p.val.all <- c()
for(i in 1:20){
  tmp <- df.list[[index[i]]]
  tmp$value.log <- log(tmp$value+1)
  tmp$rs.log <- log(tmp$rs+.5)
  m1 <- lm((tmp$rs.log)~(tmp$value.log))
  b <- coef(m1)[2];  a = coef(m1)[1]
  r.2 <- signif(summary(m1)$r.squared,2)
  plot((tmp$value.log),(tmp$rs.log),xlim=c(3,6),ylim=c(-4,7),
       xaxt='n',yaxt='n',pch=16)
  segments(x0=min(tmp$value.log),x1=max(tmp$value.log),
           y0=a+b*min(tmp$value.log),y1=a+b*max(tmp$value.log),
           lwd=2,lty='dashed')
  ifelse(i%in%c(1,6,11,16),axis(2L,las=1),NA)
  ifelse(i%in%c(16:20),axis(1L),NA)
  
  p.val.all[i] <- (summary(m1))$coefficients[,4][2]
  p.val[1] <- ifelse((summary(m1))$coefficients[,4][2]<.05/20,1,0)
  p.val <- as.numeric(p.val)
  p.val <- rep("*",sum(p.val))
  p.val <- ifelse(length(p.val)==0,"",p.val)
  text(3,-4, cex=1, bquote(R^2==.(r.2)~.(p.val)),adj=c(0,0))
  text(3,-1.5, cex=1,siteNames[index[i]] ,adj=c(0,0))
  text(3,-2.75, cex=1,paste0("b=",signif(b,2)) ,adj=c(0,0))
}

mtext("Log(per-capita reproductive success+0.5)", side = 2, outer = TRUE, line = 2)
mtext("Log(spring precipitation+1)", side = 1, outer = TRUE, line = 2.2)

dev.off()
