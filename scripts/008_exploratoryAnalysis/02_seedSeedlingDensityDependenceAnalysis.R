#################################################################################
# Script to conduct exploratory analysis of density-dependence in seed to seedling transition
# we conducted the analysis in this script after the first round of peer review
# a reviewer suggested we consider whether there was density-dependence in seedling establishment
# we followed methods from:
# Detto, M., M. D. Visser, S. J. Wright, and S. W. Pacala. 2019. Bias in the detection of negative density dependence in plant communities. Ecology Letters 22:1923–1939.
# and use code provided with the paper to conduct the analysis 
# the code is located in the Github repository:
# https://github.com/mdetto/Bias-in-the-detection-of-negative-density-dependence
#################################################################################

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

# - Libraries ----
library(tidybayes)
library(tidyverse)
library(stringr)
library(khroma)

# - Source functions for analysis ----
source("scripts/006_testHypotheses/00_utilityFunctions.R")

# - Set file path ----
dataDirectory = "data/"

# - Read in data ----
# counts of seedlings and fruiting plants in permanent plots
seedlingFruitingPlantCountsPermanentPlots <- read.csv(paste0(dataDirectory,"seedlingFruitingPlantCountsPermanentPlots.csv"),header=TRUE)
# counts of total fruit equivalents per plant in permanent plots (2006-2012)
countFruitsPerPlantFromPermanentPlots <- read.csv(paste0(dataDirectory,"countFruitsPerPlantFromPermanentPlots.csv"),header=TRUE)
# counts of undamaged and damaged fruits per plant in permanent plots (2013-2020)
countUndamagedDamagedFruitsPerPlantFromPermanentPlots <- read.csv(paste0(dataDirectory,"countUndamagedDamagedFruitsPerPlantFromPermanentPlots.csv"),header=TRUE)

# site names
siteAbiotic <- read.csv(paste0(dataDirectory,"siteAbioticData.csv"),header=TRUE)
siteDf <- data.frame(site=1:20,siteName=siteAbiotic$site)

# - Build data frame for regression ----
# the goal is to run a regression of number of seedlings in a permanent plot on 
# the number of seeds in the seed bank in that permanent plot
# we observe the number of seedlings in annual surveys 
# we DO NOT have observations of seed bank size in permanent plots so here's what we do
# 1. we calculate seed rain in the permanent plots by multiplying the total number of fruits per plot
# with the estimated number of seeds per fruit
# assumptions we make in the analysis:
# a. in some pops/years/plots, we counted fruits on ALL plants; in others, we only have fruit counts on a subset
# when we have fruit counts, we use the observed numbers
# when we don't have fruit counts for some plants, we substitute the estimated number of fruits per plant for those plants
# b. estimates for fruits per plant and seeds per fruit come from the statistical models we describe in the main text
# in other words these are population-level estimates
# 2. we use estimates of seed survival and germination from the seed bag burial experiment
# to calculate the probability of a seed remaining in the seed bank in the 1st, 2nd, and 3rd year after seed set
# 3. we combine estimates of seed rain and transition probabilities to 
# calculate the expected size of the seed bank in each plot using seed production from t-1, t-2, t-3
# because we use a 3 year seed bank, and we only have data on plot-level fruit per plant counts starting in 2007,
# we calculate the expected size of the seed bank in 2010-2020, giving us 10 years of data

# - Build dataset for fruits per plant, per plot, population and year ----

# - +Read in year-level estimates of fruit numbers ----
fecEstimate.popYear <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/combinedF-population-year-level-mat.RDS")

# empty list
df.list <- list()
# for each site, calculate the posterior mode of total fruit equivalents per plant estimate
for(i in 1:20){
  tmp <- fecEstimate.popYear[[i]][,1:15]
  df.list[[i]] <- data.frame(site=i,
                             year = c(2006:2020),
                             mode = apply(tmp,2,posterior.mode))
}
# bind to data frame
fecPopYear.summary <- do.call(rbind,df.list)

# add site names to the data frame
fecPopYear.summary <- fecPopYear.summary %>% 
  dplyr::left_join(siteDf,by=c('site')) %>%
  dplyr::select(-site) %>% 
  dplyr::rename(site=siteName)

# - +Calculate number of plants missing observations of fruits/plant ----
# those missing observations will be replaced with estimates of fruits/plant

# tally observations of fruits per plant and find the difference between 
# observations of fruits/plant and plant counts/plot 
# difference is the missing number of observations of fruits/plant
countFruitsPerPlantFromPermanentPlots.sum <- countFruitsPerPlantFromPermanentPlots %>%
  dplyr::group_by(site,year,transect,position) %>%
  dplyr::summarise(nPlants = n(), totalFruits = sum(countFruitsPerPlant)) %>%
  dplyr::left_join(seedlingFruitingPlantCountsPermanentPlots,by=c('site','year','transect','position')) %>%
  dplyr::mutate(delta=nPlants-fruitplNumber)

# - +Use posterior mode to get net fruits/plant (combine field estimate plus estimate)  ----

# - ++ 2006-2012  ----

# for plots with missing observations of fruits/plant
# join the data frame with estimates of the posterior mode for fruits/plant
# and use that estimate as a substitute to calculate the total number of fruits
# in the plot
netCountFruitsPerPlantFromPermanentPlots.sum <- countFruitsPerPlantFromPermanentPlots.sum %>%
  dplyr::filter(delta<0) %>%
  dplyr::left_join(fecPopYear.summary,by=c('site','year')) %>%
  dplyr::mutate(netFruits = totalFruits + (-1*delta)*mode) %>%
  dplyr::select(site,year,transect,position,nPlants,seedlingNumber,fruitplNumber,netFruits,delta) %>%
  dplyr::rename(totalFruits=netFruits)

# join the data with complete counts of fruit/plants
# with data that combines field counts and estimates
totalFruitEquivalentNumber <- countFruitsPerPlantFromPermanentPlots.sum %>%
  dplyr::filter(delta==0) %>%
  dplyr::bind_rows(netCountFruitsPerPlantFromPermanentPlots.sum) 

# - ++ 2013-2013  ----

# tally observations of undamaged/damaged fruits per plant and find the difference between 
# observations of fruits/plant and plant counts/plot 
# difference is the missing number of observations of fruits/plant
countUndamagedDamagedFruitsPerPlantFromPermanentPlots.sum <- countUndamagedDamagedFruitsPerPlantFromPermanentPlots %>%
  dplyr::group_by(site,year,transect,position) %>%
  dplyr::summarise(nPlants = n(), 
                   totalUndamagedFruits = sum(countUndamagedFruitsPerPlant),
                   totalDamagedFruits = sum(countDamagedFruitsPerPlant)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(seedlingFruitingPlantCountsPermanentPlots,by=c('site','year','transect','position')) %>%
  dplyr::mutate(delta=nPlants-fruitplNumber)

# for plots with missing observations of undamaged/damaged fruits/plant
# join the data frame with estimates of the posterior mode for fruits/plant
# and use that estimate as a substitute to calculate the total number of fruits
# in the plot
netCountUndamagedDamagedFruitsPerPlantFromPermanentPlots.sum <- countUndamagedDamagedFruitsPerPlantFromPermanentPlots.sum %>%
  dplyr::filter(delta<0) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(fecPopYear.summary,by=c('site','year')) %>%
  dplyr::mutate(netFruits =  (-1*delta)*mode) %>%
  dplyr::select(site,year,transect,position,nPlants,seedlingNumber,fruitplNumber,totalUndamagedFruits,totalDamagedFruits,netFruits,delta)

# - Build dataset for seeds per fruit, per population and year ----

# - +Read in posterior estimates ----

phi <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/phi-population-year-level-mat.RDS")
phiDamage <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/phi-damage-population-year-level-mat.RDS")

# empty list
df.list <- df.list2 <- list()
# for each site, calculate the posterior mode of total fruit equivalents per plant estimate
for(i in 1:20){
  tmp <- phi[[i]][,1:15]
  df.list[[i]] <- data.frame(site=i,
                             year = c(2006:2020),
                             mode = apply(tmp,2,posterior.mode))
  tmp2 <- phiDamage[[i]][,1:15]
  df.list2[[i]] <- data.frame(site=i,
                              year = c(2006:2020),
                              mode.dam = apply(boot::inv.logit(tmp2),2,posterior.mode))
  
}

# seeds per fruit
# bind to data frame
phiPopYear.summary <- do.call(rbind,df.list) %>%
  dplyr::select(site,year,mode)
# add site names to the data frame
phiPopYear.summary <- phiPopYear.summary %>% 
  dplyr::left_join(siteDf,by=c('site')) %>%
  dplyr::select(-site) %>% 
  dplyr::rename(site=siteName) 

# proportion of seeds lost to herbivory
# bind to data frame
phiDamagePopYear.summary <- do.call(rbind,df.list2) %>%
  dplyr::select(site,year,mode.dam)
# add site names to the data frame
phiDamagePopYear.summary <- phiDamagePopYear.summary %>% 
  dplyr::left_join(siteDf,by=c('site')) %>%
  dplyr::select(-site) %>% 
  dplyr::rename(site=siteName) 

# - Join fruit and seed data ----
# - + 2006-2012 ----

tmp.df <- totalFruitEquivalentNumber %>%
  dplyr::left_join(phiPopYear.summary,by=c('site','year')) %>%
  dplyr::mutate(seedRain=totalFruits*mode)

seedRain1 <- tmp.df %>%
  dplyr::ungroup() %>%
  dplyr::select(site,year,transect,position,seedRain)

# - + 2013-2020 ----

# calculate seed rain for plots with observations of fruit numbers or all plants
seedRain2 <- countUndamagedDamagedFruitsPerPlantFromPermanentPlots.sum %>%
  dplyr::filter(delta==0) %>%
  dplyr::left_join(phiPopYear.summary,by=c('site','year')) %>%
  dplyr::left_join(phiDamagePopYear.summary,by=c('site','year')) %>%
  dplyr::mutate(seedRain = totalUndamagedFruits*mode+totalDamagedFruits*mode*mode.dam) %>%
  dplyr::ungroup() %>%
  dplyr::select(site,year,transect,position,seedRain)

# calculate seed rain for plots with missing observations of fruit numbers or all plants
seedRain3 <- netCountUndamagedDamagedFruitsPerPlantFromPermanentPlots.sum %>%
  dplyr::left_join(phiPopYear.summary,by=c('site','year')) %>%
  dplyr::left_join(phiDamagePopYear.summary,by=c('site','year')) %>%
  dplyr::mutate(seedRain = (totalUndamagedFruits+netFruits)*mode+totalDamagedFruits*mode*mode.dam) %>%
  dplyr::ungroup() %>%
  dplyr::select(site,year,transect,position,seedRain)

# join seed rain for all plots from 2013-2020
seedRain4 <- seedRain2 %>%
  bind_rows(seedRain3)

# join seed rain for all plots from 2006-2020
seedRain <- seedRain1 %>%
  bind_rows(seedRain4)

# - Create data frame for 3-year seed rain history for 2010-2020 ----

# vector of years
years <- 2010:2020

# empty list
obs.list <- list()

# for each year
for(i in 1:10){
  
  # use past 3 years of data
  current.year <- years[i]-c(1,2,3)
  
  # get seed rain 1 year ago
  t1 <- seedRain %>%
    dplyr::filter(year==current.year[1]) %>%
    dplyr::select(-c(year)) %>%
    dplyr::mutate(seedRain=ifelse(is.na(seedRain),0,seedRain)) %>%
    dplyr::rename(y_1=seedRain) 
  
  # get seed rain 2 years ago
  t2 <- seedRain %>%
    dplyr::filter(year==current.year[2]) %>%
    dplyr::select(-year) %>%
    dplyr::mutate(seedRain=ifelse(is.na(seedRain),0,seedRain)) %>%
    dplyr::rename(y_2=seedRain)
  
  # get seed rain 3 years ago
  t3 <- seedRain %>%
    dplyr::filter(year==current.year[3]) %>%
    dplyr::select(-year) %>%
    dplyr::mutate(seedRain=ifelse(is.na(seedRain),0,seedRain)) %>%
    dplyr::rename(y_3=seedRain) 
  
  # build data frame of seed rain history for last 3 years
  seedRainHistory <- t1 %>%
    dplyr::full_join(t2,by=c("site","transect","position")) %>%
    dplyr::full_join(t3,by=c("site","transect","position")) %>%
    dplyr::mutate(year=years[i]) %>%
    dplyr::mutate(y_3=ifelse(is.na(y_3),0,y_3)) %>%
    dplyr::mutate(y_2=ifelse(is.na(y_2),0,y_2)) %>%
    dplyr::mutate(y_1=ifelse(is.na(y_1),0,y_1))
  
  # get seedling observations 
  seedlingObservations <- seedlingFruitingPlantCountsPermanentPlots %>% 
    dplyr::select(site,year,transect,position,seedlingNumber)
  
  # join the seedling observations for current year to seed rain history for previous 3 years
  seedRainSeedlingHistory <- seedRainHistory %>%
    dplyr::left_join(seedlingObservations,by=c('site','transect','position','year')) 
  
  # put into list
  obs.list[[i]]<-seedRainSeedlingHistory  
}

# bind list into data frame
seedHistory<-do.call(rbind,obs.list)

# - Calculate expected seed bank size ----
# Combine seed rain and seed survival/germination estimates

# - +Get estimates for seed survival and germination ----
s0 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s0-population-level.RDS")
g1 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/g1-population-level.RDS")
s1 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s1-population-level.RDS")
s2 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s2-population-level.RDS")
s3 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s3-population-level.RDS")

# calculate posterior mode
g1.hat <- apply(g1,2,posterior.mode)
s1.hat <- apply(s1,2,posterior.mode)
s2.hat <- apply(s2,2,posterior.mode)
s3.hat <- apply(s3,2,posterior.mode)
s0.hat <- apply(s0,2,posterior.mode)

# build data frame
vr.dat <- data.frame(site=siteAbiotic$site,g1=g1.hat,s0=s0.hat,s1=s1.hat,s2=s2.hat,s3=s3.hat)

# - +Calculate transition probabilities ----
seedBankTransitionProbabilities <- vr.dat %>%
  dplyr::mutate(transitionOne=s0*s1,
                transitionTwo=s0*s1*(1-g1)*s2*s3,
                transitionThree=s0*s1*(1-g1)*s2*s3*(1-g1)*s2*s3) %>%
  dplyr::select(site,transitionOne,transitionTwo,transitionThree,g1)

# - +Join seed rain history and transition probabilities ----
seedHistoryTransitionDataFrame <- seedHistory %>%
  dplyr::left_join(seedBankTransitionProbabilities,by='site')

# - +Calculate expected seed bank size ----
expectedSeedBankSizeDataFrame <- seedHistoryTransitionDataFrame %>%
  dplyr::mutate(expectSeedBank=y_1*transitionOne+y_2*transitionTwo+y_3*transitionThree) %>%
  dplyr::mutate(expectSeedlings=expectSeedBank*g1)

# - Exploratory plot ----

ggplot(data=expectedSeedBankSizeDataFrame ,aes(x=expectSeedBank,y=seedlingNumber,color=as.factor(site))) +
  geom_point() +
  facet_wrap(~site,scales='free') +
  geom_abline(intercept=0,slope=1) +
  scale_color_batlow(reverse=TRUE,discrete=TRUE) +
  xlab("Expected seeds in 3-year seed bank (January year t)") +
  ylab("Observed seedlings (year t)") +
  labs(color="Population") +
  theme_bw() +
  theme(panel.spacing=unit(1,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(legend.background = element_rect(fill="gray95"),
        legend.key = element_rect(fill="gray95")) 

# - Fit offset power-law model ----

# build data frame for regression
regressionDataFrame <- expectedSeedBankSizeDataFrame %>%
  # select relevant variables
  dplyr::select(site,year,transect,position,expectSeedBank,seedlingNumber) %>%
  # round the expected seed bank size to a whole number
  dplyr::mutate(seedBank=round(expectSeedBank)) %>%
  # if number of seedlings is greater than expected seed bank size, set the seed bank size to the number of seedlings
  dplyr::mutate(seedBank=ifelse(seedBank<seedlingNumber,seedlingNumber,seedBank)) 

# get a vector of site names
siteNames <- siteAbiotic$site

# empty list
out.list <- list()

# for each site, fit the power-offset model 
# regressing seedlings on the seed bank size
for(i in 1:20){
  df.tmp <- regressionDataFrame %>%
    dplyr::filter(site==siteNames[i])
  m1 <- lm(log(seedlingNumber+1)~log(seedBank+1),data=df.tmp)
  # get slope coefficient and confidence intervals
  out.list[[i]] <- c(coef(m1)[2],confint(m1,'log(seedBank + 1)',level=0.95))
}

# bind model coefficients into matrix
out.mat<-do.call(rbind,out.list)

# - Fit offset power-law model with correction for offset  ----
# code to run the correction is modified from the repository for Detto et al. 2019
# https://github.com/mdetto/Bias-in-the-detection-of-negative-density-dependence

# the function rbinom2 is located 
# https://github.com/mdetto/Bias-in-the-detection-of-negative-density-dependence/blob/56c0537b40201c080e33ef23ef93ef204a5bdd21/rbinom2.r
rbinom2 <- function(n,S,p){
  
  # random number generator of S trials given R number success with Bernoulli
  # probability p
  z <- runif(n)
  R <- numeric(length(S))
  
  use=S==0
  R[use] <- 0
  
  S0 <- unique(S[S>0])
  glS1 <- lgamma(S0+1)
  if (p>0 & p<1){
    for (i in 1:length(S0)){
      
      use2 <- S==S0[i]
      k <- seq.int(0,S0[i])
      C <- cumsum(exp(glS1[i] - lgamma(k+1) -lgamma(S0[i]-k+1) + k*log(p) + (S0[i]-k)*log(1-p)))
      R[use2] <- .bincode(z[use2], c(0,C))-1
      
    }
  } else {
    R <- S*p
  }
  
  return(R)
}

# create empty data frame to hold results of fitting model
H1 <-  data.frame(b1=numeric(length(20)),
                  lb1=numeric(length(20)),
                  ub1=numeric(length(20)))

# create empty vector to hold expected probability of survival
f <- numeric()

# fit model with correction for the offset for all populations
for(i in 1:20){
  
  # select site
  df.tmp <- regressionDataFrame %>%
    dplyr::filter(site==siteNames[i])
  
  # get the average probability of emergence
  f[i] <- sum(df.tmp$seedlingNumber)/sum(df.tmp$seedBank)
  
  Rep <- 1000
  mod.rnd1 <- numeric()
  
  # the code to fit the model with the correction for the offset is at
  # lines 149-155 of the script used by Detto et al. 2019 to reanalyze data from Harms et al. 2000
  # https://github.com/mdetto/Bias-in-the-detection-of-negative-density-dependence/blob/56c0537b40201c080e33ef23ef93ef204a5bdd21/reanalysisHarmsEtAl2000_v3.1.R
  for(j in 1:Rep){
    cat("\r",i,":",Rep-j,"\r")
    
    # null model 1 
    # rerun the regression on data simulated using the site-level average probability of seedling emergence
    
    # simulate data
    R.r = rbinom2(length(df.tmp$seedBank),df.tmp$seedBank,f[i])
    seedlingNumber.sim = R.r  
    seedBank.tmp = df.tmp$seedBank
    
    # fit model with observed seed bank size (seedBank.tmp) and simulated seedling number (seedlingNumber.sim)
    mod.rnd1[j] <- coef(lm(log(seedlingNumber.sim+1)~log(seedBank.tmp+1)))[2]
    
  }
  
  # get model coefficients
  H1[i,] <- c(mean(mod.rnd1,na.rm=TRUE),quantile(mod.rnd1,prob=c(0.025,0.975),na.rm=TRUE))
  
}

# - Plot figure  ----

pdf(file=paste0("products/figures/densityDependence-seedSeedling.pdf"),
    height=6,width=6,pointsize=12)

par(mar = c(6.1, 4.6, 4.1, 4.1), # change the margins
    lwd = 1, # increase the line thickness
    cex.axis = 1.2,
    xpd=F# increase default axis label size
)

index=order(siteAbiotic$easting)
plot(1:20-.1,out.mat[index,1]-1,ylim=c(-1.2,0),pch=NA,
     xaxt="n",xlab="",ylab=expression(italic(hat(b)) - italic(b)['NULL']))

rect(xleft=c(1,3,5,7,9,11,13,15,17,19)+c(.5),xright=c(1,3,5,7,9,11,13,15,17,19)+1.5,ybottom=-100,ytop=100,col='gray99',border='gray99')

points(1:20-.1,out.mat[index,1]-1,ylim=c(-1.2,0),pch=16)
abline(h=0,lty='dotted')
segments(x0=1:20-.1,y0=out.mat[index,2]-1,y1=out.mat[index,3]-1)
points(x=1:20+.1,y=out.mat[index,1]-H1$b1[index],pch=16,cex=1,col='red')
segments(y0=out.mat[index,1]-H1$ub1[index],x0=1:20+.1,y1=out.mat[index,1]-H1$lb1[index],col='red')
axis(1,at=1:20,labels=siteNames[index],las=2)
box()

## Draw the x-axis labels.
par(xpd=T)
legend(0.05, 0.25,
       c("Estimate from offset power-law model", 
         "Estimate from offset power-law model with correction for offset"),
       col = c("black", "red"),
       cex = 0.8,
       pch=16,lty=1)

dev.off()

