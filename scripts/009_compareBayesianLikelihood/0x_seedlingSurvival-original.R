
rm(list=(ls(all=TRUE))) # if using in source(script)

# - Libraries ----
library(tidyverse)
library(parallel)
library(stringr)
library(bbmle)
library(MCMCvis)
library(lme4)

# - +Read in posterior for seed bag trials & viabiility trials ----

mcmcDirectory = "outputs/002_fitStatisticalModels/mcmcSamples/"
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedlingSurvivalPosteriorSamples-Transect.RDS",mcmcSampleDirectory)]])
#mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedlingSurvivalSamplesSpatial.RDS",mcmcSampleDirectory)]])

# - +Read in site names ----
siteAbiotic <- read.csv("data/siteAbiotic.csv",header=TRUE)
siteNames <- siteAbiotic$site

# - Read in data to fit ML models ----

tmpDataDirectory = "outputs/001_prepareDataForModels/"
# 
seedlingData <- readRDS(file=paste0(tmpDataDirectory,"seedlingRef.RDS"))

modelDataDirectory = "outputs/002_fitStatisticalModels/data/"
modelDataDirectory <- paste0(modelDataDirectory,list.files(modelDataDirectory))
# seedlingData <- readRDS(modelDataDirectory[[grep("seedlingSurvivalData.RDS",modelDataDirectory)]])
seedlingData = data.frame(site=seedlingData$site,year=as.character(seedlingData$year),
                          transect=as.character(seedlingData$transect),
                          seedlingNumber=seedlingData$seedlingNumber,
                          fruitplNumber=seedlingData$fruitplNumber)

# - Set up model ----

glm.list <- glmer.list <- mean.list <- sampleSize.list <- list()

# - + fixed effects model ----

for(i in 1:20){
   site.tmp=siteNames[i]
 # site.tmp=i
  data.tmp <- seedlingData[seedlingData$site==site.tmp,]
  glm.tmp <-glm(cbind(fruitplNumber,seedlingNumber-fruitplNumber)~-1+year,data=data.tmp,family='binomial')
  glm.list[[i]] <- glm.tmp
}

# - + mixed effects effects model ----

for(i in 1:20){
   site.tmp=siteNames[i]
  # site.tmp=i
  data.tmp <- seedlingData[seedlingData$site==site.tmp,]
  glmer.tmp <-lme4::glmer(cbind(fruitplNumber,seedlingNumber-fruitplNumber)~1 + (1|year),data=data.tmp,family='binomial')
  glmer.list[[i]] <- glmer.tmp
}

var.comps<-lapply(glmer.list,f,var=TRUE)
var.comps.mat.1 <- do.call(rbind,var.comps)
alpha_y<-lapply(glmer.list,f,var=FALSE)
alpha.mat <- do.call(rbind,alpha_y)

# - + mixed effects effects model with transect ----

for(i in 1:20){
   site.tmp=siteNames[i]
  #site.tmp=i
  data.tmp <- seedlingData[seedlingData$site==site.tmp,]
  glmer.tmp <-lme4::glmer(cbind(fruitplNumber,seedlingNumber-fruitplNumber)~1 + (1|year) + (1|transect),data=data.tmp,family='binomial')
  glmer.list[[i]] <- glmer.tmp
}

f <- function(x,var=TRUE){
  if(var==FALSE){
    est  <-as.numeric(boot::inv.logit(fixef(x)[1]))
    y = est
  } else if(var==TRUE){
    yr<-VarCorr(x)[[1]][1,1]
    y = yr
   if(length(VarCorr(x))>1){ tr<-VarCorr(x)[[2]][1,1]
   y=c(yr,tr) }
    
  }
  return(y)
}
var.comps<-lapply(glmer.list,f,var=TRUE)
var.comps.mat <- do.call(rbind,var.comps)
alpha_y<-lapply(glmer.list,f,var=FALSE)
alpha.mat <- do.call(rbind,alpha_y)

# compare variance componenets for glmer fits

plot(var.comps.mat.1[,1],var.comps.mat[,1]); abline(a=0,b=1)

mu0 <- apply(boot::inv.logit(MCMCchains(mcmcSamples,params="mu0")),2,quantile,c(.025,.5,.975))
sigma_year <- apply(MCMCchains(mcmcSamples,params="sigma0"),2,quantile,c(.025,.5,.975))
sigma_transect <- apply(MCMCchains(mcmcSamples,params="sigma_transect"),2,quantile,c(.025,.5,.975))

plot(alpha.mat[,1],mu0[2,],xlim=c(0,1),ylim=c(0,1));abline(a=0,b=1)
segments(alpha.mat[,1],y0=mu0[1,],y1=mu0[3,]);abline(a=0,b=1)

plot(var.comps.mat[,1],sigma_year[2,],xlim=c(0,5),ylim=c(0,5));abline(a=0,b=1)
segments(var.comps.mat[,1],y0=sigma_year[1,],y1=sigma_year[3,]);abline(a=0,b=1)

plot(var.comps.mat[,2],sigma_transect[2,],xlim=c(0,2),ylim=c(0,2));abline(a=0,b=1)
segments(var.comps.mat[,2],y0=sigma_transect[1,],y1=sigma_transect[3,]);abline(a=0,b=1)


# - + calculate the mean ----

for(i in 1:20){
  # site.tmp=siteNames[i]
  site.tmp=i
  data.tmp <- seedlingData[seedlingData$site==site.tmp,]
  tmp = data.tmp %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(mu = sum(fruitplNumber)/sum(seedlingNumber),
                     n = sum(seedlingNumber,na.rm=TRUE))
  mean.list[[i]] <- tmp$mu
  names(mean.list[[i]]) <- (tmp$year)
  sampleSize.list[[i]] <- tmp$n
  names(sampleSize.list[[i]]) <- (tmp$year)
}

# - Extract estimates and build a matrix ----

#years = 2006:2020
years= 1:15

mu.freq = matrix(NA,nrow = 20,ncol = 15)
mu.freq.ranef = matrix(NA,nrow = 20,ncol = 15)
mu.mean = sample.mat = matrix(NA,nrow=20,ncol=15)
for(i in 1:20){
  year.tmp<-names(glm.list[[i]]$coefficients)
  year.tmp<-str_extract(year.tmp, "[[:digit:]]+")
  year.tmp.glmer<-rownames(ranef(glmer.list[[i]])$year)
  year.tmp.mean<-names((mean.list[[i]]))
  
  for(j in 1:15){
    if((years[j] %in% as.numeric(year.tmp))){
      
      mu.freq[i,j] = boot::inv.logit(glm.list[[i]]$coefficients[years[j] == as.numeric(year.tmp)])
    }
    if((years[j] %in% as.numeric(year.tmp.glmer))){
      
      mu.freq.ranef[i,j] = as.numeric(boot::inv.logit(fixef(glmer.list[[i]]) + ranef(glmer.list[[i]])$year$`(Intercept)`[years[j]==as.numeric(year.tmp.glmer)]))
    }
    
    if((years[j] %in% as.numeric(year.tmp.mean))){
      
      tmp.index = years[j] == as.numeric(year.tmp.mean)
      
      mu.mean[i,j] = mean.list[[i]][tmp.index]
      sample.mat[i,j] = sampleSize.list[[i]][tmp.index]
    }
  }
}

# - Construct a matrix for the posteriors ----

muPooled <- MCMCchains(mcmcSamples,params="mu")
sigmaPooled <- apply(muPooled,2,boot::inv.logit)
tmpDataDirectory = "outputs/001_prepareDataForModels/"
refDf <- readRDS(paste0(tmpDataDirectory,"seedlingRef.RDS"))

refDf<-refDf %>%
  dplyr::select(site,year,siteYearIndex) %>%
  unique

n.iter=dim(sigmaPooled)[1]
index.site<-unique(refDf$site)
dat.list <- list()
years = 2006:2020
#years=1:15
for(i in 1:20){
  index.tmp <- refDf[refDf$site==index.site[i],]
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  for(j in 1:15){
    if((years[j] %in% index.tmp$year)){
      index.tmp.tmp = as.numeric(index.tmp$siteYearIndex[index.tmp$year==years[j]])
      mat[,j] = sigmaPooled[,index.tmp.tmp]
    }
  }
  dat.list[[i]] <- mat
  
}

mu.bayes <- do.call(rbind,lapply(dat.list,apply,2,mean))

# - Compare estimates ----


# fixed effects vs. mixed effects
par(mfrow=c(1,1))
plot(mu.freq,mu.freq.ranef,xlim=c(0,1),ylim=c(0,1))

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i  in 1:20){
  plot(mu.freq[i,],mu.freq.ranef[i,],xlim=c(0,1),ylim=c(0,1),main=siteNames[i])
  abline(a=0,b=1)
}

# fixed effects vs. bayes
par(mfrow=c(1,1))
plot(mu.freq,mu.bayes,xlim=c(0,1),ylim=c(0,1))

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i  in 1:20){
  plot(mu.freq[i,],mu.bayes[i,],xlim=c(0,1),ylim=c(0,1),main=siteNames[i])
  abline(a=0,b=1)
}

# mixed effects vs. bayes
par(mfrow=c(1,1))
plot(mu.freq.ranef,mu.bayes,xlim=c(0,1),ylim=c(0,1))

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i  in 1:20){
  plot(mu.freq.ranef[i,],mu.bayes[i,],xlim=c(0,1),ylim=c(0,1),main=siteNames[i])
  abline(a=0,b=1)
}

par(mfrow=c(1,1))
plot(sample.mat,mu.bayes-mu.freq)

ifelse(abs(mu.bayes-mu.freq)>.05,1,0) %>% View
sample.mat
plot(sample.mat,abs(mu.bayes-mu.freq),xlim=c(0,100))

par(mfrow=c(1,1))
plot(sample.mat,abs(mu.bayes-mu.freq.ranef))
plot(sample.mat,abs(mu.bayes-mu.freq.ranef),xlim=c(0,50))
