
rm(list=(ls(all=TRUE))) # if using in source(script)

# - Libraries ----
library(tidyverse)
library(parallel)
library(stringr)
library(bbmle)
library(MCMCvis)
library(lme4)

# - +Read in posterior for seed bag trials & viabiility trials ----


# - +Read in site names ----
siteAbiotic <- read.csv("data/siteAbiotic.csv",header=TRUE)
siteNames <- siteAbiotic$site

# - Read in data to fit ML models ----

tmpDataDirectory = "outputs/001_prepareDataForModels/"

tmpDataDirectory <- paste0(tmpDataDirectory,list.files(tmpDataDirectory))
data <- readRDS(tmpDataDirectory[[grep("fruitsPerPlantAllPlants.RDS",tmpDataDirectory)]])

# - Prepare data for TFE ----

df = data.frame(site=data$site,year=data$year,y=data$y_tfe)
df$site = as.factor(df$site)
df$year = as.factor(df$year)

df_tot <- df %>%
  dplyr::group_by(site,year)

# split up data frame into lists per site
df_tot.list <- split( df_tot , f = df_tot$site )

model.list <- estimate.mat <- summary.list <-list()
matrix.aicc <- matrix(NA,nrow=20,ncol=2)
for(i in 1:20){
  site.tmp=i
  data.tmp <- df[df$site==site.tmp,]
  data.tmp$obs = as.factor(1:length(data.tmp$y))
  glmer.tmp <-lme4::glmer(y~1+(1|year),data=data.tmp,family='poisson')
  glmer.olre.tmp <-lme4::glmer(y~1+(1|year)+(1|obs),data=data.tmp,family='poisson')
  model.list[[i]] <- list(mod1 = glmer.tmp, mod2 = glmer.olre.tmp)
  year.est<-cbind(exp(fixef(glmer.tmp)[1]+ranef(glmer.tmp)$year[,1]),
        exp(fixef(glmer.olre.tmp)[1]+ranef(glmer.olre.tmp)$year[,1]))
  estimate.mat[[i]] <- year.est
  matrix.aicc[i,]<-c(MuMIn::AICc(glmer.tmp),MuMIn::AICc(glmer.olre.tmp))
  summary.list[[i]]<-data.tmp %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(mu=mean(y),var=var(y),n=n())

}

summary.mat<-do.call(rbind,summary.list)
plot(summary.mat$mu,summary.mat$var);abline(a=0,b=1)
plot(summary.mat$n,summary.mat$var)

# are the models with OLRE better?
matrix.aicc

do.call(max,estimate.mat)
plot(NA,xlim=c(0,30),ylim=c(0,30),xlab="Model: year",ylab="Model:year+OLRE")
matrix.lmer <- matrix(NA,nrow=20,ncol=2)
for(i in 1:20){
  tmp.mat<-estimate.mat[[i]]
  points(tmp.mat[,1],tmp.mat[,2])
  matrix.lmer[i,] <- c(var(tmp.mat[,1]),var(tmp.mat[,2]))
}
abline(a=0,b=1)


# - Calculate variance with White (2000) correction ----
source("~/Dropbox/chapter-3/analysis/scripts/001_geographyVariability/00_functions.R")

# - +function sourced above ----

tfe.var.est = c()
for(j in 1:20){
  df.tmp <- df_tot.list[[j]]
  d.tmp <- data.frame(df.tmp$y, df.tmp$year)
  
  optim.out=optim(par=c(1),variance.white,lower=0,upper=1000,
                  dat=d.tmp,control=list(maxit=10000),method="Brent")
  tfe.var.est[j]=optim.out$par
}

# the model with overdispersion underestimates the variance relative to 
# the method proposed in White 2000
plot(tfe.var.est,matrix.lmer[,1]);abline(a=0,b=1)
plot(tfe.var.est,matrix.lmer[,2]);abline(a=0,b=1)


# - Prepare data for total fruits ----

df = data.frame(site=data$site2,year=data$year2,y=data$y_tot)
df$site = as.factor(df$site)
df$year = as.factor(df$year)

df_tot <- df %>%
  dplyr::group_by(site,year)

# split up data frame into lists per site
df_tot.list <- split( df_tot , f = df_tot$site )

# - Calculate variance with White (2000) correction ----

# - +function sourced above ----

tot.var.est = c()
for(j in 1:20){
  df.tmp <- df_tot.list[[j]]
  d.tmp <- data.frame(df.tmp$y, df.tmp$year)
  
  optim.out=optim(par=c(1),variance.white,lower=0,upper=1000,
                  dat=d.tmp,control=list(maxit=10000),method="Brent")
  tot.var.est[j]=optim.out$par
}

tot.var.est

plot(siteAbiotic$easting,tot.var.est,ylim=c(0,150))



# - Set up model ----

glm.list <- glmer.list <- glmer.list.olre <- mean.list <- sampleSize.list <- list()

# - + fixed effects model ----



# - + mixed effects effects models ----

for(i in 1:20){
  # site.tmp=siteNames[i]
  site.tmp=i
  data.tmp <- seedsPerFruitData[seedsPerFruitData$site==site.tmp,]
  glmer.tmp <-lme4::glmer(sdno~1+ (1|year),data=data.tmp,family='poisson')
  glmer.list[[i]] <- glmer.tmp
}

for(i in 1:20){
  # site.tmp=siteNames[i]
  site.tmp=i
  data.tmp <- seedsPerFruitData[seedsPerFruitData$site==site.tmp,]
  data.tmp$observation <- as.factor(1:length(data.tmp$sdno))
  glmer.tmp <-lme4::glmer(sdno~1+ (1|year) + (1|observation),data=data.tmp,family='poisson')
  glmer.list.olre[[i]] <- glmer.tmp
}



# - + calculate the mean ----

for(i in 1:20){
  # site.tmp=siteNames[i]
  site.tmp=i
  data.tmp <- seedsPerFruitData[seedsPerFruitData$site==site.tmp,]
  tmp = data.tmp %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(mu = mean(sdno,na.rm=TRUE),
                     n = n())
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
mu.freq.ranef.olre = matrix(NA,nrow = 20,ncol = 15)
mu.mean = sample.mat = matrix(NA,nrow=20,ncol=15)
for(i in 1:20){
  year.tmp<-names(glm.list[[i]]$coefficients)
  year.tmp<-str_extract(year.tmp, "[[:digit:]]+")
  year.tmp.glmer<-rownames(ranef(glmer.list[[i]])$year)
  year.tmp.mean<-names((mean.list[[i]]))
  
  for(j in 1:15){
    if((years[j] %in% as.numeric(year.tmp))){
      
      mu.freq[i,j] = exp(glm.list[[i]]$coefficients[years[j] == as.numeric(year.tmp)])
    }
    if((years[j] %in% as.numeric(year.tmp.glmer))){
      
      mu.freq.ranef[i,j] = as.numeric(exp(fixef(glmer.list[[i]]) + ranef(glmer.list[[i]])$year$`(Intercept)`[years[j]==as.numeric(year.tmp.glmer)]))
      mu.freq.ranef.olre[i,j] = as.numeric(exp(fixef(glmer.list.olre[[i]]) + ranef(glmer.list.olre[[i]])$year$`(Intercept)`[years[j]==as.numeric(year.tmp.glmer)]))
    }
    
    if((years[j] %in% as.numeric(year.tmp.mean))){
      
     tmp.index = years[j] == as.numeric(year.tmp.mean)
      
      mu.mean[i,j] = mean.list[[i]][tmp.index]
      sample.mat[i,j] = sampleSize.list[[i]][tmp.index]
    }
  }
}

# - Compare mixed model with and without OLRE ----

plot(mu.freq.ranef,mu.freq.ranef.olre);abline(a=0,b=1)

# - Construct a matrix for the posteriors ----

muPooled <- MCMCchains(mcmcSamples,params="mu_seeds")
#sigmaPooled <- apply(muPooled,2,boot::inv.logit)
tmpDataDirectory = "outputs/001_prepareDataForModels/"
refDf <- readRDS(paste0(tmpDataDirectory,"referenceSeedsPerFruit.RDS"))

refDf<-refDf %>%
  dplyr::select(site,year,siteYearIndex_und) %>%
  unique

n.iter=dim(muPooled)[1]
index.site<-unique(refDf$site)
dat.list <- list()
years = 2006:2020
#years=1:15
for(i in 1:20){
  index.tmp <- refDf[refDf$site==index.site[i],]
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  for(j in 1:15){
    if((years[j] %in% index.tmp$year)){
      index.tmp.tmp = as.numeric(index.tmp$siteYearIndex_und[index.tmp$year==years[j]])
      mat[,j] = muPooled[,index.tmp.tmp]
    }
  }
  dat.list[[i]] <- mat
  
}

mu.bayes <- do.call(rbind,lapply(dat.list,apply,2,mean))

# - Compare estimates ----


# fixed effects vs. mixed effects
par(mfrow=c(1,1))
plot(mu.freq,mu.freq.ranef)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i  in 1:20){
  plot(mu.freq[i,],mu.freq.ranef[i,],main=siteNames[i])
  abline(a=0,b=1)
}

# fixed effects vs. bayes
par(mfrow=c(1,1))
plot(mu.freq,mu.bayes);abline(a=0,b=1)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i  in 1:20){
  plot(mu.freq[i,],mu.bayes[i,],main=siteNames[i])
  abline(a=0,b=1)
}

# mixed effects vs. bayes
par(mfrow=c(1,1))
plot(mu.freq.ranef,mu.bayes);abline(a=0,b=1)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for(i  in 1:20){
  plot(mu.freq.ranef[i,],mu.bayes[i,],main=siteNames[i])
  abline(a=0,b=1)
}

par(mfrow=c(1,1))
plot(sample.mat,mu.bayes-mu.freq)

#ifelse(abs(mu.bayes-mu.freq)>.05,1,0) %>% View
sample.mat
plot(sample.mat,abs(mu.bayes-mu.freq),xlim=c(0,100))

par(mfrow=c(1,1))
plot(sample.mat,abs(mu.bayes-mu.freq.ranef))
plot(sample.mat,abs(mu.bayes-mu.freq.ranef),xlim=c(0,50))

