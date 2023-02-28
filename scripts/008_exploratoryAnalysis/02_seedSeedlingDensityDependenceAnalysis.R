#################################################################################
# Script to conduct exploratory analysis of density-dependence in seed to seedling transition
# we conducted the analysis in this script after the first round of peer review
# a reviewer suggested we consider whether there was density-dependence in seedling establishment
# we followed methods from:
# Detto, M., M. D. Visser, S. J. Wright, and S. W. Pacala. 2019. Bias in the detection of negative density dependence in plant communities. Ecology Letters 22:1923–1939.
#################################################################################

# - Environment ----
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)

# - Source functions for analysis ----

# - Libraries ----
library(tidybayes)
library(tidyverse)
library(stringr)
library(lme4)

# - Read in data ----
dataDirectory = "~/Dropbox/clarkia-demography-projects/data/"
seedlingFruitingPlantCountsPermanentPlots <- read.csv(paste0(dataDirectory,"seedlingFruitingPlantCountsPermanentPlots.csv"),header=TRUE)

# - Prep data for analysis ----
# convert year to character

seedlingFruitingPlantCountsPermanentPlots$year <- as.character(seedlingFruitingPlantCountsPermanentPlots$year)

# if more fruiting plants than seedlings, set seedling number equal to fruiting plant number
seedlingFruitingPlantCountsPermanentPlots <- seedlingFruitingPlantCountsPermanentPlots %>%
  # more fruiting plants than seedlings
  # recode s.t. number of seedlings equals number of fruiting plants
  dplyr::mutate(seedlingNumber=ifelse(fruitplNumber>seedlingNumber,fruitplNumber,seedlingNumber))

# create index for years without observations (of seedlings, the denominator in fruitingplants/seedlings) across the entire site
tmp <- seedlingFruitingPlantCountsPermanentPlots %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n_tot = sum(seedlingNumber,na.rm=TRUE)) %>%
  dplyr::filter(n_tot>0) %>%
  dplyr::select(site,year)

# add index
tmp <- cbind(tmp,siteYearIndex = (1:(dim(tmp)[1])))

# create index for years without observations (of seedlings, the denominator in fruitingplants/seedlings) across the entire site
tmp2 <- seedlingFruitingPlantCountsPermanentPlots %>%
  dplyr::group_by(site,transect) %>%
  dplyr::summarise(n_tot = sum(seedlingNumber,na.rm=TRUE)) %>%
  dplyr::filter(n_tot>0) %>%
  dplyr::select(site,transect)

# add index
tmp2 <- cbind(tmp2,siteTransectIndex = (1:(dim(tmp2)[1])))

# filter out observations with 0 seedlings in individual plots
# or has missing data
seedlingFruitingPlantCountsPermanentPlots <- seedlingFruitingPlantCountsPermanentPlots %>%
  dplyr::filter(seedlingNumber>0)  %>%
  # NA for seedlings
  # filter these out; missing counts
  dplyr::filter(!is.na(seedlingNumber)) %>%
  # NA for fruiting plants
  #  filter these out; missing response
  dplyr::filter(!is.na(fruitplNumber))

# join the reference data frame tmp to the observations

seedlingFruitingPlantCountsPermanentPlots<- seedlingFruitingPlantCountsPermanentPlots %>%
  dplyr::left_join(tmp,by=c("site","year")) %>%
  dplyr::left_join(tmp2,by=c('site','transect'))


# - +Additional data wrangling ----

seedlingFruitingPlantCountsPermanentPlots$siteYearIndex <- as.numeric(seedlingFruitingPlantCountsPermanentPlots$siteYearIndex)
siteNames = unique(seedlingFruitingPlantCountsPermanentPlots$site)

# - Find optimal parameter for c ----

# make a grid of c from 0 to 1
c=seq(.01,1,by=.01)

# empty lists to hold model fitting output
model.list <- list()

# for each c, fit a logistic regression with covariate for seedling number^c for all populatoins
# for each c, then sum across all populations 
for(k in 1:length(c)){
  m1=numeric(0)
  for(j in 1:20){
    # filter data so that it's for the current site and only for 2006-2020 (drop 2021)
    df<-seedlingFruitingPlantCountsPermanentPlots %>%
      dplyr::filter(site==siteNames[j]) %>%
      dplyr::filter(year%in%as.character(2006:2020))# %>%

    # fit a logistic regression with density-dependence
    # fruiting plant numbers as successes, and seedling numbers as trials
    # fit model such that seedling number is raised to power of c, 
    # where c is a parameter that controls how nonlinear density-dependence is
    # model also includes a random effect for year
    m1[j]<-logLik(glmer(cbind(fruitplNumber,seedlingNumber-fruitplNumber)~1+I(seedlingNumber^c[k])+(1|year),
                        family=binomial,data=df))
  }
  model.list[[k]] <- sum(m1)
}

# - Plot figure showing c vs. log likelihood  ----

pdf(file=paste0("products/figures/densityDependence-seedlingFruitingPlantLog.pdf"),
    height=6,width=6,pointsize=12)

par(mar = c(6.1, 4.6, 4.1, 4.1), # change the margins
    lwd = 1, # increase the line thickness
    cex.axis = 1.2,
    xpd=F# increase default axis label size
)
plot(c,unlist(model.list),type='l',ylab="Log likelihood",
     cex.lab = 1.25, cex.axis = 1.25)
abline(v=c[which(unlist(model.list)==max(unlist(model.list)))],col='gray',lwd=2)
box()

dev.off()

# - Fit models with and without density-dependence ----

# use c.max for the value of c in the model with density-dependence
c.max=c[which(max(unlist(model.list))==(unlist(model.list)))]

# empty lists to hold model output
# m0 is the family of models without density-dependence
# m1 is the family of models with density-dependence
m0.list <- m1.list <- list()

# refit models separately for each population i
for(i in 1:20){
  # filter data to current population and to include years 2006-2020
  df<-seedlingFruitingPlantCountsPermanentPlots %>%
    dplyr::filter(site==siteNames[i]) %>%
    dplyr::filter(year%in%as.character(2006:2020))
  
  # fit a logistic regression without density-dependence
  # fruiting plant numbers as successes, and seedling numbers as trials
  # model also includes a random effect for year
  m0.list[[i]]<-glmer(cbind(fruitplNumber,seedlingNumber-fruitplNumber)~1+(1|year),
                          family=binomial,data=df)
  
  # fit a logistic regression with density-dependence
  # fruiting plant numbers as successes, and seedling numbers as trials
  # fit model such that seedling number is raised to power of c, 
  # where c is a parameter that controls how nonlinear density-dependence is
  # model also includes a random effect for year
  m1.list[[i]]<-glmer(cbind(fruitplNumber,seedlingNumber-fruitplNumber)~1+I(seedlingNumber^c.max)+(1|year),
                                family=binomial,data=df)
}


# - Model comparison via likelihood ratio test ----

# empty list to hold model comparison output
lrt.list <- list()
for(i in 1:20){
  model.out<-anova(m0.list[[i]],m1.list[[i]])
  lrt.list[[i]]<-c(model.out$logLik,model.out$Chisq,model.out$`Pr(>Chisq)`)
}

# write out as matrix and select relevant columns
tmp.mat<-do.call(rbind,lrt.list)
lrt.mat<-tmp.mat[,c(1:2,4,6)]

# check significance with bonferroni correction
cbind(lrt.mat,lrt.mat[,4]<(.05/20))
siteNames[lrt.mat[,4]>(.05/20)]

# get site geographic position and use to reorder data frame
siteAbiotic <- read.csv("~/Dropbox/clarkia-demography-projects/data/siteAbioticData.csv",header=TRUE)

lrt.df=data.frame(site=siteNames,lrt.mat)
lrt.df<-lrt.df[order(siteAbiotic$easting),]

# set column names and number of digits to print
names(lrt.df) = c("pop","log-lik0","log-lik1",'chi-sq','p-val')
lrt.df[,2:4]<-signif(lrt.df[,2:4],5)
lrt.df[,5]<-signif(lrt.df[,5],6)

print(xtable::xtable(lrt.df,digits=c(0,0,2,2,2,4),row.names=FALSE),
            file="~/Dropbox/clarkia-bet-hedging/manuscript/02_supplement/model-checks/seedSeedlingDensityDependenceLRT.txt",
            )

# - Plots to examine density-dependence ----

par(mfrow=c(4,5),mar=c(2,2,1,1))
for(j in 1:20){
  df.tmp<-seedlingFruitingPlantCountsPermanentPlots %>%
    dplyr::filter(site==siteNames[j]) %>%
    dplyr::mutate(z = seedlingNumber^c.max)%>%
    dplyr::filter(year%in%as.character(2006:2020))
  obsYear = unique(df.tmp$year)
  plot(NA,xlim=c(1,max(df.tmp$seedlingNumber)),ylim=c(0,1))
  legend('topright',siteNames[j],bty='n')
  for(i in 1:length(obsYear)){
    yrs=as.character(2006:2020)[as.character(2006:2020) %in% obsYear]
    df.tmp2<-seedlingFruitingPlantCountsPermanentPlots %>%
      dplyr::filter(site==siteNames[j]) %>%
      dplyr::filter(year==yrs[i])
    
    if(length(1:max(df.tmp2$seedlingNumber))>1){
      lines(1:max(df.tmp2$seedlingNumber),
            boot::inv.logit(predict(m1.list[[j]],data.frame(seedlingNumber=seq(1,max(df.tmp2$seedlingNumber),by=1),year=(obsYear)[i]))),
            col=ifelse(lrt.mat[j,4]>(.05/20),'gray','black'))
    } else {
      points(1:max(df.tmp2$seedlingNumber),
             boot::inv.logit(predict(m1.list[[j]],data.frame(seedlingNumber=seq(1,max(df.tmp2$seedlingNumber),by=1),year=(obsYear)[i]))),
             pch=16,
             col=ifelse(lrt.mat[j,4]>(.05/20),'gray','black'))
    }
  }
}

for(j in 1:20){
  par(mfrow=c(3,5),mar=c(2,2,1,1))
  
  df.tmp<-seedlingFruitingPlantCountsPermanentPlots %>%
    dplyr::filter(site==siteNames[j]) %>%
    dplyr::filter(year%in%as.character(2006:2020))
  obsYear = unique(df.tmp$year)
  
  for(i in 1:length(obsYear)){
    
    yrs=as.character(2006:2020)[as.character(2006:2020) %in% obsYear]
    df.tmp2<-seedlingFruitingPlantCountsPermanentPlots %>%
      dplyr::filter(site==siteNames[j]) %>%
      dplyr::filter(year==yrs[i])
    plot(NA,xlim=c(1,max(df.tmp2$seedlingNumber)),ylim=c(0,1))
    legend('topright',siteNames[j],bty='n')
    
    points(df.tmp2$seedlingNumber,
           df.tmp2$fruitplNumber/df.tmp2$seedlingNumber,
           pch=16,cex=1.5,col='gray90')
    
    if(length(1:max(df.tmp2$seedlingNumber))>1){
      lines(1:max(df.tmp2$seedlingNumber),
            boot::inv.logit(predict(m1.list[[j]],data.frame(seedlingNumber=seq(1,max(df.tmp2$seedlingNumber),by=1),year=(obsYear)[i]))))
      lines(0:max(df.tmp2$seedlingNumber),boot::inv.logit(predict(m0.list[[j]],data.frame(seedlingNumber=seq(0,max(df.tmp2$seedlingNumber),by=1),year=(obsYear)[i]))),
            col='red',lty='dotted')
      abline(v=mean(df.tmp2$seedlingNumber),col='red',lty='dotted')
    } else {
      points(1:max(df.tmp2$seedlingNumber),
             boot::inv.logit(predict(m1.list[[j]],data.frame(seedlingNumber=seq(1,max(df.tmp2$seedlingNumber),by=1),year=(obsYear)[i]))))
    }
  }
}

