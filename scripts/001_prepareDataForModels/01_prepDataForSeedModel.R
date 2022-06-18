####
####
# Script to prepare field observations from seed bag burial experiment
# for fitting in JAGS
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

dataDirectory = "data/"
tmpDataDirectory = "outputs/001_prepareDataForModels/"

# - Libraries ----
library(tidybayes)
library(tidyverse)
library(stringr)

# - Read in data ----

# - +Seed bag burial experiment ----
seedBagsData <- read.csv(paste0(dataDirectory,"seedBagsData.csv"),header=TRUE)

# - +Seed production and seedlings in plots ----
seedlingFruitingPlantCountsPermanentPlots <- read.csv(paste0(dataDirectory,"seedlingFruitingPlantCountsPermanentPlots.csv"),header=TRUE)
countFruitsPerPlantFromPermanentPlots <- read.csv(paste0(dataDirectory,"countFruitsPerPlantFromPermanentPlots.csv"),header=TRUE)
countSeedPerFruit <- read.csv(paste0(dataDirectory,"countSeedPerFruit.csv"),header=TRUE) 

# - +Seed production in all plots ----
countFruitsPerPlantAllPlants <- read.csv(paste0(dataDirectory,"countFruitsPerPlantAllPlants.csv"),header=TRUE)

# - Select seedling observations ----
## remove rows with missing data in January for germination data
## this keeps only the data that have observations of seedlings and totals in January

seedlingData <- seedBagsData %>%
  dplyr::filter(!is.na(totalJan)&!is.na(seedlingJan))

# - +Rename datasets ----
seedBagsDataIntactSeeds <- seedBagsData 
seedBagsDataSeedlings <- seedlingData

# - Prep data for JAGS ----

# - +Seed bags: intact seeds ----

# Create a variable for the number of seeds at the start of the trial
seedBagsDataIntactSeeds$seedStart<-as.double(100)

# rename variables and organize
seedBagsDataIntactSeeds = seedBagsDataIntactSeeds %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag))   %>%
  dplyr::mutate(year = as.factor(yearStart),
                age = as.factor(age)) %>%
  dplyr::select(siteBag,site,year,age,totalJan,seedStart,intactOct) %>%
  dplyr::rename(siteBags = site,
                yearBags = year,
                ageBags = age)

# create a reference dataframe
refIndex <- data.frame(name=rep(c("totalJan","intactOct"),3),
                       ageBags=as.factor(c(1,1,2,2,3,3)),
                       months=c(3,12,15,24,27,36)) %>%
  dplyr::mutate(months=months/36)

# pivot to storing data in long format indexed by response type
seedBagsDataIntactSeeds <- seedBagsDataIntactSeeds %>%
  tidyr::pivot_longer(cols=c(totalJan,intactOct)) %>%
  dplyr::left_join(refIndex,by=c("name","ageBags")) %>%
  dplyr::rename(response = name, y = value, siteSurvival = siteBags) %>%
  dplyr::mutate(yearSurvival = yearBags)

# create a reference data frame for event histories
referenceDataFrameEvents <- data.frame(yearSurvival = as.factor(c(rep(2005,6),rep(2006,4),rep(2007,2))),
                                       ageBags=as.factor(c(1,1,2,2,3,3,1,1,2,2,1,1)),
                                       compIndex=c(1:12),
                                       response=rep(c("totalJan","intactOct"),6)) 

# join intact seed data frame with event history data frame
seedBagsDataIntactSeeds <- seedBagsDataIntactSeeds %>% 
  dplyr::left_join(referenceDataFrameEvents,by=c('yearSurvival','ageBags',"response"))

# - +Seed bags: seedlings ----

# rename variables and organize
seedBagsDataSeedlings <- seedBagsDataSeedlings %>%
  tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
  tidyr::unite(col='siteBag', c(site,bagNo), sep="-", remove=FALSE) %>%
  dplyr::mutate(siteBag = as.factor(siteBag))   %>%
  dplyr::mutate(year = as.factor(yearStart),
                age = as.factor(age)) %>%
  dplyr::select(siteBag,site,year,age,totalJan,seedlingJan) %>%
  dplyr::rename(siteGermination = site,
                yearGermination = year,
                ageBags = age) %>%
  dplyr::mutate(gIndex=ageBags)

# create reference data frame to resolve indexing
# index to avoid issue of dealing with ragged arrays in JAGS
germinationAge = c(1,2,3,1,2,1)
germinationYear = c(2005,2005,2005,2006,2006,2007)
germinationIndex = 1:6
reducedGerminationIndex = c(1,2,3,4,3,5)

referenceDataFrame=data.frame(ageBags=as.factor(germinationAge),
                              yearGermination=as.factor(germinationYear),
                              germinationIndex=germinationIndex,
                              reducedGerminationIndex = reducedGerminationIndex)

# join observation and reference dataframes
seedBagsDataSeedlings <- seedBagsDataSeedlings %>%
  dplyr::left_join(referenceDataFrame,by=c("ageBags","yearGermination")) %>%
  dplyr::mutate(yearGermination=as.factor(yearGermination),ageBags=as.factor(ageBags))

# - +Seed production and seedlings in plots ----

# - ++Seedling numbers in permanent plots ----
# data from 2007, 2008
# calculate total number of seedlings per plot each year

# rename variables and organize
censusSeedlings <- seedlingFruitingPlantCountsPermanentPlots %>%
  dplyr::mutate(site=as.character(site),transect=as.character(transect),position=as.character(position)) %>%
  dplyr::filter(year%in%c(2007,2008)) %>%
  dplyr::mutate(yearRef=year-1) %>%
  dplyr::rename(t1=year) %>%
  dplyr::select(-fruitplNumber)

# - ++Fruit totals in permanent plots ----
# data from 2006, 2007
# calculate total number of fruits per plot each year

countFruitsPerPlotTransects <- countFruitsPerPlantFromPermanentPlots %>%
  dplyr::filter(year%in%c(2006,2007)) %>%
  dplyr::group_by(site,year,transect,position) %>%
  dplyr::summarise(totalFruitsPerPlot = sum(countFruitsPerPlant))

# - ++Average fruit per plant from all surveyed plants ----
# data from 2006, 2007
# calculate average number of fruits per plant each year
averageCountFruitsPerPlantAllPlants <- countFruitsPerPlantAllPlants %>%
  dplyr::filter(year%in%c(2006,2007)) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(averageFruitsPerPlant = mean(countFruitNumberPerPlant))

# - ++Seeds per fruit ----
# data from 2006, 2007
# calculate average number of seeds per fruit

countSeedPerTotalFruitEquivalent <- countSeedPerFruit %>%
  dplyr::filter(year %in% c(2006,2007)) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(meanSeedPerTotalFruitEquivalent = mean(sdno))

# - ++Seed rain ----

# calculate average seed rain per plot using counts of fruits per plant from plots
# round to whole number for use in binomial 
seedRainPlot <- countFruitsPerPlotTransects %>%
  dplyr::left_join(countSeedPerTotalFruitEquivalent,by=c("site","year")) %>%
  dplyr::mutate(fecundityEstimate = round(totalFruitsPerPlot*meanSeedPerTotalFruitEquivalent)) %>%
  dplyr::mutate(t=year) %>% 
  dplyr::rename(yearRef=year) %>%
  dplyr::mutate(site=as.character(site),transect=as.character(transect),position=as.character(position))

# calculate number of fruiting plants per plot  then 
# calculate the seed rain per plot

n_2006 <- seedlingFruitingPlantCountsPermanentPlots %>%
  dplyr::filter(year==2006) %>%
  dplyr::select(site,transect,position,year,fruitplNumber) %>%
  dplyr::left_join(averageCountFruitsPerPlantAllPlants,by=c('site','year')) %>%
  dplyr::left_join(countSeedPerTotalFruitEquivalent,by=c('site','year')) %>%
  dplyr::mutate(n_2006 = round(fruitplNumber*averageFruitsPerPlant*meanSeedPerTotalFruitEquivalent))

n_2006 <- n_2006 %>%
  dplyr::select(site,transect,position,n_2006)

n_2007 <- seedlingFruitingPlantCountsPermanentPlots %>%
  dplyr::filter(year==2007) %>%
  dplyr::select(site,transect,position,year,fruitplNumber) %>%
  dplyr::left_join(averageCountFruitsPerPlantAllPlants,by=c('site','year')) %>%
  dplyr::left_join(countSeedPerTotalFruitEquivalent,by=c('site','year')) %>%
  dplyr::mutate(n_2007 = round(fruitplNumber*averageFruitsPerPlant*meanSeedPerTotalFruitEquivalent))

n_2007 <- n_2007 %>%
  dplyr::select(site,transect,position,n_2007)

# - ++Seedling numbers in permanent plots ----
# data from 2007, 2008
# calculate total number of seedlings per plot each year

# rename variables and organize
censusSeedlings2 <- seedlingFruitingPlantCountsPermanentPlots %>%
  dplyr::mutate(site=as.character(site),transect=as.character(transect),position=as.character(position)) %>%
  dplyr::filter(year%in%c(2008)) %>%
  dplyr::select(-year,-fruitplNumber)

# calculate average seed rain per plot using counts of fruiting plants per plot
seedRainSeedlings <- n_2006 %>%
  dplyr::left_join(n_2007,by=c("site","transect","position")) %>%
  dplyr::mutate(position=as.character(position)) %>%
  dplyr::left_join(censusSeedlings2) %>%
  dplyr::rename(n_t2 = n_2006, n_t1 = n_2007)

# aggregate seed rain/seedlings to the transect level 
# for this part of the analysis
seedRainSeedlingsTransect = seedRainSeedlings %>%
  dplyr::group_by(site,transect) %>%
  dplyr::summarise(n_t2 = sum(n_t2), 
                   n_t1 = sum(n_t1), 
                   seedlingNumber = sum(seedlingNumber))

# Rename variables and prepare data for JAGS
seedRainSeedlingsTransect = seedRainSeedlingsTransect %>% 
  dplyr::rename(sitePlot = site) %>%
  dplyr::rename(plotSeedlings = seedlingNumber) 

# create reference data frame for indexing
# indexing to avoid ragged array
df.index=data.frame(ageBags=as.factor(germinationAge),yearGermination=as.factor(germinationYear),germinationIndex=germinationIndex,reducedGerminationIndex=as.factor(reducedGerminationIndex))
df.index=df.index %>%
  dplyr::filter(yearGermination==2006|yearGermination==2007) %>%
  droplevels() %>%
  dplyr::rename(yearPlot=yearGermination) %>%
  dplyr::rename(fecIndex = germinationIndex) %>%
  dplyr::select(yearPlot,fecIndex)

# join reference data frame for indexing
referenceDataFrameEvents.index=referenceDataFrameEvents %>%
  dplyr::filter(yearSurvival==2006|yearSurvival==2007) %>%
  droplevels() %>%
  dplyr::rename(yearPlot=yearSurvival) %>%
  dplyr::select(yearPlot,compIndex) %>%
  dplyr::filter(compIndex%in%c(7:9,11)) %>%
  dplyr::rename(fecCompIndex = compIndex)

seedRainSeedlingsTransect <- seedRainSeedlingsTransect[!(seedRainSeedlingsTransect$n_t2==0&seedRainSeedlingsTransect$n_t1==0&seedRainSeedlingsTransect$plotSeedlings>0),]

# - +Data to list ----
# JAGS takes a list of data
data <- tidybayes::compose_data(seedBagsDataIntactSeeds,seedBagsDataSeedlings,seedRainSeedlingsTransect)

# variable with length of each dataset
data$n1 <- dim(seedBagsDataIntactSeeds)[1]
data$n2 <- dim(seedBagsDataSeedlings)[1]
data$n3 <- dim(seedRainSeedlingsTransect)[1]

# reassign variables from numeric to integer, or remove 1D notation
# numeric is 8-byte; integer is 4-byte
data$siteSurvival = as.integer(data$siteSurvival)
data$seedStart = as.integer(data$seedStart)
data$y = as.integer(data$y)
data$yearSurvival = as.integer(data$yearSurvival)
data$compIndex = as.integer(data$compIndex)
data$siteGermination = as.integer(data$siteGermination)
data$yearGermination = as.integer(data$yearGermination)
data$totalJan = as.integer(data$totalJan)
data$seedlingJan = as.integer(data$seedlingJan)
data$germinationIndex = as.integer(data$germinationIndex)
data$n_t2 = as.integer(data$n_t2)
data$n_t1 = as.integer(data$n_t1)
data$plotSeedlings = as.integer(data$plotSeedlings)

# - +Save data for model fitting ----
saveRDS(data,file=paste0(tmpDataDirectory,"seedData.RDS"))

saveRDS(seedRainSeedlingsTransect,file=paste0(tmpDataDirectory,"seedBagsMatrixForInitialization.RDS"))