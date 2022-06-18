####
####
# Script to prepare data from lab viability trials
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

viabilityData <- read.csv(paste0(dataDirectory,"viabilityData.csv"),header=TRUE)

# - Prep data for JAGS ----

viabilityData$bag <-as.integer(as.numeric(viabilityData$bagNo))

viabilityData = viabilityData %>%
  dplyr::mutate(year = as.factor(round)) %>%
  dplyr::select(site, year, age, germStart, germCount, viabStart, viabStain, bagNo) %>%
  dplyr::rename(siteViab = site,
                yearViab = year,
                ageViab = age,
                bag = bagNo) %>%
  dplyr::mutate(bag = as.factor(bag)) %>%
  dplyr::group_by(siteViab,yearViab, ageViab, bag) %>%
  # sum observations in each bag
  dplyr::summarise(germStart = sum(germStart,na.rm=TRUE),
                   germCount = sum(germCount,na.rm=TRUE),
                   viabStart = sum(viabStart,na.rm=TRUE),
                   viabStain = sum(viabStain,na.rm=TRUE)) 

# indexing to avoid issue of dealing with ragged arrays in JAGS
experimentalAge = c(1,1,1,2,2,3)
experimentalRound = c(1,2,3,1,2,1)
experimentalIndex = 1:6

df.index=data.frame(ageViab=experimentalAge,yearViab=experimentalRound,indexViab=experimentalIndex)

# viability data frame with event history data frame
viabilityData <- viabilityData %>%
  dplyr::mutate(yearViab=as.numeric(yearViab)) %>%
  dplyr::left_join(df.index,by=c("ageViab","yearViab")) %>%
  dplyr::mutate(yearViab=as.factor(yearViab))

# - +drop trials where viability trials start with 0 seeds ----

germinationTrialData <- viabilityData %>%
  dplyr::select(siteViab,yearViab,ageViab,bag,germStart,germCount,indexViab) 

viabilityTrialData <- viabilityData %>%
  dplyr::select(siteViab,yearViab,ageViab,bag,viabStart,viabStain,indexViab) %>%
  dplyr::filter(viabStart>0) %>%
  dplyr::rename(siteViab_v = siteViab,
                yearViab_v = yearViab, 
                ageViab_v = ageViab,
                bag_v = bag,
                indexViab_v = indexViab)

# - +Data to list ----
# JAGS takes a list of data

data <- tidybayes::compose_data(viabilityTrialData,germinationTrialData)
data$n <- dim(germinationTrialData)[1]
data$n_v <- dim(viabilityTrialData)[1]

# reassign variables from numeric to integer, or remove 1D notation
# numeric is 8-byte; integer is 4-byte
data$siteViab_v = as.integer(data$siteViab_v)
data$siteViab = as.integer(data$siteViab)
data$viabStart = as.integer(data$viabStart)
data$viabStain = as.integer(data$viabStain)
data$indexViab_v = as.integer(data$indexViab_v)
data$germStart = as.integer(data$germStart)
data$germCount = as.integer(data$germCount)
data$indexViab = as.integer(data$indexViab)

# - +Save data for model fitting ----
saveRDS(data,file=paste0(tmpDataDirectory,"viabilityData.RDS"))