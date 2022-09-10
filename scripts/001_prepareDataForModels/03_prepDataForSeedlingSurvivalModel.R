####
####
# Script to prepare observations of seedlings and fruiting plants
# from permanent plots for fitting in JAGS
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

seedlingFruitingPlantCountsPermanentPlots <- read.csv(paste0(dataDirectory,"seedlingFruitingPlantCountsPermanentPlots.csv"),header=TRUE)

# - Prep data for JAGS ----
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


# - +Data to list ----
# JAGS takes data as a list;
# tidybayes transforms 'tidy' data to list form for JAGS
# for more information see: http://mjskay.github.io/tidybayes/articles/tidybayes.html

seedlingFruitingPlantCountsPermanentPlots$siteYearIndex <- as.numeric(seedlingFruitingPlantCountsPermanentPlots$siteYearIndex)
data <- tidybayes::compose_data(seedlingFruitingPlantCountsPermanentPlots)

# convert all vectors to lists for more efficient storage vs. numeric
data$site = as.integer(data$site)
data$transect = as.integer(data$transect)
data$position = as.integer(data$position)
data$year = as.integer(data$year)
data$seedlingNumber = as.integer(data$seedlingNumber)
data$fruitplNumber = as.integer(data$fruitplNumber)
data$siteYearIndex = as.integer(data$siteYearIndex)
data$siteTransectIndex = as.integer(data$siteTransectIndex)

# finish adding the ragged reference array to the list
tmp$siteYearIndex=as.integer(tmp$siteYearIndex)
tmp = tidybayes::compose_data(tmp)

data$site_observed = as.integer(tmp$site)
data$year_observed = as.integer(tmp$year)
data$siteYearIndex_observed = as.integer(tmp$siteYearIndex)
data$n_siteYearIndex = tmp$n

# add array for transects
tmp2 = tidybayes::compose_data(tmp2)

data$site_transect_observed = as.integer(tmp2$site)
data$transect_observed = as.integer(tmp2$transect)
data$siteTransectIndex_observed = as.integer(tmp2$siteTransectIndex)
data$n_siteTransectIndex = tmp2$n

# check arrays
check.n = sample(1:data$n_siteYearIndex,1)
data$site_observed[check.n]
data$site[data$siteYearIndex==check.n]

check.n = sample(1:data$n_siteTransectIndex,1)
data$site_transect_observed[check.n]
data$site[data$siteTransectIndex==check.n]

# - +Save data for model fitting ----
saveRDS(data,file=paste0(tmpDataDirectory,"seedlingFruitingPlantCountsPermanentPlots.RDS"))

# - +Save reference data frame ----
saveRDS(seedlingFruitingPlantCountsPermanentPlots,file=paste0(tmpDataDirectory,"seedlingRef.RDS"))