####
####
# Script to prepare data on fruits per plant
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

countFruitsPerPlantAllPlants <- read.csv(paste0(dataDirectory,"countFruitsPerPlantAllPlants.csv"),header=TRUE)
countUndamagedDamagedFruitsPerPlantAllPlants <- read.csv(paste0(dataDirectory,"countUndamagedDamagedFruitsPerPlantAllPlants.csv"),header=TRUE) 

# - Prep data for JAGS ----

# - +Total fruit equivalents (2006-2012) ----

countFruitsPerPlantAllPlants <- countFruitsPerPlantAllPlants %>%
  dplyr::rename(y_tfe = countFruitNumberPerPlant) %>%
  dplyr::select(site,year,y_tfe) %>%
  dplyr::filter(!is.na(y_tfe))

countFruitsPerPlantAllPlants$year <- as.character(countFruitsPerPlantAllPlants$year)

# create index for years with observations across the entire site
tmp <- countFruitsPerPlantAllPlants %>%  
  dplyr::select(site,year) %>%
  unique

# add index
tmp_tfe <- cbind(tmp,siteYearIndex_tfe = (1:(dim(tmp)[1])))

# join the reference data frame tmp to the observations
countFruitsPerPlantAllPlants<- countFruitsPerPlantAllPlants %>%
  dplyr::left_join(tmp_tfe,by=c("site","year"))

# check reference df
check.n = sample(1:139,1)
as.integer(as.factor(countFruitsPerPlantAllPlants$site))[countFruitsPerPlantAllPlants$siteYearIndex_tfe==check.n]
as.integer(as.factor(tmp_tfe$site))[tmp_tfe$siteYearIndex_tfe==check.n]

# - +Undamaged/damaged fruits (2013-2020) ----

countUndamagedDamagedFruitsPerPlantAllPlants <- countUndamagedDamagedFruitsPerPlantAllPlants %>%
  dplyr::rename(y_und = countUndamagedFruitNumberPerPlant) %>%
  dplyr::rename(y_dam = countDamagedFruitNumberPerPlant) %>%
  dplyr::rename(site2 = site) %>%
  dplyr::rename(year2 = year) %>%
  dplyr::select(site2,year2,y_und,y_dam) %>%
  dplyr::filter(!is.na(y_und)) %>%
  # create variable for total number of fruits
  dplyr::mutate(y_tot = y_und + y_dam)

countUndamagedDamagedFruitsPerPlantAllPlants$year2 <- as.character(countUndamagedDamagedFruitsPerPlantAllPlants$year2)

# create index for years with observations across the entire site
tmp2 <- countUndamagedDamagedFruitsPerPlantAllPlants %>%
  dplyr::select(site2,year2) %>%
  unique

# add index
tmp_tot <- cbind(tmp2,siteYearIndex_tot = (1:(dim(tmp2)[1])))

# join the reference data frame tmp to the observations
countUndamagedDamagedFruitsPerPlantAllPlants<- countUndamagedDamagedFruitsPerPlantAllPlants %>%
  dplyr::left_join(tmp_tot,by=c("site2","year2"))

# - +Data to list ----
# JAGS takes a list of data
data <- tidybayes::compose_data(countFruitsPerPlantAllPlants,
                                countUndamagedDamagedFruitsPerPlantAllPlants)

data$n = dim(countFruitsPerPlantAllPlants)[1]
data$n2 = dim(countUndamagedDamagedFruitsPerPlantAllPlants)[1]

og.size = object.size(data)
# reassign variables from numeric to integer, or remove 1D notation
# numeric is 8-byte; integer is 4-byte
data$site = as.integer(data$site)
data$year = as.integer(data$year)
data$y_tfe = as.integer(data$y_tfe)
data$siteYearIndex_tfe = as.integer(data$siteYearIndex_tfe)

data$site2 = as.integer(data$site2)
data$year2 = as.integer(data$year2)
data$y_tot = as.integer(data$y_tot)
data$y_dam = as.integer(data$y_dam)
data$y_und = as.integer(data$y_und)
data$siteYearIndex_tot = as.integer(data$siteYearIndex_tot)

new.size = object.size(data)
# the new object is 30% smaller!
new.size/og.size

# finish adding the ragged reference array to the list
tmp_tfe = tidybayes::compose_data(tmp_tfe)

data$site_tfe_observed = as.integer(tmp_tfe$site)
data$year_tfe_observed = as.integer(tmp_tfe$year)
data$siteYearIndex_tfe_observed = as.integer(tmp_tfe$siteYearIndex_tfe)
data$n_siteYearIndex_tfe = tmp_tfe$n

# finish adding the ragged reference array to the list (damaged)
tmp_tot = tidybayes::compose_data(tmp_tot)

data$site_tot_observed = as.integer(tmp_tot$site2)
data$year_tot_observed = as.integer(tmp_tot$year2)
data$siteYearIndex_tot_observed = as.integer(tmp_tot$siteYearIndex_tot)
data$n_siteYearIndex_tot = tmp_tot$n

# check arrays
check.n = sample(1:139,1)
data$site_tfe_observed[check.n]
data$site[data$siteYearIndex_tfe==check.n]

check.n = sample(1:149,1)
data$site_tot_observed[check.n]
data$site2[data$siteYearIndex_tot==check.n]

# - +Save data for model fitting ----
saveRDS(data,file=paste0(tmpDataDirectory,"fruitsPerPlantAllPlants.RDS"))

# - +Save data for working with output ----
names(countFruitsPerPlantAllPlants)[1:2] = c('site','year')
names(countUndamagedDamagedFruitsPerPlantAllPlants)[1:2] = c('site','year')

ref_tfe<-countFruitsPerPlantAllPlants %>%
  dplyr::select(site,year,siteYearIndex_tfe) %>%
  unique

ref_tot<-countUndamagedDamagedFruitsPerPlantAllPlants %>%
  dplyr::select(site,year,siteYearIndex_tot) %>%
  unique

saveRDS(ref_tot,file=paste0(tmpDataDirectory,"referenceTOTFruitsPerPlant.RDS"))
saveRDS(ref_tfe,file=paste0(tmpDataDirectory,"referenceTFEFruitsPerPlant.RDS"))

