####
####
# Script to prepare data on seeds per fruit
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

countSeedPerFruit <- read.csv(paste0(dataDirectory,"countSeedPerFruit.csv"),header=TRUE)

# - Prep data for JAGS ----

# seeds per undamaged fruit
countSeedPerUndamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(damaged==0) %>%
  dplyr::rename(site3 = site) %>%
  dplyr::rename(year3 = year) %>%
  dplyr::select(site3,year3,sdno) %>% 
  dplyr::filter(!is.na(sdno))

countSeedPerUndamagedFruit$year3 <- as.character(countSeedPerUndamagedFruit$year3)

# create index for years with observations across the entire site
tmp <- countSeedPerUndamagedFruit %>%  
  dplyr::select(site3,year3) %>%
  unique
  
# add index
tmp_und <- cbind(tmp,siteYearIndex_und = (1:(dim(tmp)[1])))

# join the reference data frame tmp to the observations
countSeedPerUndamagedFruit<- countSeedPerUndamagedFruit %>%
  dplyr::left_join(tmp_und,by=c("site3","year3"))

# seeds per damaged fruit

countSeedPerDamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(damaged==1) %>%
  dplyr::rename(site4 = site) %>%
  dplyr::rename(year4 = year) %>%
  dplyr::rename(sdno_dam = sdno) %>%
  dplyr::select(site4,year4,sdno_dam) %>%
  dplyr::filter(!is.na(sdno_dam))

countSeedPerDamagedFruit$year4 <- as.character(countSeedPerDamagedFruit$year4)

# create index for years with observations across the entire site
tmp2 <- countSeedPerDamagedFruit %>%
  dplyr::select(site4,year4) %>%
  unique

# add index
tmp_dam <- cbind(tmp2,siteYearIndex_dam = (1:(dim(tmp2)[1])))

# join the reference data frame tmp to the observations
countSeedPerDamagedFruit<- countSeedPerDamagedFruit %>%
  dplyr::left_join(tmp_dam,by=c("site4","year4"))

# - +Data to list ----
# JAGS takes a list of data

data <- tidybayes::compose_data(countSeedPerUndamagedFruit,
                                countSeedPerDamagedFruit)

data$n3 = dim(countSeedPerUndamagedFruit)[1]
data$n4 = dim(countSeedPerDamagedFruit)[1]

og.size = object.size(data)
# reassign variables from numeric to integer, or remove 1D notation
# numeric is 8-byte; integer is 4-byte
data$site3 = as.integer(data$site3)
data$year3 = as.integer(data$year3)
data$sdno = as.integer(data$sdno)
data$siteYearIndex_und = as.integer(data$siteYearIndex_und)

data$site4 = as.integer(data$site4)
data$year4 = as.integer(data$year4)
data$sdno_dam = as.integer(data$sdno_dam)
data$siteYearIndex_dam = as.integer(data$siteYearIndex_dam)

new.size = object.size(data)
# the new object is 30% smaller!
new.size/og.size

# finish adding the ragged reference array to the list
tmp_und$siteYearIndex_und=as.integer(tmp_und$siteYearIndex_und)

tmp_und = tidybayes::compose_data(tmp_und)
data$site_und_observed = as.integer(tmp_und$site3)
data$year_und_observed = as.integer(tmp_und$year3)
data$siteYearIndex_und_observed = as.integer(tmp_und$siteYearIndex_und)
data$n_siteYearIndex_und = tmp_und$n

# damaged fruits
tmp_dam$siteYearIndex_dam=as.integer(tmp_dam$siteYearIndex_dam)

tmp_dam = tidybayes::compose_data(tmp_dam)
data$site_dam_observed = as.integer(tmp_dam$site4)
data$year_dam_observed = as.integer(tmp_dam$year4)
data$siteYearIndex_dam_observed = as.integer(tmp_dam$siteYearIndex_dam)
data$n_siteYearIndex_dam = tmp_dam$n

# - +Save data for model fitting ----
saveRDS(data,file=paste0(tmpDataDirectory,"seedsPerFruit.RDS"))

# - +Save data for working with output ----
names(countSeedPerUndamagedFruit)[1:2] = c('site','year')
names(countSeedPerDamagedFruit)[1:2] = c('site','year')

ref_und<-countSeedPerUndamagedFruit %>%
  dplyr::select(site,year,siteYearIndex_und) %>%
  unique

ref_dam<-countSeedPerDamagedFruit %>%
  dplyr::select(site,year,siteYearIndex_dam) %>%
  unique

df<-ref_und %>% dplyr::left_join(ref_dam,by=c('site','year'))

saveRDS(df,file=paste0(tmpDataDirectory,"referenceSeedsPerFruit.RDS"))