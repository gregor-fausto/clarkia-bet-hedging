####
####
# Script to summarize sample sizes of data used for analysis
# Produces txt files that are used to produce tables in the supplement
####
####

# - Libraries ----
library(tidyverse)
library(readxl)
library(knitr)
library(lubridate)
library(tidyxl)
library(readr)

# - Directories ----

dataDirectory <- "data/"
tableDirectory <- "outputs/008_sampleSizeSummaries/"

# - Prep site names for plotting ----

siteAbiotic<-read.csv(paste0(dataDirectory,"siteAbioticData.csv"),header=TRUE )

sitePosition <- siteAbiotic %>% dplyr::select(site,easting)
# get sitenames
# siteNames=unique(siteAbiotic$site)
# 
# # get easting
# position<-siteAbiotic %>% 
#   dplyr::select(site,easting) %>%
#   dplyr::mutate(easting=easting/1000)
# 
# # reorder by easting
# position=position[order(position$easting,decreasing=FALSE),]
# sitePosition=as.numeric(rownames(position))
# siteNames=position$site
# siteNames=(rev(siteNames))

# - Seed bag burial experiment ----

seedBagsData = read.csv(paste0(dataDirectory,"seedBagsData.csv"),header=TRUE) 

seedlingData <- seedBagsData %>%
  dplyr::filter(!is.na(totalJan)&!is.na(seedlingJan))

# - +Rename datasets ----
seedBagsDataIntactSeeds <- seedBagsData 
seedBagsDataSeedlings <- seedlingData

# - +Seedling counts ----

# seedling data (january)
seedBagsDataSeedlings<-seedBagsDataSeedlings  %>%
  tidyr::unite("yearAge", yearData,age,sep="-age") %>%
  dplyr::group_by(yearAge,site) %>%
  dplyr::summarise(count = sum(!is.na(seedlingJan))) %>%
  tidyr::spread(key=c("yearAge"),value="count") 

seedBagsDataSeedlings <- seedBagsDataSeedlings %>%
  dplyr::select(site,`2006-age1`,`2007-age1`,`2008-age1`,`2007-age2`,`2008-age2`,`2008-age3`)

seedBagsDataSeedlings <- seedBagsDataSeedlings %>%
  dplyr::left_join(sitePosition,by="site") %>%
  dplyr::arrange(easting) %>%
  dplyr::select(-easting)

kable(seedBagsDataSeedlings, caption="Summary table of the number of seed bags with seedling seed counts in January")
print(xtable::xtable(seedBagsDataSeedlings, type = "latex"), include.rownames=FALSE, NA.string = "--",
      file=paste0(tableDirectory,"seedBagsSeedlings.txt"))

# - +Intact seed counts ----

# january data
seedBagsDataIntactSeedsJan<-seedBagsDataIntactSeeds  %>%
  tidyr::unite("yearAge", yearData,age,sep="-age") %>%
  dplyr::group_by(yearAge,site) %>%
  dplyr::summarise(count = sum(!is.na(seedlingJan))) %>%
  tidyr::spread(key=c("yearAge"),value="count") 

seedBagsDataIntactSeedsJan <- seedBagsDataIntactSeedsJan %>%
  dplyr::select(site,`2006-age1`,`2007-age1`,`2008-age1`,`2007-age2`,`2008-age2`,`2008-age3`)

# october intacts
seedBagsDataIntactSeedsOct <-seedBagsDataIntactSeeds  %>%
  tidyr::unite("yearAge", yearData,age,sep="-age") %>%
  dplyr::group_by(yearAge,site) %>%
  dplyr::summarise(count = sum(!is.na(intactOct))) %>%
  tidyr::spread(key=c("yearAge"),value="count")

seedBagsDataIntactSeedsOct <- seedBagsDataIntactSeedsOct %>%
  dplyr::select(site,`2006-age1`,`2007-age1`,`2008-age1`,`2007-age2`,`2008-age2`,`2008-age3`)

summarySeedBagsIntactSeeds <- seedBagsDataIntactSeedsJan %>% 
  dplyr::right_join(seedBagsDataIntactSeedsOct, by=c("site")) %>% 
  mutate(`2006-age1` = pmax(`2006-age1.x`, `2006-age1.y`)) %>%
  mutate(`2007-age1` = pmax(`2007-age1.x`, `2007-age1.y`)) %>%
  mutate(`2008-age1` = pmax(`2008-age1.x`, `2008-age1.y`)) %>%
  mutate(`2007-age2` = pmax(`2007-age2.x`, `2007-age2.y`)) %>%
  mutate(`2008-age2` = pmax(`2008-age2.x`, `2008-age2.y`)) %>%
  mutate(`2008-age3` = pmax(`2008-age3.x`, `2008-age3.y`)) %>%
  dplyr::select(site,`2006-age1`,`2007-age1`,`2008-age1`,`2007-age2`,`2008-age2`,`2008-age3`)

summarySeedBagsIntactSeeds <- summarySeedBagsIntactSeeds %>%
  dplyr::left_join(sitePosition,by="site") %>%
  dplyr::arrange(easting) %>%
  dplyr::select(-easting)

kable(summarySeedBagsIntactSeeds, caption="Summary table of the number of seed bags with intact seed counts in January and/or October")
print(xtable::xtable(summarySeedBagsIntactSeeds, type = "latex"), include.rownames=FALSE, NA.string = "--",
      file=paste0(tableDirectory,"seedBagsIntactSeeds.txt"))


# - +Bags available for lab trials ----

seedBagsData = read.csv(paste0(dataDirectory,"seedBagsData.csv"),header=TRUE) 

seedlingData <- seedBagsData %>%
  dplyr::filter(!is.na(totalJan)&!is.na(seedlingJan))

# - +Rename datasets ----
seedBagsDataIntactSeeds <- seedBagsData 
seedBagsDataSeedlings <- seedlingData

# january data
seedBagsDataIntactSeedsOct<-seedBagsDataIntactSeeds  %>%
  tidyr::unite("yearAge", yearData,age,sep="-age") %>%
  dplyr::group_by(yearAge,site) %>%
  dplyr::summarise(count = sum(!is.na(intactOct))) %>%
  tidyr::spread(key=c("yearAge"),value="count") 

seedBagsDataIntactSeedsOct <- seedBagsDataIntactSeedsOct %>%
  dplyr::select(site,`2006-age1`,`2007-age1`,`2008-age1`,`2007-age2`,`2008-age2`,`2008-age3`)

seedBagsDataIntactSeedsOct <- seedBagsDataIntactSeedsOct %>%
  dplyr::left_join(sitePosition,by="site") %>%
  dplyr::arrange(easting) %>%
  dplyr::select(-easting)

kable(seedBagsDataIntactSeedsOct, caption="Summary table of the number of seed bags with intact seed counts in October")
print(xtable::xtable(seedBagsDataIntactSeedsOct, type = "latex"), include.rownames=FALSE, NA.string = "--",
      file=paste0(tableDirectory,"seedBagsIntactSeedsOct.txt"))


# - Viability trials ----

# - Seedlings and fruits per plant ----

# - +Sample size of seedling survival ----

censusSeedlingsFruitingPlants = read.csv(file=paste0(dataDirectory,"seedlingFruitingPlantCountsPermanentPlots.csv"),header=TRUE)

censusSeedlingsFruitingPlants <- censusSeedlingsFruitingPlants %>%
  dplyr::filter(year<2021)

censusSeedlingsFruitingPlantsTmp <- censusSeedlingsFruitingPlants %>%
  dplyr::filter(!(seedlingNumber==0&fruitplNumber==0)) %>%
  dplyr::group_by(year,site) %>%
  dplyr::summarise(count = sum(!is.na(seedlingNumber))) %>%
  tidyr::pivot_wider(names_from=year,values_from=count)

censusSeedlingsFruitingPlantsTmpTable <- censusSeedlingsFruitingPlantsTmp %>%
  dplyr::left_join(sitePosition,by="site") %>%
  dplyr::arrange(easting) %>%
  dplyr::select(-easting) 

kable(censusSeedlingsFruitingPlantsTmpTable, caption="Summary table of the number of estimates for seedling survival to fruiting")
print(xtable::xtable(censusSeedlingsFruitingPlantsTmpTable, type = "latex"), include.rownames=FALSE, NA.string = "--",
      file=paste0(tableDirectory,"censusSeedlingsFruitingPlants.txt"))

# - +Undercounting in seedling survival ----

censusSeedlingsFruitingPlants = read.csv(file=paste0(dataDirectory,"seedlingFruitingPlantCountsPermanentPlots.csv"),header=TRUE)

censusSeedlingsFruitingPlants <- censusSeedlingsFruitingPlants %>%
  dplyr::filter(year<2021)

undercountingSeedlingsFruitingPlants <- censusSeedlingsFruitingPlants %>%
  dplyr::filter(!(seedlingNumber==0&fruitplNumber==0)) %>%
  dplyr::group_by(year,site) %>%
  dplyr::filter(fruitplNumber<=seedlingNumber) %>%
  dplyr::summarise(count = sum(!is.na(seedlingNumber))) %>%
  tidyr::pivot_wider(names_from=year,values_from=count)

vals<-signif((1-undercountingSeedlingsFruitingPlants[,-1]/censusSeedlingsFruitingPlantsTmp[,-1])*100,2)

undercountingSeedlingsFruitingPlants = cbind(censusSeedlingsFruitingPlantsTmp[,1],vals)

undercountingSeedlingsFruitingPlants <- undercountingSeedlingsFruitingPlants %>%
  dplyr::left_join(sitePosition,by="site") %>%
  dplyr::arrange(easting) %>%
  dplyr::select(-easting)

kable(undercountingSeedlingsFruitingPlants, caption="Summary table of the percentage of plots with fruiting plant counts exceeding seedling counts")
print(xtable::xtable(undercountingSeedlingsFruitingPlants, type = "latex",digits=1), include.rownames=FALSE, NA.string = "--",
      file=paste0(tableDirectory,"undercountingSeedlingsFruitingPlants.txt"))

# - Fruits per plant from transects ----

# not run because not used in this analysis

# - +Total fruit equivalents from transects  ----
# countFruitsPerPlantTransects = readRDS(paste0(cleanDataDirectory,"countFruitsPerPlantTransects.RDS"))
# 
# countFruitsPerPlantTransects <- countFruitsPerPlantTransects %>%
#   dplyr::group_by(year,site) %>%
#   dplyr::summarise(count = sum(!is.na(countFruitsPerPlant))) %>%
#   tidyr::pivot_wider(names_from=year,values_from=count)
# 
# countFruitsPerPlantTransects <- arrange(countFruitsPerPlantTransects,countFruitsPerPlantTransects$site)
# 
# countFruitsPerPlantTransects <- countFruitsPerPlantTransects %>%
#   dplyr::left_join(sitePosition,by="site") %>%
#   dplyr::arrange(easting) %>%
#   dplyr::select(-easting)
# 
# kable(countFruitsPerPlantTransects, caption="Summary of dataset on total fruit equivalents per plant from transects")
# print(xtable::xtable(countFruitsPerPlantTransects, type = "latex"), include.rownames=FALSE, NA.string = "--",
#       file=paste0(tableDirectory,"countFruitsPerPlantTransects.txt"))
# 
# # - +Undamaged and damaged fruits from transects  ----
# countUndamagedDamagedFruitsPerPlantTransects = readRDS(paste0(cleanDataDirectory,"countUndamagedDamagedFruitsPerPlantTransects.RDS"))
# 
# countFruitsPerPlantTransects <- countFruitsPerPlantTransects %>%
#   dplyr::group_by(year,site) %>%
#   dplyr::summarise(count = sum(!is.na(countFruitsPerPlant))) %>%
#   tidyr::pivot_wider(names_from=year,values_from=count)
# 
# countFruitsPerPlantTransects <- arrange(countFruitsPerPlantTransects,countFruitsPerPlantTransects$site)
# 
# countFruitsPerPlantTransects <- countFruitsPerPlantTransects %>%
#   dplyr::left_join(sitePosition,by="site") %>%
#   dplyr::arrange(easting) %>%
#   dplyr::select(-easting)
# 
# kable(countFruitsPerPlantTransects, caption="Summary of dataset on total fruit equivalents per plant from transects")
# print(xtable::xtable(countFruitsPerPlantTransects, type = "latex"), include.rownames=FALSE, NA.string = "--",
#       file=paste0(tableDirectory,"countFruitsPerPlantTransects.txt"))

# - Fruits per plant from all plants ----

# - +Total fruit equivalents from all plants  ----
countFruitsPerPlantAllPlants <- read.csv(paste0(dataDirectory,"countFruitsPerPlantAllPlants.csv"),header=TRUE)

countFruitsPerPlantAllPlants <- countFruitsPerPlantAllPlants %>%
  #dplyr::filter(permanentPlot==0) %>%
  dplyr::group_by(year,site) %>%
  dplyr::summarise(count = sum(!is.na(countFruitNumberPerPlant))) %>%
  tidyr::pivot_wider(names_from=year,values_from=count)

countFruitsPerPlantAllPlants<-arrange(countFruitsPerPlantAllPlants,countFruitsPerPlantAllPlants$site)

countFruitsPerPlantAllPlants <- countFruitsPerPlantAllPlants %>%
  dplyr::left_join(sitePosition,by="site") %>%
  dplyr::arrange(easting) %>%
  dplyr::select(-easting)

kable(countFruitsPerPlantAllPlants, caption="Summary of dataset on total fruit equivalents per plant from all plants")
print(xtable::xtable(countFruitsPerPlantAllPlants, type = "latex"), include.rownames=FALSE, NA.string = "--",
      file=paste0(tableDirectory,"countFruitsPerPlantAllPlants.txt"))

# - +Undamaged and damaged fruits per plant from all plants  ----
countUndamagedDamagedFruitsPerPlantAllPlants <- read.csv(paste0(dataDirectory,"countUndamagedDamagedFruitsPerPlantAllPlants.csv"),header=TRUE) 

countUndamagedDamagedFruitsPerPlantAllPlants <- countUndamagedDamagedFruitsPerPlantAllPlants %>%
 # dplyr::filter(permanentPlot==0) %>%
  dplyr::group_by(year,site) %>%
  dplyr::summarise(count = sum(!is.na(countUndamagedFruitNumberPerPlant))) %>%
  tidyr::pivot_wider(names_from=year,values_from=count)

countUndamagedDamagedFruitsPerPlantAllPlants<-arrange(countUndamagedDamagedFruitsPerPlantAllPlants,countUndamagedDamagedFruitsPerPlantAllPlants$site)

countUndamagedDamagedFruitsPerPlantAllPlants <- countUndamagedDamagedFruitsPerPlantAllPlants %>%
  dplyr::left_join(sitePosition,by="site") %>%
  dplyr::arrange(easting) %>%
  dplyr::select(-easting)

kable(countUndamagedDamagedFruitsPerPlantAllPlants, caption="Summary of dataset on undamaged and damaged fruits per plant from all plants")
print(xtable::xtable(countUndamagedDamagedFruitsPerPlantAllPlants, type = "latex"), include.rownames=FALSE, NA.string = "--",
      file=paste0(tableDirectory,"countUndamagedDamagedFruitsPerPlantAllPlants.txt"))

# - Seeds per fruit ----

countSeedPerFruit <- read.csv(paste0(dataDirectory,"countSeedPerFruit.csv"),header=TRUE)

# - +Seeds per undamaged fruit ----

countSeedPerUndamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(damaged==0) %>%
  # dplyr::filter(demography==1) %>%
  dplyr::group_by(year,site) %>%
  dplyr::summarise(count = sum(!is.na(sdno))) %>%
  tidyr::pivot_wider(names_from=year,values_from=count)

countSeedPerUndamagedFruit <- countSeedPerUndamagedFruit %>%
  dplyr::left_join(sitePosition,by="site") %>%
  dplyr::arrange(easting) %>%
  dplyr::select(-easting)

kable(countSeedPerUndamagedFruit, caption="Summary of dataset on seeds per undamaged fruit")
print(xtable::xtable(countSeedPerUndamagedFruit, type = "latex"), include.rownames=FALSE, NA.string = "--",
      file=paste0(tableDirectory,"countSeedPerUndamagedFruit.txt"))

# - +Seeds per damaged fruit ----

countSeedPerDamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(damaged==1) %>%
 # dplyr::filter(demography==1) %>%
  dplyr::group_by(year,site) %>%
  dplyr::summarise(count = sum(!is.na(sdno))) %>%
  tidyr::pivot_wider(names_from=year,values_from=count)

countSeedPerDamagedFruit <- countSeedPerDamagedFruit %>%
  dplyr::left_join(sitePosition,by="site") %>%
  dplyr::arrange(easting) %>%
  dplyr::select(-easting)

kable(countSeedPerDamagedFruit, caption="Summary of dataset on seeds per damaged fruit")
print(xtable::xtable(countSeedPerDamagedFruit, type = "latex"), include.rownames=FALSE, NA.string = "--",
      file=paste0(tableDirectory,"countSeedPerDamagedFruit.txt"))
