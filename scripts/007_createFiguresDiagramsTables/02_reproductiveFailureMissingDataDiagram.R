# --
#  Script checks individual datasets to look for missing data values 
#  and create figure summarizing reproductive failure
# Produces Figure 1B
# --

# ---
# - Set up environment ----
# ---

# - +remove unused objects ----
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

# - Libraries ----
library(tidyverse)
library(khroma) # iridescent color scheme

# - Directories ----

dataDirectory = "data/"
figureDirectory <- "products/figures/"
outputDirectory <- "outputs/007_figuresDiagrams/"

# - Read data ----

seedlingFruitingPlantCountsPermanentPlots <- read.csv(paste0(dataDirectory,"seedlingFruitingPlantCountsPermanentPlots.csv"),header=TRUE)
countFruitsPerPlantAllPlants <- read.csv(paste0(dataDirectory,"countFruitsPerPlantAllPlants.csv"),header=TRUE)
countUndamagedDamagedFruitsPerPlantAllPlants <- read.csv(paste0(dataDirectory,"countUndamagedDamagedFruitsPerPlantAllPlants.csv"),header=TRUE)
countSeedPerFruit <- read.csv(paste0(dataDirectory,"countSeedPerFruit.csv"),header=TRUE)
siteAbiotic <- read.csv(paste0(dataDirectory,"siteAbioticData.csv"), header=TRUE )

# - Seedling survival to fruiting ----

# calculate proportion of seedlings that survive
sigmaData = seedlingFruitingPlantCountsPermanentPlots %>%
  dplyr::mutate( p = fruitplNumber/seedlingNumber )

# get total number of observations, number of NAs, and number of zeros or NA per site/year
sigmaNoObs = sigmaData %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na = sum(is.na(p)),
                   n.zero = sum(p==0|is.na(p)))

# get number of NAs and zeros
sigmaSummary = sigmaNoObs %>%
  dplyr::mutate(trueNA = ifelse(n.obs==n.na,1,0)) %>%
  dplyr::mutate(obsZero = ifelse(n.obs==n.zero,1,0))

# - Fruits per plant ----

# fruits per plant 2006-2012
countFruitsPerPlantAllPlants <- countFruitsPerPlantAllPlants %>%
  dplyr::rename(y_tfe = countFruitNumberPerPlant) %>%
  dplyr::select(site,year,y_tfe)
countFruitsPerPlantAllPlants$year <- as.character(countFruitsPerPlantAllPlants$year)

# fruits per plant 2013-2020
countUndamagedDamagedFruitsPerPlantAllPlants <- countUndamagedDamagedFruitsPerPlantAllPlants %>%
  dplyr::rename(y_und = countUndamagedFruitNumberPerPlant) %>%
  dplyr::rename(y_dam = countDamagedFruitNumberPerPlant) %>%
  dplyr::rename(site2 = site) %>%
  dplyr::rename(year2 = year) %>%
  dplyr::select(site2,year2,y_und,y_dam) 
countUndamagedDamagedFruitsPerPlantAllPlants$year2 <- as.character(countUndamagedDamagedFruitsPerPlantAllPlants$year2)

# summarize no observations for 2006-2012
fruitData = countFruitsPerPlantAllPlants

fruitNoObsOne = fruitData %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na = sum(is.na(y_tfe)),
                   n.zero = sum(y_tfe==0|is.na(y_tfe)))

fruitSummaryOne = fruitNoObsOne %>%
  dplyr::mutate(trueNA = ifelse(n.obs==n.na,1,0)) %>%
  dplyr::mutate(obsZero = ifelse(n.obs==n.zero,1,0)) %>%
  dplyr::select(site,year,n.obs,trueNA,obsZero)

# summarize no observations for 2013-2020
fruitDataTwo = countUndamagedDamagedFruitsPerPlantAllPlants %>%
  dplyr::rename(site=site2,year=year2)

fruitNoObsTwo = fruitDataTwo %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na1 = sum(is.na(y_und)),n.na2=sum(is.na(y_dam)),
                   n.zero1 = sum(y_und==0|is.na(y_und)),n.zero2 = sum(y_dam==0|is.na(y_dam)) )

fruitSummaryTwo = fruitNoObsTwo %>%
  dplyr::mutate(trueNA = ifelse(n.obs==((n.na1+n.na2)/2),1,0)) %>% 
  dplyr::mutate(obsZero = ifelse(n.obs==((n.zero1+n.zero2)/2),1,0)) %>%
  dplyr::select(site,year,n.obs,trueNA,obsZero)

# merge summaries
fruitSummary <- fruitSummaryOne %>%
  dplyr::bind_rows(fruitSummaryTwo)

# - Seeds per fruit ----

countSeedPerUndamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(damaged==0) %>%
  dplyr::select(site,year,sdno)

countSeedPerDamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(damaged==1) %>%
  dplyr::rename(sdno_dam = sdno) %>%
  dplyr::select(site,year,sdno_dam)

# summarize seeds per undamaged fruit
seedsUndamagedNoObs = countSeedPerUndamagedFruit %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na = sum(is.na(sdno)),
                   n.zero = sum(sdno==0|is.na(sdno)))

seedsUndamagedSummary = seedsUndamagedNoObs %>%
  dplyr::mutate(trueNA = ifelse(n.obs==n.na,1,0)) %>%
  dplyr::mutate(obsZero = ifelse(n.obs==n.zero,1,0)) %>%
  dplyr::select(site,year,n.obs,trueNA,obsZero)

# summarize seeds per damaged fruit
seedsDamagedNoObs = countSeedPerDamagedFruit %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na = sum(is.na(sdno_dam)),
                   n.zero = sum(sdno_dam==0|is.na(sdno_dam)))

seedsSummaryTwo = seedsDamagedNoObs %>%
  dplyr::mutate(trueNA = ifelse(n.obs==n.na,1,0)) %>%
  dplyr::mutate(obsZero = ifelse(n.obs==n.zero,1,0)) %>%
  dplyr::select(site,year,n.obs,trueNA,obsZero)

# join undamaged and damaged datasets
seedsSummary <- seedsSummaryTwo %>%
  dplyr::left_join(seedsUndamagedSummary,by=c('site','year')) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(trueNA=ifelse((trueNA.x+trueNA.y)==2,1,0),
                obsZero=ifelse((obsZero.x+obsZero.y)==2,1,0),
                n.obs = n.obs.x+n.obs.y) %>% 
  dplyr::select(site,year,n.obs,trueNA,obsZero) %>%
  dplyr::bind_rows(seedsUndamagedSummary %>% dplyr::filter(year<2013))

# - Prep site names for plotting ----

# get sitenames
siteNames=unique(siteAbiotic$site)

# get easting
position<-siteAbiotic %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

# reorder by easting
position=position[order(position$easting,decreasing=TRUE),]
siteNames=position$site

# ---
# - Manuscript figure ----
# ---


# - +set font sizes ----
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12

tiff(filename=paste0("products/figures/zero-fitness.tif"),
     units='px',height = 5600/3, width = 5400/2,res=800,pointsize=12)

par(mfrow=c(1,1),mar=c(.5,1,0,0),oma=c(0,.5,0,0)+.1,mgp=c(3,.1,0),xpd=TRUE)
plot(NA,NA,xlim=c(2006,2020),ylim=c(0,23),
     axes=FALSE,frame=FALSE,xaxt='n',yaxt='n',
     xlab='',ylab='')
polygon(x=c(2006,2007,2007,2006)+.5,
        y=c(0,0,21,21),
        col='gray97',border='gray97')
polygon(x=c(2008,2009,2009,2008)+.5,
        y=c(0,0,21,21),
        col='gray97',border='gray97')
polygon(x=c(2010,2011,2011,2010)+.5,
        y=c(0,0,21,21),
        col='gray97',border='gray97')
polygon(x=c(2012,2013,2013,2012)+.5,
        y=c(0,0,21,21),
        col='gray97',border='gray97')
polygon(x=c(2014,2015,2015,2014)+.5,
        y=c(0,0,21,21),
        col='gray97',border='gray97')
polygon(x=c(2016,2017,2017,2016)+.5,
        y=c(0,0,21,21),
        col='gray97',border='gray97')
polygon(x=c(2018,2019,2019,2018)+.5,
        y=c(0,0,21,21),
        col='gray97',border='gray97')

text(x=2007.5,y=21.5,cex=pt6,
     "Seedling survival")
text(x=2011,y=21.5,cex=pt6,
     "Fruits per plant")
text(x=2014.25,y=21.5,cex=pt6,
     "Seeds per fruit")

segments(x0=2009.1,x1=2011-.25,y0=21.15,y1=20.15,lty='solid')
segments(x0=2011,y0=21.25,y1=20)
segments(x0=2011.25,x1=2012.85,y0=20.15,y1=21.15,lty='solid')

y.pt=20:1
for(i in 20:1){
  tmp=sigmaSummary[sigmaSummary$site==siteNames[i],]
  points(tmp$year-.25,y=rep(y.pt[i],length(tmp$year)),
         pch=ifelse(tmp$trueNA==0&tmp$obsZero==1,21,ifelse(tmp$trueNA==0&tmp$obsZero!=1,21,4)),
         col=ifelse(tmp$obsZero==1,"#EECC66","black"),
         bg=ifelse(tmp$obsZero==1,"white","black"), 
         cex=.45)
  
  tmp=fruitSummary[fruitSummary$site==siteNames[i],]
  points(tmp$year,y=rep(y.pt[i],length(tmp$year)),
         pch=ifelse(tmp$trueNA==0,21,4),
         col=ifelse(tmp$obsZero==1,"#EECC66","black"),cex=.45,bg='black')
  
  tmp=seedsSummary[seedsSummary$site==siteNames[i],]
  points(tmp$year+.25,y=rep(y.pt[i],length(tmp$year)),
         pch=ifelse(tmp$trueNA==0,21,4),
         col=ifelse(tmp$obsZero==1,"#EECC66","black"),cex=.45,bg='black')
}

year.labs = c("2006","","2008","","2010","","2012","","2014","","2016","","2018","","2020")

axis(1,  at=2006:2020, labels = year.labs, 
     col.ticks = 1,las=1, 
     cex.axis = pt7, padj=-1,  tck=-0.01)
axis(2, (1:20), lwd.ticks = 0 ,
     labels = rev(siteNames), las = 2, 
     col = NA, col.ticks = 1, cex.axis = pt7)

legend(x=2004.5,y=24,
       legend = c("Observations of fitness component",
                  "0 seedlings survive in plots",
                  "No observations"),
       pch=c(19,21,4),
       col=c("black","#EECC66","#EECC66"),
       cex=pt6, pt.cex = .45,
       horiz=TRUE,
       bty='n',
      x.intersp=0.5,
      text.width=c(1,6.4,5.7),
      xjust = 0)



dev.off()