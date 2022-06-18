# --
#  Script checks individual datasets to look for missing data values 
# --

# ---
# - Set up environment ----
# ---

# - +remove unused objects ----
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)

# - +load packages ----
library(rjags) 
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(HDInterval)
library(bayesplot)
library(stringr)
library(khroma) # iridescent color scheme

# ---
# - Seedling survival to fruiting ----
# ---
censusSeedlingsFruitingPlants <- readRDS("~/Dropbox/dataLibrary/postProcessingData-2021/censusSeedlingsFruitingPlants.RDS")

censusSeedlingsFruitingPlants <- censusSeedlingsFruitingPlants %>% 
  # more fruiting plants than seedlings
  # recode s.t. number of seedlings equals number of fruiting plants
  dplyr::mutate(seedlingNumber=ifelse(fruitplNumber>seedlingNumber,fruitplNumber,seedlingNumber)) %>% 
  # NA for seedlings
  # filter these out; these are true missing data
  # recode s.t. number of seedlings equals number of fruiting plants
  #dplyr::mutate(seedlingNumber=ifelse(is.na(seedlingNumber),fruitplNumber,seedlingNumber)) %>%
  dplyr::filter(!is.na(seedlingNumber)) %>% 
  # NA for fruiting plants
  # still filter these out; missing response
  dplyr::filter(!is.na(fruitplNumber)) 

sigma.df = censusSeedlingsFruitingPlants %>%
  dplyr::mutate( p = fruitplNumber/seedlingNumber )# %>% 
 # dplyr::filter( p <= 1 | is.na(p))

sigma.noObs = sigma.df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na = sum(is.na(p)),
                   n.zero = sum(p==0|is.na(p)))

sigmaSummary = sigma.noObs %>%
  dplyr::mutate(trueNA = ifelse(n.obs==n.na,1,0)) %>%
  dplyr::mutate(obsZero = ifelse(n.obs==n.zero,1,0))

siteNames=unique(sigmaSummary$site)

plot(NA,NA,xlim=c(2006,2020),ylim=c(0,20),
     axes=FALSE,frame=FALSE,xaxt='n',yaxt='n',
     xlab='',ylab='')
y.pt=20:1
for(i in 20:1){
  tmp=sigmaSummary[sigmaSummary$site==siteNames[i],]
  points(tmp$year,y=rep(y.pt[i],length(tmp$year)),
         pch=ifelse(tmp$trueNA==0,19,4),
         col=ifelse(tmp$obsZero==1,"gray","black"))
}
axis(1,  2006:2020, col.ticks = 1,las=2)
axis(2, (1:20),
     labels = rev(siteNames), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)

zeroFitness = sigmaSummary %>% dplyr::filter(obsZero==1) %>%
  dplyr::mutate(site.year = paste0(site,year))
zeroFitness.df=sigma.df %>%
  dplyr::mutate(site.year = paste0(site,year)) %>%
  dplyr::filter(site.year %in% zeroFitness$site.year) %>% 
  dplyr::filter(p==0) 

f = function(x){
  max(seq(0,1,length.out=1000)[dbinom(0,x,seq(0,1,length.out=1000))>.5])
}

combos=unique(zeroFitness.df$site.year)
combos.list = list()
for(i in 1:length(combos)){
  tmp = zeroFitness.df[zeroFitness.df$site.year==combos[i],]
  combos.list[[i]] = sapply(tmp$seedlingNumber,f)
}
prob.sigma=lapply(combos.list, prod )



# ---
# - Fruits per plant ----
# ---
countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData-2021/countFruitsPerPlantAllPlots.RDS")

countFruitsPerPlantAllPlots <- countFruitsPerPlantAllPlots %>%
  dplyr::rename(y_tfe = countFruitNumberPerPlant) %>%
  dplyr::select(site,year,y_tfe)

countFruitsPerPlantAllPlots$year <- as.character(countFruitsPerPlantAllPlots$year)

countUndamagedDamagedFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData-2021/countUndamagedDamagedFruitsPerPlantAllPlots.RDS")

countUndamagedDamagedFruitsPerPlantAllPlots <- countUndamagedDamagedFruitsPerPlantAllPlots %>%
  dplyr::rename(y_und = countUndamagedFruitNumberPerPlant) %>%
  dplyr::rename(y_dam = countDamagedFruitNumberPerPlant) %>%
  dplyr::rename(site2 = site) %>%
  dplyr::rename(year2 = year) %>%
  dplyr::select(site2,year2,y_und,y_dam) 

countUndamagedDamagedFruitsPerPlantAllPlots$year2 <- as.character(countUndamagedDamagedFruitsPerPlantAllPlots$year2)

fec.df = countFruitsPerPlantAllPlots

fec.noObs = fec.df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na = sum(is.na(y_tfe)),
                   n.zero = sum(y_tfe==0|is.na(y_tfe)))

fecSummary1 = fec.noObs %>%
  dplyr::mutate(trueNA = ifelse(n.obs==n.na,1,0)) %>%
  dplyr::mutate(obsZero = ifelse(n.obs==n.zero,1,0)) %>%
  dplyr::select(site,year,n.obs,trueNA,obsZero)

fec.df = countUndamagedDamagedFruitsPerPlantAllPlots %>%
  dplyr::rename(site=site2,year=year2)

fec.noObs = fec.df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na1 = sum(is.na(y_und)),n.na2=sum(is.na(y_dam)),
                   n.zero1 = sum(y_und==0|is.na(y_und)),n.zero2 = sum(y_dam==0|is.na(y_dam)) )

fecSummary2 = fec.noObs %>%
  dplyr::mutate(trueNA = ifelse(n.obs==((n.na1+n.na2)/2),1,0)) %>% 
  dplyr::mutate(obsZero = ifelse(n.obs==((n.zero1+n.zero2)/2),1,0)) %>%
  dplyr::select(site,year,n.obs,trueNA,obsZero)

fecSummary <- fecSummary1 %>%
  dplyr::bind_rows(fecSummary2)

siteNames=unique(fecSummary$site)

plot(NA,NA,xlim=c(2006,2020),ylim=c(0,20),
     axes=FALSE,frame=FALSE,xaxt='n',yaxt='n',
     xlab='',ylab='')
y.pt=20:1
for(i in 20:1){
  tmp=fecSummary[fecSummary$site==siteNames[i],]
  points(tmp$year,y=rep(y.pt[i],length(tmp$year)),
         pch=ifelse(tmp$trueNA==0,19,4),
         col=ifelse(tmp$obsZero==1,"gray","black"))
}
axis(1,  2006:2020, col.ticks = 1,las=2)
axis(2, (1:20),
     labels = rev(siteNames), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)

# ---
# - Seeds per fruit ----
# ---

countSeedPerFruit <- readRDS("~/Dropbox/dataLibrary/postProcessingData-2021/countSeedPerFruit.RDS")

countSeedPerUndamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(damaged==0) %>%
  dplyr::select(site,year,sdno)

countSeedPerDamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(damaged==1) %>%
  dplyr::rename(sdno_dam = sdno) %>%
  dplyr::select(site,year,sdno_dam)

seeds.df = countSeedPerUndamagedFruit

seeds.noObs = seeds.df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na = sum(is.na(sdno)),
                   n.zero = sum(sdno==0|is.na(sdno)))

seedsSummary1 = seeds.noObs %>%
  dplyr::mutate(trueNA = ifelse(n.obs==n.na,1,0)) %>%
  dplyr::mutate(obsZero = ifelse(n.obs==n.zero,1,0)) %>%
  dplyr::select(site,year,n.obs,trueNA,obsZero)

countSeedPerDamagedFruit <- countSeedPerFruit %>%
  dplyr::filter(demography==1) %>%
  dplyr::filter(damaged==1) %>%
  dplyr::rename(sdno_dam = sdno) %>%
  dplyr::select(site,year,sdno_dam)

seeds.df = countSeedPerDamagedFruit

seeds.noObs = seeds.df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(n.obs = n(),
                   n.na = sum(is.na(sdno_dam)),
                   n.zero = sum(sdno_dam==0|is.na(sdno_dam)))

seedsSummary2 = seeds.noObs %>%
  dplyr::mutate(trueNA = ifelse(n.obs==n.na,1,0)) %>%
  dplyr::mutate(obsZero = ifelse(n.obs==n.zero,1,0)) %>%
  dplyr::select(site,year,n.obs,trueNA,obsZero)

# join undamaged and damaged datasets
seedsSummary <- seedsSummary2 %>%
  dplyr::left_join(seedsSummary1,by=c('site','year')) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(trueNA=ifelse((trueNA.x+trueNA.y)==2,1,0),
               obsZero=ifelse((obsZero.x+obsZero.y)==2,1,0),
               n.obs = n.obs.x+n.obs.y) %>% 
  dplyr::select(site,year,n.obs,trueNA,obsZero) %>%
  dplyr::bind_rows(seedsSummary1 %>% dplyr::filter(year<2013))

siteNames=unique(seedsSummary1$site)

plot(NA,NA,xlim=c(2006,2020),ylim=c(0,20),
     axes=FALSE,frame=FALSE,xaxt='n',yaxt='n',
     xlab='',ylab='')
y.pt=20:1
for(i in 20:1){
  tmp=seedsSummary[seedsSummary$site==siteNames[i],]
  points(tmp$year,y=rep(y.pt[i],length(tmp$year)),
         pch=ifelse(tmp$trueNA==0,19,4),
         col=ifelse(tmp$obsZero==1,"gray","black"))
}
axis(1,  2006:2020, col.ticks = 1,las=2)
axis(2, (1:20),
     labels = rev(siteNames), las = 2, 
     col = NA, col.ticks = 1, cex.axis = 1)

# ---
# - Save data object with missing data/zero fitness years ----
# ---
# Write out data objects with missing data/zero fitness years
# Specifically no plants 

missing.list = list()
for(i in 1:20){
  tmp.sigma=sigmaSummary[sigmaSummary$site==siteNames[i],]
  tmp.fec=fecSummary[fecSummary$site==siteNames[i],]
  tmp.seeds=seedsSummary[seedsSummary$site==siteNames[i],]
  tmp.seeds= tmp.seeds[ order( tmp.seeds$year),]
  na.obs=tmp.sigma$trueNA==1&tmp.fec$trueNA==1&tmp.seeds$trueNA==1
  na.obs2=tmp.sigma$obsZero==1&tmp.fec$trueNA==1&tmp.seeds$trueNA==1
  tmp=tmp.sigma[na.obs|na.obs2,1:2]
  missing.list[[i]] = tmp
}
lowFitnessYears=do.call(rbind,missing.list)
saveRDS(lowFitnessYears,"/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/007_hypothesisTesting/objects/lowFitnessYears.RDS")


# ---
# - Reorder site by position ----
# ---
position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

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

tiff(filename=paste0("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/products/figures/zero-fitness.tif"),
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
  
  tmp=fecSummary[fecSummary$site==siteNames[i],]
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


# 
# missing.list = list()
# for(i in 1:20){
#   tmp.sigma=sigmaSummary[sigmaSummary$site==siteNames[i],]
#   tmp.fec=fecSummary[fecSummary$site==siteNames[i],]
#   tmp.seeds=seedsSummary[seedsSummary$site==siteNames[i],]
#   tmp.seeds= tmp.seeds[ order( tmp.seeds$year),]
#   na.obs=tmp.sigma$trueNA==1
#   na.obs2=tmp.sigma$obsZero==1
#   tmp=tmp.sigma[na.obs|na.obs2,1:2]
#   missing.list[[i]] = tmp
# }
# lowFitnessYearsPlots=do.call(rbind,missing.list)
# saveRDS(lowFitnessYearsPlots,"/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/output/lowFitnessYearsPlots.RDS")
