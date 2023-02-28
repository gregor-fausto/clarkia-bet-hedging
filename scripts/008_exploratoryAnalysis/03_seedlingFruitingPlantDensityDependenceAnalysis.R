#################################################################################
# Script to conduct exploratory analysis of density-dependence in seedling survival to fruiting
# we conducted the analysis in this script after the first round of peer review
# a reviewer suggested we consider whether there was density-dependence in seedling establishment
# we followed methods from:
# Detto, M., M. D. Visser, S. J. Wright, and S. W. Pacala. 2019. Bias in the detection of negative density dependence in plant communities. Ecology Letters 22:1923–1939.
# and use code provided with the paper to conduct the analysis 
#################################################################################

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

dataDirectory = "~/Dropbox/clarkia-demography-projects/data/"
tmpDataDirectory = "outputs/001_prepareDataForModels/"

# - Libraries ----
library(tidybayes)
library(tidyverse)
library(stringr)
library(khroma)

# - Read in data ----

seedlingFruitingPlantCountsPermanentPlots <- read.csv(paste0(dataDirectory,"seedlingFruitingPlantCountsPermanentPlots.csv"),header=TRUE)
countFruitsPerPlantFromPermanentPlots <- read.csv(paste0(dataDirectory,"countFruitsPerPlantFromPermanentPlots.csv"),header=TRUE)
countUndamagedDamagedFruitsPerPlantFromPermanentPlots <- read.csv(paste0(dataDirectory,"countUndamagedDamagedFruitsPerPlantFromPermanentPlots.csv"),header=TRUE)

# - Read in posterior estimates ----
source("scripts/006_testHypotheses/00_utilityFunctions.R")

fec <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/combinedF-population-year-level-mat.RDS")

vr.popYear <- fec

df.list <- list()
for(i in 1:20){
  tmp<-vr.popYear[[i]][,1:15]
  df.list[[i]]<-data.frame(site=i,
                           year = c(2006:2020),
                           median=apply(tmp,2,median),
                           mode = apply(tmp,2,posterior.mode),
                           mean = apply(tmp,2,mean),
                           ci.lo = apply(tmp,2,function(x) hdi(x,credMass=0.95))[1,], 
                           ci.hi = apply(tmp,2,function(x) hdi(x,credMass=0.95))[2,],
                           ci.lo2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[1,], 
                           ci.hi2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[2,])
}
df.summary<-do.call(rbind,df.list)
siteAbiotic <- read.csv("~/Dropbox/clarkia-demography-projects/data/siteAbioticData.csv",header=TRUE)
siteDf <- data.frame(site=1:20,siteName=siteAbiotic$site)
df.summary<-df.summary %>% dplyr::left_join(siteDf,by=c('site')) %>%
  dplyr::select(-site) %>% dplyr::rename(site=siteName)


# - Tally data ----

countFruitsPerPlantFromPermanentPlots.sum<-countFruitsPerPlantFromPermanentPlots %>%
  dplyr::group_by(site,year,transect,position) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::left_join(seedlingFruitingPlantCountsPermanentPlots,by=c('site','year','transect','position')) %>%
  dplyr::mutate(delta=n-fruitplNumber) 

df1<-countFruitsPerPlantFromPermanentPlots %>%
  dplyr::group_by(site,year,transect,position) %>%
  dplyr::summarise(nPlants = n(), totalFruits = sum(countFruitsPerPlant)) %>%
  dplyr::left_join(seedlingFruitingPlantCountsPermanentPlots,by=c('site','year','transect','position')) %>%
  dplyr::mutate(delta=nPlants-fruitplNumber)

df2<-df1 %>%
  dplyr::filter(delta<0) %>%
  dplyr::left_join(df.summary,by=c('site','year')) %>%
  dplyr::mutate(netFruits = totalFruits + (-1*delta)*mode) %>%
  dplyr::select(site,year,transect,position,nPlants,seedlingNumber,fruitplNumber,netFruits,delta) %>%
  dplyr::rename(totalFruits=netFruits)

fruitNumber<-df1 %>%
  dplyr::filter(delta==0) %>%
  dplyr::bind_rows(df2) 

# second round
countUndamagedDamagedFruitsPerPlantFromPermanentPlots.sum<-countUndamagedDamagedFruitsPerPlantFromPermanentPlots %>%
  dplyr::group_by(site,year,transect,position) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::left_join(seedlingFruitingPlantCountsPermanentPlots,by=c('site','year','transect','position')) %>%
  dplyr::mutate(delta=n-fruitplNumber) 

df1.1<-countUndamagedDamagedFruitsPerPlantFromPermanentPlots %>%
  dplyr::group_by(site,year,transect,position) %>%
  dplyr::summarise(nPlants = n(), 
                   totalUndamagedFruits = sum(countUndamagedFruitsPerPlant),
                   totalDamagedFruits = sum(countDamagedFruitsPerPlant)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(seedlingFruitingPlantCountsPermanentPlots,by=c('site','year','transect','position')) %>%
  dplyr::mutate(delta=nPlants-fruitplNumber)

df2.1<-df1.1 %>%
  dplyr::filter(delta<0) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(df.summary,by=c('site','year')) %>%
  dplyr::mutate(netFruits =  (-1*delta)*mode) %>%
  dplyr::select(site,year,transect,position,nPlants,seedlingNumber,fruitplNumber,totalUndamagedFruits,totalDamagedFruits,netFruits,delta)



# - Read in posterior estimates ----
source("scripts/006_testHypotheses/00_utilityFunctions.R")

phi <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/phi-population-year-level-mat.RDS")
phiDamage <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/phi-damage-population-year-level-mat.RDS")

vr.popYear <- phi

df.list <- df.list2 <- list()
for(i in 1:20){
  tmp<-vr.popYear[[i]][,1:15]
  df.list[[i]]<-data.frame(site=i,
                           year = c(2006:2020),
                           mode = apply(tmp,2,posterior.mode))
  tmp2<-phiDamage[[i]][,1:15]
  df.list2[[i]]<-data.frame(site=i,
                            year = c(2006:2020),
                            mode.dam = apply(boot::inv.logit(tmp2),2,posterior.mode))
  
}
df.summary<-do.call(rbind,df.list) %>%
  dplyr::select(site,year,mode)
df.summary2<-do.call(rbind,df.list2) %>%
  dplyr::select(site,year,mode.dam)
siteAbiotic <- read.csv("~/Dropbox/clarkia-demography-projects/data/siteAbioticData.csv",header=TRUE)
siteDf <- data.frame(site=1:20,siteName=siteAbiotic$site)
df.summary<-df.summary %>% dplyr::left_join(siteDf,by=c('site')) %>%
  dplyr::select(-site) %>% dplyr::rename(site=siteName) 
df.summary2<-df.summary2 %>% dplyr::left_join(siteDf,by=c('site')) %>%
  dplyr::select(-site) %>% dplyr::rename(site=siteName) 
# first round join

tmp.df<-fruitNumber %>%
  dplyr::left_join(df.summary,by=c('site','year')) %>%
  dplyr::mutate(seedRain=totalFruits*mode)

seedRain <- tmp.df %>%
  dplyr::ungroup() %>%
  dplyr::select(site,year,transect,position,seedRain)

# second round join

seedRain2<-df1.1 %>%
  dplyr::filter(delta==0) %>%
  dplyr::left_join(df.summary,by=c('site','year')) %>%
  dplyr::left_join(df.summary2,by=c('site','year')) %>%
  dplyr::mutate(seedRain = totalUndamagedFruits*mode+totalDamagedFruits*mode*mode.dam) %>%
  dplyr::ungroup() %>%
  dplyr::select(site,year,transect,position,seedRain)

seedRain3<-df2.1 %>%
  dplyr::left_join(df.summary,by=c('site','year')) %>%
  dplyr::left_join(df.summary2,by=c('site','year')) %>%
  dplyr::mutate(seedRain = (totalUndamagedFruits+netFruits)*mode+totalDamagedFruits*mode*mode.dam) %>%
  dplyr::ungroup() %>%
  dplyr::select(site,year,transect,position,seedRain)

seedRain2 <- seedRain2 %>%
  bind_rows(seedRain3)

seedRain <- seedRain %>%
  bind_rows(seedRain2)
# net
obs.list<-list()
years=2010:2020
for(i in 1:10){
  
  current.year = years[i]-c(1,2,3)
  t1<-seedRain %>%
    dplyr::filter(year==current.year[1]) %>%
    dplyr::select(-c(year)) %>%
    dplyr::mutate(seedRain=ifelse(is.na(seedRain),0,seedRain)) %>%
    dplyr::rename(y_1=seedRain) 
  t2<-seedRain %>%
    dplyr::filter(year==current.year[2]) %>%
    dplyr::select(-year) %>%
    dplyr::mutate(seedRain=ifelse(is.na(seedRain),0,seedRain)) %>%
    dplyr::rename(y_2=seedRain)
  t3<-seedRain %>%
    dplyr::filter(year==current.year[3]) %>%
    dplyr::select(-year) %>%
    dplyr::mutate(seedRain=ifelse(is.na(seedRain),0,seedRain)) %>%
    dplyr::rename(y_3=seedRain) 
  seedRainHistory<-t1 %>%
    dplyr::full_join(t2,by=c("site","transect","position")) %>%
    dplyr::full_join(t3,by=c("site","transect","position")) %>%
    dplyr::mutate(year=years[i]) %>%
    dplyr::mutate(y_3=ifelse(is.na(y_3),0,y_3)) %>%
    dplyr::mutate(y_2=ifelse(is.na(y_2),0,y_2)) %>%
    dplyr::mutate(y_1=ifelse(is.na(y_1),0,y_1))
  
  
  seedlingObservations<- seedlingFruitingPlantCountsPermanentPlots %>% 
    dplyr::select(site,year,transect,position,seedlingNumber)
  seedRainSeedlingHistory<-seedRainHistory %>%
    dplyr::left_join(seedlingObservations,by=c('site','transect','position','year')) 
  obs.list[[i]]<-seedRainSeedlingHistory  
}
seedHistory<-do.call(rbind,obs.list)

# - Read in posterior estimates ----
s0 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s0-population-level.RDS")
g1 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/g1-population-level.RDS")
s1 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s1-population-level.RDS")
s2 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s2-population-level.RDS")
s3 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s3-population-level.RDS")
g1.hat  <- apply(g1,2,posterior.mode)
s1.hat  <- apply(s1,2,posterior.mode)
s2.hat  <- apply(s2,2,posterior.mode)
s3.hat  <- apply(s3,2,posterior.mode)
s0.hat  <- apply(s0,2,posterior.mode)

vr.dat <- data.frame(site=siteAbiotic$site,g1=g1.hat,s0=s0.hat,s1=s1.hat,s2=s2.hat,s3=s3.hat)
vr.summary=vr.dat %>%
  dplyr::mutate(transitionOne=s0*s1,
                transitionTwo=s0*s1*(1-g1)*s2*s3,
                transitionThree=s0*s1*(1-g1)*s2*s3*(1-g1)*s2*s3) %>%
  dplyr::select(site,transitionOne,transitionTwo,transitionThree,g1)
out=seedHistory %>%
  dplyr::left_join(vr.summary,by='site')
head(out)
df<-out %>%
  dplyr::mutate(expectSeedBank=y_1*transitionOne+y_2*transitionTwo+y_3*transitionThree) %>%
  dplyr::mutate(expectSeedlings=expectSeedBank*g1)

plot(df$expectSeedBank,df$seedlingNumber,xlim=c(0,6500),ylim=c(0,6500))
abline(a=0,b=1)

ggplot(data=df ,aes(x=expectSeedBank,y=seedlingNumber,color=as.factor(site))) +
  geom_point() +
  facet_wrap(~site,scales='free') +
  geom_abline(intercept=0,slope=1) +
  scale_color_batlow(reverse=TRUE,discrete=TRUE) +
  xlab("Expected seeds in 3-year seed bank (January year t)") +
  ylab("Observed seedlings (year t)") +
  labs(color="Population") +
  theme_bw() +
  theme(panel.spacing=unit(1,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(
    legend.background = element_rect(fill="gray95"),
    legend.key = element_rect(fill="gray95")) 

# run regression
df.reg<-df %>%
  dplyr::select(site,year,transect,position,expectSeedBank,seedlingNumber) %>%
  dplyr::mutate(seedBank=round(expectSeedBank)) %>%
  dplyr::mutate(seedBank=ifelse(seedBank<seedlingNumber,seedlingNumber,seedBank)) 

siteNames=siteAbiotic$site
out.list<-list()
for(i in 1:20){
  df.tmp<-df.reg %>%
    dplyr::filter(site==siteNames[i])
  m1 <- lm(log(seedlingNumber+1)~log(seedBank+1),data=df.tmp)
  out.list[[i]]<-c(coef(m1)[2],confint(m1,'log(seedBank + 1)',level=0.95))
}

rbinom2 <- function(n,S,p){
  
  # random number generator of S trials given R number success with Bernoulli
  # probability p
  z = runif(n)
  R=numeric(length(S))
  
  use=S==0
  R[use] <- 0
  
  S0=unique(S[S>0])
  glS1 <- lgamma(S0+1)
  if (p>0 & p<1){
    for (i in 1:length(S0)){
      
      use2=S==S0[i]
      k<-seq.int(0,S0[i])
      C <- cumsum(exp(glS1[i] - lgamma(k+1) -lgamma(S0[i]-k+1) + k*log(p) + (S0[i]-k)*log(1-p)))
      R[use2] <- .bincode(z[use2], c(0,C))-1
      
    }
  } else {
    R = S*p
  }
  
  return(R)
}

H1 <-  data.frame(b1=numeric(length(20)),
                  lb1=numeric(length(20)),
                  ub1=numeric(length(20)))
f <- numeric()
for(i in 1:20){
  df.tmp<-df.reg %>%
    dplyr::filter(site==siteNames[i])
  f[i] = sum(df.tmp$seedlingNumber)/sum(df.tmp$seedBank)
  Rep<-1000
  mod.rnd1 <- numeric()
  
  for(j in 1:Rep){
    cat("\r",i,":",Rep-j,"\r")
    
    # null model 1
    R.r = rbinom2(length(df.tmp$seedBank),df.tmp$seedBank,f[i])
    seedlingNumber.sim = R.r  
    seedBank.tmp = df.tmp$seedBank
    mod.rnd1[j] <- coef(lm(log(seedlingNumber.sim+1)~log(seedBank.tmp+1)))[2]
    
  }
  
  
  H1[i,] <- c(mean(mod.rnd1,na.rm=TRUE),quantile(mod.rnd1,prob=c(0.025,0.975),na.rm=TRUE))
  
}

out.mat<-do.call(rbind,out.list)

pdf(file=paste0("products/figures/densityDependence-seedSeedling.pdf"),
    height=6,width=6,pointsize=12)

par(mar = c(6.1, 4.6, 4.1, 4.1), # change the margins
    lwd = 1, # increase the line thickness
    cex.axis = 1.2,
    xpd=F# increase default axis label size
)
index=order(siteAbiotic$easting)
plot(1:20-.1,out.mat[index,1]-1,ylim=c(-1.2,0),pch=NA,
     xaxt="n",xlab="",ylab=expression(italic(hat(b)) - italic(b)['NULL']))

rect(xleft=c(1,3,5,7,9,11,13,15,17,19)+c(.5),xright=c(1,3,5,7,9,11,13,15,17,19)+1.5,ybottom=-100,ytop=100,col='gray99',border='gray99')

points(1:20-.1,out.mat[index,1]-1,ylim=c(-1.2,0),pch=16)
abline(h=0,lty='dotted')
segments(x0=1:20-.1,y0=out.mat[index,2]-1,y1=out.mat[index,3]-1)
points(x=1:20+.1,y=out.mat[index,1]-H1$b1[index],pch=16,cex=1,col='red')
segments(y0=out.mat[index,1]-H1$ub1[index],x0=1:20+.1,y1=out.mat[index,1]-H1$lb1[index],col='red')
axis(1,at=1:20,labels=siteNames[index],las=2)
box()
## Draw the x-axis labels.
par(xpd=T)
legend(0.05, 0.25,
       c("Estimate from offset power-law model", 
         "Estimate from offset power-law model with correction for offset"),
       col = c("black", "red"),
       cex = 0.8,
       pch=16,lty=1)

dev.off()

d<-df %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(cor.out=cor(expectSeedBank,seedlingNumber)) %>%
  dplyr::left_join(siteAbiotic,by='site') 
d$site = with(d, reorder(site, easting))

ggplot(d,aes(x=site,y=cor.out,fill=site)) + geom_boxplot()+
  scale_fill_batlow(reverse=TRUE,discrete=TRUE)

ggplot(d,aes(x=year,y=cor.out,color=site)) + geom_point()+geom_line() +
  scale_color_batlow(reverse=TRUE,discrete=TRUE)

countUndamagedDamagedFruitsPerPlantFromPermanentPlots.sum<-countUndamagedDamagedFruitsPerPlantFromPermanentPlots %>%
  dplyr::group_by(site,year,transect,position) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::left_join(seedlingFruitingPlantCountsPermanentPlots,by=c('site','year','transect','position')) %>%
  dplyr::mutate(delta=n-fruitplNumber) 


