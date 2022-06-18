
# - Libraries ----
library(tidyverse)
library(parallel)
library(stringr)
library(bbmle)

# - Read in estimates ----

estimates <- readRDS(paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","seedSamplesLikelihood.rds"))
modelList <- readRDS(paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","seedSamplesLikelihoodModels.rds"))

# - Read in additional data ----

climate <- readRDS("~/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
position<-climate %>% 
  dplyr::ungroup() %>%
  dplyr::select(site,easting,intenseDemography) %>%
  dplyr::mutate(easting=easting/1000)

position = position %>% dplyr::filter(intenseDemography==1) %>% unique
estimates <- estimates %>% dplyr::left_join(data.frame(siteName=position$site,site = 1:20),by='site')

# - Plot estimates ----

# - +germination ----

# plot estimates for germination
plot(1:6,rep(NA,6),ylim=c(0,1))

for(i in 1:20){
  m.tmp = modelList[[i]]
  points(1:6,m.tmp@coef[1:6])
}

# plot interannual estimates

par(mfrow=c(1,3))

# - ++g1 ----

plot(rep(NA,3),c(1,2,3),ylim=c(0,1),xlim=c(340,375),xlab="Easting",ylab="Probability of germination")

for(i in 1:20){
  
  m.tmp = modelList[[i]]
  points(rep(position$easting[i],3),m.tmp@coef[c(1,4,6)],pch=16,col=c("black","red","pink"))
  points(c(position$easting[i]),mean(m.tmp@coef[c(1,4,6)]),pch=5,col=c("black","red","pink"))
}

legend("topleft",c("Round 1","Round 2", "Round 3"),pch=16,col=c("black","red","pink"))

# - ++g2 ----

plot(rep(NA,3),c(1,2,3),ylim=c(0,1),xlim=c(340,375),xlab="Easting",ylab="Probability of germination")

for(i in 1:20){
  
  m.tmp = modelList[[i]]
  points(rep(position$easting[i],2),m.tmp@coef[c(2,5)],pch=16,col=c("black","red","pink"))
}

legend("topleft",c("Round 1","Round 2"),pch=16,col=c("black","red"))


# - ++g3 ----

plot(rep(NA,3),c(1,2,3),ylim=c(0,1),xlim=c(340,375),xlab="Easting",ylab="Probability of germination")

for(i in 1:20){
  
  m.tmp = modelList[[i]]
  points(rep(position$easting[i],1),m.tmp@coef[c(3)],pch=16,col=c("black","red","pink"))
}

legend("topleft",c("Round 1"),pch=16,col=c("black"))

# - ++g1: time ----

par(mfrow=c(1,1))
plot(rep(NA,3),c(1,2,3),ylim=c(0,.6),xlim=c(1,20),xlab="Population",ylab="Probability of germination")

orderedIndex = order(position$easting,decreasing=FALSE)

for(i in 1:20){
  
  m.tmp = modelList[[orderedIndex[i]]]
  points(c(i-.25,i,i+.25),m.tmp@coef[c(1,4,6)],pch=16,col=c("black","red","pink"),type='o')
  points(c(i),mean(m.tmp@coef[c(1,4,6)]),pch=5,col=c("black"))
  
}

abline(v=seq(.5,20.5,by=1),lty='dotted')

legend("topleft",c("Round 1","Round 2", "Round 3"),pch=16,col=c("black","red","pink"))



# - +survival ----

# plot estimates for survival
plot(1:12,rep(NA,12),ylim=c(0,1))

for(i in 1:20){
  m.tmp = modelList[[i]]
  points(1:12,m.tmp@coef[7:18])
}

# year 1
plot(rep(NA,1),c(NA),ylim=c(0,1),xlim=c(340,375))

for(i in 1:20){
  
  m.tmp = modelList[[i]]
  m.surv = m.tmp@coef[c(7:18)]
  points(rep(position$easting[i],3),m.surv[c(1,7,11)],pch=16,col=c("black","red","pink"))
}

par(mfrow=c(1,1))
plot(rep(NA,3),c(1,2,3),ylim=c(0,1),xlim=c(1,20),xlab="Population",ylab="Probability of s1")

for(i in 1:20){
  
  m.tmp = modelList[[orderedIndex[i]]]
  m.surv = m.tmp@coef[c(7:18)]
  points(c(i-.25,i,i+.25),m.surv[c(1,7,11)],pch=16,col=c("black","red","pink"),type='o')
}
abline(v=seq(.5,20.5,by=1),lty='dotted')

legend("bottomleft",c("Round 1","Round 2", "Round 3"),pch=16,col=c("black","red","pink"))


# year 1
plot(rep(NA,1),c(NA),ylim=c(0,1),xlim=c(340,375))

for(i in 1:20){
  
  m.tmp = modelList[[i]]
  m.surv = m.tmp@coef[c(7:18)]
  points(rep(position$easting[i],3),m.surv[c(1,7,11)+1],pch=16,col=c("black","red","pink"))
}


par(mfrow=c(1,1))
plot(rep(NA,3),c(1,2,3),ylim=c(0,1),xlim=c(1,20),xlab="Population",ylab="Probability of s2")

for(i in 1:20){
  
  m.tmp = modelList[[orderedIndex[i]]]
  m.surv = m.tmp@coef[c(7:18)]
  points(c(i-.25,i,i+.25),m.surv[c(1,7,11)+1],pch=16,col=c("black","red","pink"))
}
abline(v=seq(.5,20.5,by=1),lty='dotted')

legend("bottomleft",c("Round 1","Round 2", "Round 3"),pch=16,col=c("black","red","pink"))



# year 2
plot(rep(NA,1),c(NA),ylim=c(0,1),xlim=c(340,375))

for(i in 1:20){
  
  m.tmp = modelList[[i]]
  m.surv = m.tmp@coef[c(7:18)]
  points(rep(position$easting[i],2),m.surv[c(1,7)+2],pch=16,col=c("black","red","pink"))
}


par(mfrow=c(1,1))
plot(rep(NA,3),c(1,2,3),ylim=c(0,1),xlim=c(1,20),xlab="Population",ylab="Probability of s3")

for(i in 1:20){
  
  m.tmp = modelList[[orderedIndex[i]]]
  m.surv = m.tmp@coef[c(7:18)]
  points(c(i-.25,i),m.surv[c(1,7)+2],pch=16,col=c("black","red","pink"))
}
abline(v=seq(.5,20.5,by=1),lty='dotted')

legend("bottomleft",c("Round 1","Round 2", "Round 3"),pch=16,col=c("black","red","pink"))



# year 2
plot(rep(NA,1),c(NA),ylim=c(0,1),xlim=c(340,375))

for(i in 1:20){
  
  m.tmp = modelList[[i]]
  m.surv = m.tmp@coef[c(7:18)]
  points(rep(position$easting[i],2),m.surv[c(1,7)+3],pch=16,col=c("black","red","pink"))
}

par(mfrow=c(1,1))
plot(rep(NA,3),c(1,2,3),ylim=c(0,1),xlim=c(1,20),xlab="Population",ylab="Probability of s4")

for(i in 1:20){
  
  m.tmp = modelList[[orderedIndex[i]]]
  m.surv = m.tmp@coef[c(7:18)]
  points(c(i-.25,i),m.surv[c(1,7)+3],pch=16,col=c("black","red","pink"))
}
abline(v=seq(.5,20.5,by=1),lty='dotted')

legend("bottomleft",c("Round 1","Round 2", "Round 3"),pch=16,col=c("black","red","pink"))



# year 3
plot(rep(NA,1),c(NA),ylim=c(0,1),xlim=c(340,375))

for(i in 1:20){
  
  m.tmp = modelList[[i]]
  m.surv = m.tmp@coef[c(7:18)]
  points(rep(position$easting[i],1),m.surv[c(1)+4],pch=16,col=c("black","red","pink"))
}

plot(rep(NA,1),c(NA),ylim=c(0,1),xlim=c(340,375))

for(i in 1:20){
  
  m.tmp = modelList[[orderedIndex[i]]]
  m.surv = m.tmp@coef[c(7:18)]
  points(rep(position$easting[i],1),m.surv[c(1)+6],pch=16,col=c("black","red","pink"))
}


