rm(list=(ls(all=TRUE))) # if using in source(script)

# - Libraries ----
library(tidyverse)
library(parallel)
library(stringr)
library(bbmle)

# - Read in estimates ----

likelihoodEstimates <- readRDS(paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","seedRateEstimatesFrequentistCombined.rds"))
bayesianEstimates <- readRDS(paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","seedRateEstimatesBayesianCombined.rds"))

# - Functions ----

posterior.mode = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]
  return(modeParam)
}


# - Bayesian estimates ----

germinationList <- bayesianEstimates[[1]]
survivalList <- bayesianEstimates[[2]]

# - + Compute bayesian medians ----

# - ++ g1 ----

g1 <- germinationList[[1]]
g1.summary <- apply(g1,2,quantile,c(0.025,.25,.5,.75,.975))
g1.mode <- apply(g1,2,posterior.mode)

# - ++ g2 ----

g2 <- germinationList[[2]]
g2.summary <- apply(g2,2,quantile,c(0.025,.25,.5,.75,.975))
g2.mode <- apply(g2,2,posterior.mode)

# - ++ g3 ----

g3 <- germinationList[[3]]
g3.summary <- apply(g3,2,quantile,c(0.025,.25,.5,.75,.975))
g3.mode <- apply(g3,2,posterior.mode)

# - ++ s1 ----

s1 <- survivalList[[1]]
s1.summary <- apply(s1,2,quantile,c(0.025,.25,.5,.75,.975))
s1.mode <- apply(s1,2,posterior.mode)

# - ++ s2 ----

s2 <- survivalList[[2]]
s2.summary <- apply(s2,2,quantile,c(0.025,.25,.5,.75,.975))
s2.mode <- apply(s2,2,posterior.mode)

# - ++ s3 ----

s3 <- survivalList[[3]]
s3.summary <- apply(s3,2,quantile,c(0.025,.25,.5,.75,.975))
s3.mode <- apply(s3,2,posterior.mode)

# - ++ s4 ----

s4 <- survivalList[[4]]
s4.summary <- apply(s4,2,quantile,c(0.025,.25,.5,.75,.975))
s4.mode <- apply(s4,2,posterior.mode)

# - ++ s5 ----

s5 <- survivalList[[5]]
s5.summary <- apply(s5,2,quantile,c(0.025,.25,.5,.75,.975))
s5.mode <- apply(s5,2,posterior.mode)

# - ++ s6 ----

s6 <- survivalList[[6]]
s6.summary <- apply(s6,2,quantile,c(0.025,.25,.5,.75,.975))
s6.mode <- apply(s6,2,posterior.mode)

# - + Summarize estimates by medians ----

germinationBayesian <- cbind(g1.summary[3,1:20],g2.summary[3,1:20],g3.summary[3,],
                             g1.summary[3,21:40],g2.summary[3,21:40],
                             g1.summary[3,41:60])

germinationBayesian.mode <- cbind(g1.mode[1:20],g2.mode[1:20],g3.mode[],
                             g1.mode[21:40],g2.mode[21:40],
                             g1.mode[41:60])

survivalBayesian <- cbind(s1.summary[3,1:20],s2.summary[3,1:20],s3.summary[3,1:20],
                          s4.summary[3,1:20],s5.summary[3,1:20],s6.summary[3,1:20],
                          s1.summary[3,21:40],s2.summary[3,21:40],s3.summary[3,21:40],
                          s4.summary[3,21:40],
                          s1.summary[3,41:60],s2.summary[3,41:60])

survivalBayesian.mode <- cbind(s1.mode[1:20],s2.mode[1:20],s3.mode[1:20],
                          s4.mode[1:20],s5.mode[1:20],s6.mode[1:20],
                          s1.mode[21:40],s2.mode[21:40],s3.mode[21:40],
                          s4.mode[21:40],
                          s1.mode[41:60],s2.mode[41:60])

# - Frequentist estimates ----

germinationLikelihood <- likelihoodEstimates[[1]]
survivalLikelihood <- likelihoodEstimates[[2]]

# - Compare estimates ----

par(mfrow=c(1,2))
plot(germinationLikelihood,germinationBayesian,xlim=c(0,.6),ylim=c(0,.6))
abline(a=0,b=1)

plot(germinationLikelihood,germinationBayesian.mode,xlim=c(0,.6),ylim=c(0,.6))
abline(a=0,b=1)

plot(survivalLikelihood,survivalBayesian,xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)

plot(survivalLikelihood,survivalBayesian.mode,xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)

# - + Compare germination ----

par(mfrow=c(3,3),mar=c(3,3,2,2))

plot(germinationLikelihood[,1],germinationBayesian[,1],xlim=c(0,.4),ylim=c(0,.4),
     main="Round 1: g1")
abline(a=0,b=1)
segments(x0=germinationLikelihood[,1],y0=g1.summary[1,1:20],y1=g1.summary[5,1:20])

plot(germinationLikelihood[,2],germinationBayesian[,2],xlim=c(0,.6),ylim=c(0,.6),
     main="Round 1: g2")
abline(a=0,b=1)
segments(x0=germinationLikelihood[,2],y0=g2.summary[1,1:20],y1=g2.summary[5,1:20])

plot(germinationLikelihood[,3],germinationBayesian[,3],xlim=c(0,.8),ylim=c(0,.8),
     main="Round 1: g3")
abline(a=0,b=1)
segments(x0=germinationLikelihood[,3],y0=g3.summary[1,1:20],y1=g3.summary[5,1:20])

plot(germinationLikelihood[,4],germinationBayesian[,4],xlim=c(0,.3),ylim=c(0,.3),
     main="Round 2: g1")
abline(a=0,b=1)
segments(x0=germinationLikelihood[,4],y0=g1.summary[1,21:40],y1=g1.summary[5,21:40])

plot(germinationLikelihood[,5],germinationBayesian[,5],xlim=c(0,.45),ylim=c(0,.45),
     main="Round 2: g2")
abline(a=0,b=1)
segments(x0=germinationLikelihood[,5],y0=g2.summary[1,21:40],y1=g2.summary[5,21:40])

plot(NA,xlim=c(0,.6),ylim=c(0,.6),type='n',bty='n',xaxt='n',yaxt='n',
     main="")

plot(germinationLikelihood[,6],germinationBayesian[,6],xlim=c(0,.7),ylim=c(0,.7),
     main="Round 3: g1")
abline(a=0,b=1)
segments(x0=germinationLikelihood[,6],y0=g1.summary[1,41:60],y1=g1.summary[5,41:60])

plot(NA,xlim=c(0,.6),ylim=c(0,.6),type='n',bty='n',xaxt='n',yaxt='n',
     main="")

plot(NA,xlim=c(0,.6),ylim=c(0,.6),type='n',bty='n',xaxt='n',yaxt='n',
     main="")

# - Calculate probability of persistence ----

# - +persistence & viability: round 1 ----

viableAfterThreeYears.freq <-survivalLikelihood[,1] * (1-germinationLikelihood[,1]) * survivalLikelihood[,2]  * survivalLikelihood[,3] * (1-germinationLikelihood[2]) * survivalLikelihood[,4]  * survivalLikelihood[,5] * (1-germinationLikelihood[3]) * survivalLikelihood[,6]
viableAfterThreeYears.bayes <- s1[,1:20] * (1-g1[,1:20]) * s2[,1:20]  * s3[,1:20] * (1-g2[,1:20]) * s4[,1:20]  * s5[,1:20] * (1-g3[,1:20]) * s6[,1:20]
viableAfterThreeYears.bayes.summary <- apply(viableAfterThreeYears.bayes,2,quantile,c(.025,.25,.5,.75,.975))

par(mfrow=c(1,1),mar=c(4,4,4,4))

plot(viableAfterThreeYears.freq,viableAfterThreeYears.bayes.summary[3,],xlim=c(0,.3),ylim=c(0,.3),
     main="Persistent & viable after 3 years",xlab="Frequentist estimate",ylab="Bayesian estimate")
abline(a=0,b=1)
segments(x0=viableAfterThreeYears.freq,y0=viableAfterThreeYears.bayes.summary[1,],y1=viableAfterThreeYears.bayes.summary[5,])

# - Read in additional data ----

climate <- readRDS("~/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData.RDS")
position<-climate %>% 
  dplyr::ungroup() %>%
  dplyr::select(site,easting,intenseDemography) %>%
  dplyr::mutate(easting=easting/1000)

position = position %>% dplyr::filter(intenseDemography==1) %>% unique
siteNames=droplevels(position$site)
siteNames = siteNames[order(position$easting)]

plot(NA,NA,type='n',xlim=c(0,.35),ylim=c(0,20),
     axes=FALSE,frame=FALSE,
     xlab="",ylab="")
y.pt = 20:1
for(i in 1:20){
  index=order(position$easting)[i]
  tmp<-viableAfterThreeYears.bayes.summary[,index]
  segments(x0=tmp[1],x1=tmp[5],y0=y.pt[i]-.2)
  segments(x0=tmp[2],x1=tmp[4],y0=y.pt[i]-.2,lwd=3)
  points(x=tmp[3],y=y.pt[i]-.2,pch=21,bg='white')
  tmp2 <- viableAfterThreeYears.freq[index]
  points(x=tmp2,y=y.pt[i]-.2,pch=16,col='purple')
  }
axis(1,  seq(0,1,by=.1), col.ticks = 1)
axis(2, (1:20),
     labels = rev(siteNames), las = 1, 
     col = NA, col.ticks = 1, cex.axis = 1)

mtext("Probability seed is persistent & viable after 3 years",side=1,line=2)



# - Plot against space ----


orderedIndex = order(position$easting,decreasing=FALSE)
siteNames = unique(position$site)

plot(position$easting, viableAfterThreeYears.freq,ylim=c(0,.3),col='black',pch=16,type='n',
     xlab="Easting (km)", ylab="Probability seeds are intact and viable in October at end of year 3")
points(position$easting, viableAfterThreeYears.bayes.summary[2,],col='red',pch=16)
segments(x0=position$easting,y0=viableAfterThreeYears.bayes.summary[1,],y1=viableAfterThreeYears.bayes.summary[3,],col='red')
points(position$easting, viableAfterThreeYears.freq,col='black',pch=16)


legend("topleft",pch=16,legend=c("MLE estimates","Bayes estimates"),col=c("black",'red'),bty='n')


seedBagsDataIntactSeeds = readRDS(paste0(dataDirectory,"seedBagsDataIntactSeeds.RDS"))

seedBagsDataIntactSeeds %>%
  dplyr::filter(round==1&age==3) %>%
  View()
