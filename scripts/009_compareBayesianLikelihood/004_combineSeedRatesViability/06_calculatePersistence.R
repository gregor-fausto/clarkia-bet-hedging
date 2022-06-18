rm(list=(ls(all=TRUE))) # if using in source(script)

# - Libraries ----
library(tidyverse)
library(parallel)
library(stringr)
library(bbmle)
library(MCMCvis)

# - Read in estimates ----

likelihoodEstimates <- readRDS(paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","seedRateEstimatesFrequentistCombinedPopulationLevel.rds"))
bayesianEstimates <- readRDS(paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","seedRateEstimatesBayesianCombinedPopulationLevel.rds"))

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

# - Frequentist estimates ----

germinationLikelihood <- likelihoodEstimates[[1]]
survivalLikelihood <- likelihoodEstimates[[2]]

# - Calculate probability of persistence ----

# - +persistence & viability ----

viableAfterOneYear.freq <-survivalLikelihood[,1] * (1-germinationLikelihood[,1]) * survivalLikelihood[,2]  
viableAfterOneYear.bayes <- s1[,1:20] * (1-g1[,1:20]) * s2[,1:20]  

viableAfterOneYear.bayes.summary <- apply(viableAfterOneYear.bayes,2,quantile,c(.025,.25,.5,.75,.975))



# - +persistence & viability ----

viableAfterTwoYears.freq <-survivalLikelihood[,1] * (1-germinationLikelihood[,1]) * survivalLikelihood[,2]  * survivalLikelihood[,3] * (1-germinationLikelihood[2]) * survivalLikelihood[,4]  
viableAfterTwoYears.bayes <- s1[,1:20] * (1-g1[,1:20]) * s2[,1:20]  * s3[,1:20] * (1-g2[,1:20]) * s4[,1:20]  

viableAfterTwoYears.bayes.summary <- apply(viableAfterTwoYears.bayes,2,quantile,c(.025,.25,.5,.75,.975))


# - +persistence & viability ----

viableAfterThreeYears.freq <-survivalLikelihood[,1] * (1-germinationLikelihood[,1]) * survivalLikelihood[,2]  * survivalLikelihood[,3] * (1-germinationLikelihood[2]) * survivalLikelihood[,4]  * survivalLikelihood[,5] * (1-germinationLikelihood[3]) * survivalLikelihood[,6]
viableAfterThreeYears.bayes <- s1[,1:20] * (1-g1[,1:20]) * s2[,1:20]  * s3[,1:20] * (1-g2[,1:20]) * s4[,1:20]  * s5[,1:20] * (1-g3[,1:20]) * s6[,1:20]

viableAfterThreeYears.bayes.summary <- apply(viableAfterThreeYears.bayes,2,quantile,c(.025,.25,.5,.75,.975))



# - Calculate probability of persistence ----

par(mfrow=c(1,3),mar=c(4,4,4,4))

plot(viableAfterOneYear.freq,viableAfterOneYear.bayes.summary[3,],xlim=c(0,1),ylim=c(0,1),
     main="Persistent & viable after 1 year",xlab="Frequentist estimate",ylab="Bayesian estimate")
abline(a=0,b=1)
segments(x0=viableAfterOneYear.freq,y0=viableAfterOneYear.bayes.summary[1,],y1=viableAfterOneYear.bayes.summary[5,])

plot(viableAfterTwoYears.freq,viableAfterTwoYears.bayes.summary[3,],xlim=c(0,1),ylim=c(0,1),
     main="Persistent & viable after 2 years",xlab="Frequentist estimate",ylab="Bayesian estimate")
abline(a=0,b=1)
segments(x0=viableAfterTwoYears.freq,y0=viableAfterTwoYears.bayes.summary[1,],y1=viableAfterTwoYears.bayes.summary[5,])

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

par(mfrow=c(1,1),mar=c(4,4,4,4))

orderedIndex = order(position$easting,decreasing=FALSE)
siteNames = unique(position$site)

plot(position$easting, viableAfterThreeYears.freq,ylim=c(0,.3),col='black',pch=16,type='n',
     xlab="Easting (km)", ylab="Probability seeds are intact and viable in October at end of year 3")
points(position$easting, viableAfterThreeYears.bayes.summary[3,],col='red',pch=16)
segments(x0=position$easting,y0=viableAfterThreeYears.bayes.summary[1,],y1=viableAfterThreeYears.bayes.summary[5,],col='red')
points(position$easting, viableAfterThreeYears.freq,col='black',pch=16)

legend("topleft",pch=16,legend=c("MLE estimates","Bayes estimates"),col=c("black",'red'),bty='n')


# summarize

position = position %>% dplyr::filter(intenseDemography==1) %>% unique
siteNames=droplevels(position$site)

freq.summary = c(viableAfterOneYear.freq,viableAfterTwoYears.freq,viableAfterThreeYears.freq)
bayes.summary = c(viableAfterOneYear.bayes.summary[3,],viableAfterTwoYears.bayes.summary[3,],viableAfterThreeYears.bayes.summary[3,])
bayes.lo = c(viableAfterOneYear.bayes.summary[1,],viableAfterTwoYears.bayes.summary[1,],viableAfterThreeYears.bayes.summary[1,])
bayes.hi = c(viableAfterOneYear.bayes.summary[5,],viableAfterTwoYears.bayes.summary[5,],viableAfterThreeYears.bayes.summary[5,])

summaryDf <- data.frame(site=rep(siteNames,3),persistentAndViableAfterNYears=rep(1:3,each=20),
           maximumLikelihoodEstimate=freq.summary,bayesEstimatePosteriorMedian=bayes.summary,
           bayesEstimate95CredibleIntervalLower=bayes.lo,
           bayesEstimate95CredibleIntervalHigh=bayes.hi)

write.csv(summaryDf,file="~/Dropbox/chapter-3/analysis/outputs/seedPersistenceViability.csv")

subsetDF <- summaryDf %>%
  dplyr::filter(site%in%c("KYE","MC","S22")) 

pdf("~/Dropbox/chapter-3/analysis/outputs/seedPersistenceViability.pdf",width=5,height=5)

par(mar=c(4,5,4,1),mfrow=c(1,1))

plot(subsetDF$persistentAndViableAfterNYears,subsetDF$bayesEstimatePosteriorMedian,type='n',ylim=c(0,.7),
     xlab="Years after seed bags were buried",ylab="Probability seeds remain intact & viable in seed bags")
for(i in 1:3){
  tmp.site <- unique(subsetDF$site)[i]
  tmp.df <- subsetDF[subsetDF$site==tmp.site,]
  colors = c("black","coral","cornflowerblue")
  offset = c(-.05,0,.05)
  points(tmp.df$persistentAndViableAfterNYears+offset[i],tmp.df$bayesEstimatePosteriorMedian,type='b', 
         col=colors[i],pch=16)
  segments(x0=tmp.df$persistentAndViableAfterNYears+offset[i],
           y0=tmp.df$bayesEstimate95CredibleIntervalLower,
           y1=tmp.df$bayesEstimate95CredibleIntervalHigh,
           type='b', 
         col=colors[i])
  
  points(tmp.df$persistentAndViableAfterNYears+offset[i],tmp.df$maximumLikelihoodEstimate, 
         col=colors[i],pch=2)
}

legend("topright",c("KYE","MC","S22"),pch=16,col=colors,bty='n')
dev.off()
