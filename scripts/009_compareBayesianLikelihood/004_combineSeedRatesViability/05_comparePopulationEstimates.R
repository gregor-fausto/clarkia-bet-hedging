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

# - Compare estimates ----

germinationBayesian <- cbind(g1.summary[3,1:20],g2.summary[3,1:20],g3.summary[3,1:20])

germinationBayesian.mode <- cbind(g1.mode[1:20],g2.mode[1:20],g3.mode[1:20])

survivalBayesian <- cbind(s1.summary[3,1:20],s2.summary[3,1:20],s3.summary[3,1:20],
                          s4.summary[3,1:20],s5.summary[3,1:20],s6.summary[3,1:20])

survivalBayesian.mode <- cbind(s1.mode[1:20],s2.mode[1:20],s3.mode[1:20],
                               s4.mode[1:20],s5.mode[1:20],s6.mode[1:20])

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

par(mfrow=c(1,3),mar=c(3,3,2,2))

plot(germinationLikelihood[,1],germinationBayesian[,1],xlim=c(0,.4),ylim=c(0,.4),
     main="g1")
abline(a=0,b=1)
segments(x0=germinationLikelihood[,1],y0=g1.summary[1,1:20],y1=g1.summary[5,1:20])

plot(germinationLikelihood[,2],germinationBayesian[,2],xlim=c(0,.6),ylim=c(0,.6),
     main="g2")
abline(a=0,b=1)
segments(x0=germinationLikelihood[,2],y0=g2.summary[1,1:20],y1=g2.summary[5,1:20])

plot(germinationLikelihood[,3],germinationBayesian[,3],xlim=c(0,.8),ylim=c(0,.8),
     main="g3")
abline(a=0,b=1)
segments(x0=germinationLikelihood[,3],y0=g3.summary[1,1:20],y1=g3.summary[5,1:20])
