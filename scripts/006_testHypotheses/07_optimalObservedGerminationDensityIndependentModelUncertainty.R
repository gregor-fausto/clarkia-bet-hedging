# -------------------------------------------------------------------
# Density-independent model of germination + uncertainty
# -------------------------------------------------------------------

# - Environment ----
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)

# - Libraries ----
library(MCMCvis)
library(dplyr)
library(reshape2)
library(tidyr)
library(HDInterval)
library(bayesplot)

# - Source functions for analysis ----
source("scripts/006_testHypotheses/00_utilityFunctions.R")

# - Function to calculate one step population growth rate ----

fitness <- function(g=g1,s0=s0,s1=s1,s2=s2,s3=s3,rs=rs){
  p1 = g*rs*s0*s1
  p2 = (1-g)*(s2*s3)
  return(as.numeric(p1+p2))
}

# - Site names ----
position<-read.csv(file="~/Dropbox/clarkia-demography-projects/data/siteAbioticData.csv",header=TRUE) %>%
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteNames <- unique(position$site)

# - Create data frame to match sites and years  ----
siteIndex <- data.frame(site=siteNames,siteIndex=1:20)
yearIndex <- data.frame(year=2006:2020,yearIndex=1:15)
index=expand.grid(1:20,1:15)


# - Samples from the posterior distribution ----

s0 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s0-population-level.RDS")
g1 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/g1-population-level.RDS")
s1 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s1-population-level.RDS")
s2 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s2-population-level.RDS")
s3 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s3-population-level.RDS")
perCapitaRS <- readRDS("outputs/005_calculatePopulationModelParameters/04_reproductiveSuccess/reproductiveSuccessWithCorrectionForMissingness-populationYear-mat.RDS")

# - Calculate the mode from the full posterior distribution ----

# get mode from full posterior and calculate var in RS based on modes
g1.hat  <- apply(g1,2,posterior.mode)
s1.hat  <- apply(s1,2,posterior.mode)
s2.hat  <- apply(s2,2,posterior.mode)
s3.hat  <- apply(s3,2,posterior.mode)
s0.hat  <- apply(s0,2,posterior.mode)
perCapitaRS.hat <- lapply(perCapitaRS,apply,2,posterior.mode)

# - Analyze effects of parameter uncertainty ----
# not sure how to intepret this given that 
# there is a LOT of uncertainty in parameters

fr <- function(x,pars,pars_rs){
  x1 = x[1]
  p1 = pars[1]
  p2 = pars[2]
  p3 = pars[3]
  p4 = pars[4]
  p5 = pars_rs
  fit <- fitness(g=x1,s0=p1,s1=p2,s2=p3,s3=p4,rs=p5)
  -gm_mean(fit)
}

fr1 <- function(x,pars,pars_rs){
  if(x<0|x>1){
    return(1000)
  } else {
    fr(x,pars,pars_rs)
  }
}

# code for revision

vec.list = list()

# for each population 
for(k in 1:20){
  
  # - ++get the population index  ----
  pop.index = k
  
  # - ++get full posterior distribution of parameter estimates  ----
  y_t = perCapitaRS[[pop.index]]
  y_t = y_t[,!is.na(y_t[1,])]
  s0_unc = s0[,k]
  s1_unc = s1[,k]
  s2_unc = s2[,k]
  s3_unc = s3[,k]
  
  # draw samples from the full posterior distribution  
  n.reps = 1000
  n.iter = length(s3_unc)
  index.draws=sample(1:n.iter,n.reps)
  
  s0_est = s0_unc[sample(1:n.iter,n.reps)]
  s1_est = s1_unc[sample(1:n.iter,n.reps)]
  s2_est = s2_unc[sample(1:n.iter,n.reps)]
  s3_est = s3_unc[sample(1:n.iter,n.reps)]
  y_t_est = y_t[sample(1:n.iter,n.reps),]
  
  # run a single, long optimization
  y_t.draws = sample(1:(sum(!is.na(y_t[1,]))),100000,replace=TRUE)
  
  vec.vec = c()

  # for each of these samples
  for(j in 1:n.reps){
    
    # draw the estimates for seed survival
    params = c(s0_est[j],s1_est[j],s2_est[j],s3_est[j])
    
    # - ++just draw a 100,000 year sequence with replacement  ----
    y_t.resample =  y_t_est[j,][y_t.draws]
    
    # calculate the optimal germination fraction for this draw
    tmp = optim(c(.5),fr1,pars=params,pars_rs=y_t.resample,method='Brent',lower=0,upper=1)
    vec.vec[j] = tmp$par
  }
  vec.list[[k]] = vec.vec
}



# - Quantify parameter uncertainty  ----

# - ++Calculate the 50% and 95% highest posterior density interval ----
HPDI.50 <- lapply(vec.list,FUN = function(x) hdi(x, .50))
HPDI.95 <- lapply(vec.list,FUN = function(x) hdi(x, .95))
# - ++Calculate the posterior mode ----
mode.optg <- lapply(vec.list, FUN = posterior.mode)

#

# # - Calculate optima with posterior mode  ----

fr <- function(x){
  x1 = x[1]
  fit <- fitness(g=x1,s0.hat[k],s1.hat[k],s2.hat[k],s3.hat[k],rs=y_t.resample)
  -gm_mean(fit)
}

fr1 <- function(x){
  if(x<0|x>1){
    return(1000)
  } else {
    fr(x)
  }
}

iterations = 1
optima.mat.method2 = matrix(NA,nrow=20,ncol=iterations)
for(j in 1:iterations){
  vec=c()
  for(k in 1:20){
    # - ++get the population index  ----
    pop.index = k
    
    # - ++draw the 15 years of reproductive success estimates  ----
    y_t = perCapitaRS.hat[[pop.index]]
    y_t = y_t[!is.na(y_t)]
    
    # - ++draw 1000 samples for reproductive success with replacement  ----
    y_t.resample = sample(y_t,100000,replace=TRUE)
    
    tmp = optim(c(.5),fr1,method='Brent',lower=0,upper=1)
    vec[k]=tmp$par
    
  }
  optima.mat.method2[,j] = vec
}


# # - PLOT  ----

pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12
pt5 = 5/12


pdf("products/figures/optimalGerminationParameterUncertainty-revised.pdf",
    height = 7, width = 7)

index = order(position$easting,decreasing=FALSE)

par(mfrow=c(1,1),mar=c(0,.25,.25,0),
    oma=c(4,4,1,1))

plot(x = NA,
     y = NA,
     xlim=c(1-.5,20+.5),ylim=c(0,1),
     pch=16, cex = 0.5,
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt10)

rect(xleft=c(1,3,5,7,9,11,13,15,17,19)+c(.5),xright=c(1,3,5,7,9,11,13,15,17,19)+1.5,ybottom=-100,ytop=100,col='gray99',border='gray99')

tmp<-do.call(rbind,HPDI.50)
segments(x0=1:20,x1=1:20,
         y0=tmp[index,1],y1=tmp[index,2], lwd=2)
tmp<-do.call(rbind,HPDI.95)
segments(x0=1:20,x1=1:20,
         y0=tmp[index,1],y1=tmp[index,2], lwd=1)

points(1:20,as.vector(do.call(rbind,mode.optg))[index],
       pch=21,cex=1,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 1),lwd=0)
points(1:20,as.vector(do.call(rbind,mode.optg))[index],
       pch=21,col='black',cex=1,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 0.5))

points(1:20,as.vector(optima.mat.method2)[index],
       pch=4,col='red',cex=1,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 0.5))

box()
axis(1,at=1:20,labels=siteNames[index],las=2)
axis(2, seq(0,1,by=.1),
     labels = seq(0,1,by=.1), las = 1, line = 0, hadj= .8,
     col = NA, col.ticks = 1, cex.axis = pt10)
axis(2, seq(.05,1,by=.1),labels=FALSE,tck=-0.0125)
mtext("Predicted optimal germination fraction", 
      side = 2, outer = TRUE, line = 2,cex=pt12)

dev.off()


# # Code for the original submission
# 
# vec.mat = matrix(NA,nrow=100,ncol=50)
# vec.list = vec.list2= list()
# 
# # for each population 
# for(k in 1:20){
#   
#   # - ++get the population index  ----
#   pop.index = k
#   
#   # - ++get full posterior distribution of parameter estimates  ----
#   y_t = perCapitaRS[[pop.index]]
#   y_t = y_t[,!is.na(y_t[1,])]
#   s0_unc = s0[,k]
#   s1_unc = s1[,k]
#   s2_unc = s2[,k]
#   s3_unc = s3[,k]
#   
#   # draw the full posterior distribution  
#   n.reps = 100
#   n.iter = length(s3_unc)
#   index.draws=sample(1:n.iter,n.reps)
#   
#   s0_est = s0_unc[sample(1:n.iter,n.reps)]
#   s1_est = s1_unc[sample(1:n.iter,n.reps)]
#   s2_est = s2_unc[sample(1:n.iter,n.reps)]
#   s3_est = s3_unc[sample(1:n.iter,n.reps)]
#   y_t_est = y_t[sample(1:n.iter,n.reps),]
#   
#   y_t.draws = matrix(NA,nrow=50,ncol=1000)
#   for(tmp in 1:50){
#     y_t.draws[tmp,] = sample(1:15,1000,replace=TRUE)
#   }
#   
#   vec.mat = matrix(NA,nrow=n.reps,ncol=50)
#   
#   vec.mat2 = list()
#   # for each of these samples
#   for(j in 1:n.reps){
#     
#     vec.vec = c()
#     # draw the estimates for seed survival
#     params = c(s0_est[j],s1_est[j],s2_est[j],s3_est[j])
#     
#     # - ++draw 1000 samples for reproductive success with replacement  ----
#     
#     
#     
#     for(h in 1:50){
#       
#       # for each replicate draw from the posterior, draw the same 50 year sequences
#       y_t.resample =  y_t_est[j,][y_t.draws[h,]]# sample(y_t_est[j,],1000,replace=TRUE)
#       
#       # calculate the optimal germination fraction for this draw
#       tmp = optim(c(.5),fr1,pars=params,pars_rs=y_t.resample,method='Brent',lower=0,upper=1)
#       vec.mat[j,h]=tmp$par
#       vec.vec[h] = c(tmp$par)
#     }
#     vec.mat2[[j]] = cbind(optg1=vec.vec,repYear=rep(1:50,each=1),repDraw=j)
#   }
#   vec.list[[k]] = vec.mat
#   vec.list2[[k]] = do.call(rbind,vec.mat2)
# }
# 
# # - Quantify parameter uncertainty: ANOVA  ----
# 
# pdf("products/figures/optimalGerminationParameterUncertainty.pdf",
#     height = 7, width = 7)
# 
# index = order(position$easting,decreasing=FALSE)
# 
# par(mfrow=c(4,5),mar=c(0,.25,.25,0),
#     oma=c(4,4,1,1))
# for(i in 1:20){
#   i.tmp = index[i]
#   df=data.frame(vec.list2[[i.tmp]])
#   names(df) = c("optg1","parmDraw","repDraw")
#   boxplot(optg1~repDraw,data=df,main='',
#           ylim=c(0,1),type='l',xaxt='n',yaxt='n')
#   
#   legend("bottomleft",siteNames[i.tmp],bty='n',cex=.9,inset=c(-.01,0))
#   
#   ifelse(i%in%c(1,6,11,16),axis(2L,las=1),NA)
#   ifelse(i%in%c(16:20),axis(1L),NA)
# }
# 
# mtext("Predicted optimal germination fraction", side = 2, outer = TRUE, line = 2)
# mtext('Index for sample from posterior distribution', side = 1, outer = TRUE, line = 2.3)
# 
# dev.off()
# 
# # - Quantify parameter uncertainty  ----
# # vec.list is length 20 for 20 pops
# # each pop has a 45000x50 matrix
# # where each draw from the posterior was used to calculate the optimal germination fraction 50 times
# 
# # calculate the optimal germination fraction for each row in the posterior
# # by taking row means
# optimalGerminationFraction <- lapply(vec.list,apply,1,mean)
# 
# # calculate the uncertainty in optimal germination fraction by taking posterior mode + HDI
# 
# # - +Summary statistics: g1 ----
# # - ++Calculate the 68% highest posterior density interval ----
# HPDI.g1.68 <- apply(as.matrix(g1),2,FUN = function(x) quantile(x,c( .5-.34,.5+.34)))
# #HPDI.g1.95 <- apply(as.matrix(g1),2,FUN = function(x) hdi(x, .95))
# # - ++Calculate the posterior mode ----
# mode.g1 <- apply(g1,2, FUN = posterior.mode)
# 
# # - ++Construct data frame ----
# # use highest posterior density interval and mode
# g1PosteriorSummary <- data.frame(cbind(t(HPDI.g1.68),mode.g1))
# names(g1PosteriorSummary) <- c("lo.g1","hi.g1","mode.g1")
# 
# # - +Summary statistics: optimal_g1 ----
# # - ++Calculate the 68% highest posterior density interval ----
# HPDI.optG1.68 <- lapply(optimalGerminationFraction,FUN = function(x) quantile(x,c( .5-.34,.5+.34)))
# HPDI.optG1.68 <- t(do.call(rbind,HPDI.optG1.68))
# #HPDI.optG1.95 <- apply(optimalGerminationFraction,2,FUN = function(x) hdi(x, .95))
# # - ++Calculate the posterior mode ----
# mode.optG1 <- unlist(lapply(optimalGerminationFraction, FUN = posterior.mode))
# 
# # - ++Construct data frame ----
# # use highest posterior density interval and mode
# optimalG1PosteriorSummary <- data.frame(cbind(t(HPDI.optG1.68),mode.optG1))
# names(optimalG1PosteriorSummary) <- c("lo.g1","hi.g1","mode.g1")
# 
# 
# # - PLOT  ----
# 
# tiff(filename=paste0("products/figures/optimalGerminationFractionPlusUncertainty.tif"),
#      height=3.75,width=3.5,units="in",res=300,compression="lzw",pointsize=12)
# 
# par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(2,2.2,.7,1)+.1,mgp=c(3,.45,0))
# 
# plot(NA,NA,xlim=c(.4,1),ylim=c(0,.5),
#      xlab = "",
#      ylab = "",bty='n',
#      xaxt= "n", yaxt="n",)
# 
# 
# segments(x0=optimalG1PosteriorSummary[,1],x1=optimalG1PosteriorSummary[,2],
#          y0=g1PosteriorSummary[,3],pch=21,col='black',bg='white',cex=2.5,xlim=c(0,1),ylim=c(0,1))
# segments(y0=g1PosteriorSummary[,1],y1=g1PosteriorSummary[,2],
#          x0=optimalG1PosteriorSummary[,3],pch=21,col='black',bg='white',cex=2.5,xlim=c(0,1),ylim=c(0,1))
# 
# 
# points(x=optimalG1PosteriorSummary[,3],y=g1PosteriorSummary[,3],
#        pch=21,cex=1.5,
#        bg=rgb(red = 1, green = 1, blue = 1, alpha = 1),lwd=0)
# points(x=optimalG1PosteriorSummary[,3],y=g1PosteriorSummary[,3],
#        lwd = 0.5,
#        pch=21,col='black',cex=1.5,
#        bg=rgb(red = 1, green = 1, blue = 1, alpha = 0.5))
# 
# box()
# 
# axis(1, seq(0,1,by=.1), padj = -.5,
#      labels = seq(0,1,by=.1), line = 0,
#      col = NA, col.ticks = 1, cex.axis = pt8)
# axis(1, seq(.25,1,by=.1),labels=FALSE)
# axis(2, seq(0,1,by=.1),
#      labels = seq(0,1,by=.1), las = 1, line = 0, hadj= 1.2,
#      col = NA, col.ticks = 1, cex.axis = pt8)
# axis(2, seq(.05,1,by=.1),labels=FALSE)
# 
# mtext("Observed germination fraction",
#       side=2,line=1.5,adj=.5,col='black',cex=pt10)
# mtext("Predicted germination fraction",
#       side=1,line=1,adj=.5,col='black',cex=pt10)
# 
# dev.off()


