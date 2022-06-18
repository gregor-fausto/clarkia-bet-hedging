# -------------------------------------------------------------------
# Density-independent model of germination + uncertainty
# -------------------------------------------------------------------

# - Environment ----
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)

# - Function to calculate one step population growth rate ----

fitness <- function(g=g1,s0=s0,s1=s1,s2=s2,s3=s3,rs=rs){
  p1 = g*rs*s0*s1
  p2 = (1-g)*(s2*s3)
  return(as.numeric(p1+p2))
}

# - Site names ----
position<-read.csv(file="data/siteAbioticData.csv",header=TRUE) %>%
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteNames <- unique(position$site)

# - Create data frame to match sites and years  ----
siteIndex <- data.frame(site=siteNames,siteIndex=1:20)
yearIndex <- data.frame(year=2006:2020,yearIndex=1:15)
index=expand.grid(1:20,1:15)

# - Libraries ----
library(MCMCvis)
library(dplyr)
library(reshape2)
library(tidyr)
library(HDInterval)
library(bayesplot)

# - Samples from the posterior distribution ----

s0 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s0-population-level.RDS")
g1 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/g1-population-level.RDS")
s1 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s1-population-level.RDS")
s2 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s2-population-level.RDS")
s3 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/s3-population-level.RDS")
rs <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/reproductiveSuccess-population-year-level-mat.RDS")

# - Calculate the mode from the full posterior distribution ----

# get mode from full posterior and calculate var in RS based on modes
g1.hat  <- apply(g1,2,posterior.mode)
s1.hat  <- apply(s1,2,posterior.mode)
s2.hat  <- apply(s2,2,posterior.mode)
s3.hat  <- apply(s3,2,posterior.mode)
s0.hat  <- apply(s0,2,posterior.mode)
rs.hat <- lapply(rs,apply,2,posterior.mode)



# - Compare 2 methods for calculating the optimal value of germination ----

# - + Grid search----
# for each population k 
# get the mode of per-capita reproductive success in each year
# sample a sequence of n=1000 years to use in the analysis for each value of g
# for each g, calculate the growth rate for n=1000 years

g = seq(0,1,by=.001)
g.seq = g
fit = c()
g.mat = matrix(NA,nrow=1000,ncol=length(g))
g.sites = list()
yt = c()
n.iter=dim(rs)[1]


# - + One dimensional optimization----
# use the functions below to calculate the optimal germination fraction

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

# - + Compare grid search vs. one dimensional optimization ----

vec=c()
for(k in 1:20){
  
  # START WITH GRID SEARCH
  # - ++get the population index  ----
  pop.index = k
  
  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = rs.hat[[pop.index]]
  y_t = y_t[!is.na(y_t)]
  
  # - ++draw 1000 samples for reproductive success with replacement  ----
  y_t.resample = sample(y_t,1000,replace=TRUE)
  
  for( i in 1:length(g)){
    fit<-fitness(g=g[i],s0.hat[k],s1.hat[k],s2.hat[k],s3.hat[k],rs=y_t.resample)
    g.mat[,i] <- fit
  }
  g.sites[[k]] <- g.mat
  
  # ONE DIMENSIONAL OPTIMIZATION
  tmp = optim(c(.5),fr1,method='Brent',lower=0,upper=1)
  vec[k]=tmp$par
  
}

# - ++calculate optimal germination fraction based on grid search  ----

gmean <- function(x){apply(x,2,gm_mean)}
maxfun <- function(x){g[which(x %in% max(x))]}

# construct a list of the population growth rate at each value of g
site.optima<-lapply(g.sites,gmean)
# unlist and calculate the value of g that maximizes population growth
optima<-unlist(lapply(site.optima,maxfun))

# - ++calculate optimal germination fraction based on 1D optimization  ----
vec

par(mfrow=c(1,1))
# compare optima with optim vs grid search
plot(optima,vec,ylab="optima grid",xlab="optima optim",xlim=c(0,1.1),ylim=c(0,1.1));
abline(a=0,b=1)
text(optima+.05,vec,1:20)

# - Calculate optimal g ----
# repeat the optimization 1000 times and calculate
# the mean and variance of optimal germination fractions
# for each repetition, we resample the sequence y_t

vec=c()
iterations = 1000
optima.mat = matrix(NA,nrow=20,ncol=iterations)
for(j in 1:iterations){
  for(k in 1:20){
    # - ++get the population index  ----
    pop.index = k
    
    # - ++draw the 15 years of reproductive success estimates  ----
    y_t = rs.hat[[pop.index]]
    y_t = y_t[!is.na(y_t)]
    
    # - ++draw 1000 samples for reproductive success with replacement  ----
    y_t.resample = sample(y_t,1000,replace=TRUE)
    
    tmp = optim(c(.5),fr1,method='Brent',lower=0,upper=1)
    vec[k]=tmp$par
    
  }
  optima.mat[,j] = vec
}

# - + summarize optimal g ----

# calculate mean and variance of optimal germination fractions
optima.mean <-apply(optima.mat,1,mean)
optima.var <- apply(optima.mat,1,var)

plot(x=optima.mean,y=g1.hat,pch=21,col='black',bg='white',cex=2.5,xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1,lty='dotted')

# mean and variance of optimal germination fractions are negatively correlated
plot(x=optima.mean,y=optima.var)
cor(optima.mean,optima.var)

# calculate confidence intervals on estimate
# these represent the uncertainty due to resampling the modes for reproductive success
optima.est = apply(optima.mat,1,quantile,c(.025,.5,.975))
plot(NA,pch=21,col='black',bg='white',cex=2.5,xlim=c(0,1),ylim=c(0,1))
segments(x0=optima.est[1,],x1=optima.est[3,],y0=g1.hat,pch=21,col='black',bg='white',cex=2.5,xlim=c(0,1),ylim=c(0,1))
points(x=optima.est[2,],y=g1.hat,pch=21,bg='white')
abline(a=0,b=1,lty='dotted')




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



vec.mat = matrix(NA,nrow=500,ncol=50)
vec.list = list()

# for each population 
for(k in 1:20){
  
  # - ++get the population index  ----
  pop.index = k
  
  # - ++get full posterior distribution of parameter estimates  ----
  y_t = rs[[pop.index]]
  y_t = y_t[,!is.na(y_t[1,])]
  s0_unc = s0[,k]
  s1_unc = s1[,k]
  s2_unc = s2[,k]
  s3_unc = s3[,k]
  
  # draw 500 samples from the posterior distribution  
  n.reps = 500
  n.iter = length(s3_unc)
  index.draws=sample(1:n.iter,n.reps)
  
  s0_est = s0_unc[sample(1:n.iter,n.reps)]
  s1_est = s1_unc[sample(1:n.iter,n.reps)]
  s2_est = s2_unc[sample(1:n.iter,n.reps)]
  s3_est = s3_unc[sample(1:n.iter,n.reps)]
  y_t_est = y_t[sample(1:n.iter,n.reps),]
  
  vec.mat = matrix(NA,nrow=500,ncol=50)
  # for each of these samples
  for(j in 1:n.reps){
    
    # draw the estimates for seed survival
    params = c(s0_est[j],s1_est[j],s2_est[j],s3_est[j])
    
    for(h in 1:50){
      # - ++draw 1000 samples for reproductive success with replacement  ----
      y_t.resample = sample(y_t_est[j,],1000,replace=TRUE)
      
      # calculate the optimal germination fraction for this draw
      tmp = optim(c(.5),fr1,pars=params,pars_rs=y_t.resample,method='Brent',lower=0,upper=1)
      vec.mat[j,h]=tmp$par
    }
  }
  vec.list[[k]] = vec.mat
}

par(mfrow=c(2,5))
for(i in 1:10){
hist(apply(vec.list[[i]],1,mean))
}

