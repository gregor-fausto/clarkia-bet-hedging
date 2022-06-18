# -------------------------------------------------------------------
# Density-independent model of germination
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(HDInterval)
library(bayesplot)

# -------------------------------------------------------------------
# Functions for use when analyzing data
# -------------------------------------------------------------------
temporal_variance <- function(x,fun=var){
  apply(x,1,fun)
}

cols_fun <- function(x,fun=var){
  apply(x,2,fun)
}

# geometric var
gsd <- function(x){
  y <- exp(sd(log(x)))
  return(y)
}


f<-function(x="parm",chain){
  chain<-MCMCchains(chain=belowground,params = x)
  p<-boot::inv.logit(chain)
  BCI <- t(apply(p,2,FUN = function(x) quantile(x, c(.025, .5, .975))))
  return(BCI)
}

# geometric mean from 
# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x), na.rm=na.rm) / length(x))
}

# geometric mean from +.5 to 0
# gm_mean = function(x, na.rm=TRUE){
#   exp(sum(log(x+.5)) / length(x))
# }


# -------------------------------------------------------------------

# read in samples from posterior distributions
s0 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/s0-pop.RDS")
s1 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/s1-pop.RDS")
g1 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/g1-pop.RDS")
s2 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/s2-pop.RDS")
s3 <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/s3-pop.RDS")
rs <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/006_derivedQuantities/estimates/rsPosterior.RDS")

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) 
siteNames = unique(position$site)

siteIndex <- data.frame(site=siteNames,siteIndex=1:20)
yearIndex <- data.frame(year=2006:2020,yearIndex=1:15)
index=expand.grid(1:20,1:15)

# -------------------------------------------------------------------

fitness <- function(g=g1,s0=s0,s1=s1,s2=s2,s3=s3,rs=rs){
  p1 = g*rs*s1*s0
  p2 = (1-g)*(s2*s3)
  return(as.numeric(p1+p2))
}

posterior.median <- function(x){return(apply(x,2,median))}

posterior.mode = function(x){
  x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]
  return(modeParam)
}

modeEst <- function(x){return(apply(x,2,posterior.mode))}

# use mode as Bayesian estimator
s0.hat = modeEst(s0)
g1.hat = modeEst(g1)
s1.hat = modeEst(s1)
s2.hat = modeEst(s2)
s3.hat = modeEst(s3)
rs.hat = modeEst(rs)

# -------------------------------------------------------------------
# calculate autocorrelation of per-capita reproductive success (mode)
# no significant autocorrelation, include in appendix
# -------------------------------------------------------------------

allPopulations = order((position %>% dplyr::select(site,easting))$easting)

par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
for( k in allPopulations){
  pop.index=index[,1]==k
  rs.tmp = rs.hat[pop.index]
  yt = rs.tmp
  acf(yt)
  text(5,.9,(position %>% dplyr::select(site,easting))$site[k])
}

# -------------------------------------------------------------------
# calculate harmonic mean for each population
# -------------------------------------------------------------------

f = function(x){1/mean(1/x)}

pop.harmonicMean =c ()
for( k in 1:20){
  pop.index=index[,1]==k
  rs.tmp = c(rs.hat[pop.index])
  hm = f(rs.tmp)
  pop.harmonicMean[k] = hm
}

par(mfrow=c(1,1))
plot(pop.harmonicMean,s2.hat*s3.hat,type='n');
abline(a=0,b=1);
text(pop.harmonicMean,s2.hat*s3.hat,siteNames)

# -------------------------------------------------------------------
# Simulation to calculate optimal value of germination
# -------------------------------------------------------------------

g = seq(0,1,by=.001)
g.seq = g
fit = c()
g.mat = matrix(NA,nrow=10000,ncol=length(g))
g.sites = list()
yt = c()
n.iter=dim(rs)[1]

# -------------------------------------------------------------------
# for each population k (use mode)
# get the median of per-capita reproductive success in each year
# sample a sequence of 1000 years to use in the analysis for each value of g
# for each g, calculate the growth rate for 1000 years
# -------------------------------------------------------------------


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

vec=c()

for(k in 1:20){
  # - ++get the population index  ----
  pop.index=index[,1]==k
  
  # - ++draw the 15 years of reproductive success estimates  ----
  y_t = rs.hat[pop.index]
  
  # - ++draw 10000 samples for reproductive success with replacement  ----
  y_t.resample = sample(y_t,10000,replace=TRUE)
  
  for( i in 1:length(g)){
    fit<-fitness(g=g[i],s0.hat[k],s1.hat[k],s2.hat[k],s3.hat[k],rs=y_t.resample)
    #logfit<-log(fit)
    g.mat[,i] <- fit
    #g.mat[,i] <- logfit
  }
g.sites[[k]] <- g.mat

tmp = optim(c(.5),fr1,method='Brent',lower=0,upper=1)
vec[k]=tmp$par

}
  gmean <- function(x){apply(x,2,gm_mean)}
  site.optima<-lapply(g.sites,gmean)
  
  maxfun <- function(x){g[which(x %in% max(x))]}
  optima<-unlist(lapply(site.optima,maxfun))
  optima
  
  # for(k in 1:20){
  #   # - ++get the population index  ----
  #   pop.index=index[,1]==k
  #   
  #   # - ++draw the 15 years of reproductive success estimates  ----
  #   y_t = rs.hat[pop.index]
  #   
  #   # - ++draw 10000 samples for reproductive success with replacement  ----
  #   y_t.resample = sample(y_t,10000,replace=TRUE)
  #   
  #   tmp = optim(c(.5),fr1,method='Brent',lower=0,upper=1)
  #   vec[k]=tmp$par
  # 
  # }
  
  dev.off()
# compare optima with optim vs grid search
    plot(optima,vec,ylab="optima grid",xlab="optima optim");abline(a=0,b=1)
    text(optima+.05,vec,1:20)
    
    # now look at effect of uncertainty on one site
    
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
    
    
    
    vec=c()
    vec.list = list()
      
    for(k in 1:20){
      
      # perturbation
     
      
      # - ++get the population index  ----
      pop.index=index[,1]==k
   #   y_t_unc = rs[,pop.index]
    
      # - ++draw the 15 years of reproductive success estimates  ----
      y_t = rs.hat[pop.index]
    
      # - ++draw 10000 samples for reproductive success with replacement  ----
    #  y_t.resample = sample(y_t,10000,replace=TRUE)
    
     # s0_unc = s0.hat[k]
     # s1_unc = s1.hat[k]
     # s2_unc = s2.hat[k]
     # s3_unc = s3.hat[k]

    s0_unc = s0[,k]
    s1_unc = s1[,k]
    s2_unc = s2[,k]
     s3_unc = s3[,k]

      for(j in 1:500){
        
        # - ++draw the 15 years of reproductive success estimates  ----
        
     #   draws=sample(1:45000,15)
        
      #  rows <- draws
      #  cols <- 1:15
        call <- cbind(rows,cols)
      #  
        # - ++draw the 15 years of reproductive success estimates  ----
     #   y_t = y_t_unc[call]
      #  y_t2 = y_t_unc[draws,]
        
        y_t2 = y_t#*runif(n=15,.8,1.2)
        
          draws=sample(1:45000,1)
         # 
          s0_est = s0_unc[draws]
          s1_est = s1_unc[draws]
          s2_est = s2_unc[draws]
          s3_est = s3_unc[draws]
        
          perturb = runif(n=1,.9,1.1)
       #  s0_est = ifelse(s0_unc*runif(n=1,.9,1.1)>1,1,s0_unc*runif(n=1,.9,1.1))
       #  s1_est = ifelse(s1_unc*perturb>1,1,s1_unc*perturb)
       #  s2_est = ifelse(s2_unc*perturb>1,1,s2_unc*perturb)
       #  s3_est = ifelse(s3_unc*perturb>1,1,s3_unc*perturb)
         
        # - ++draw 10000 samples for reproductive success with replacement  ----
        y_t.resample = sample(y_t2,10000,replace=TRUE)
        
        params = c(s0_est,s1_est,s2_est,s3_est)
        
        tmp = optim(c(.5),fr1,pars=params,pars_rs=y_t.resample,method='Brent',lower=0,upper=1)
        vec[j]=tmp$par
       # if (tmp$par>.9999) stop("Urgh, the iphone is in the blender !")
      }
      vec.list[[k]] = vec
}
 #hist(vec,breaks=100)
 #abline(v=optima[1],col='red',lwd=2)
    

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>%
  dplyr::select(site,easting,northing,elevation) %>%
  dplyr::mutate(easting=easting/1000,northing=northing/1000)
siteNames=position$site
geo=order(position$easting)


#dev.off()

vec.list.mode <- (lapply(vec.list,hdi,credMass=.68))
hdis=do.call(rbind,vec.list.mode)

post.mode.optima = unlist(lapply(vec.list,posterior.mode))
plot(post.mode.optima,g1.hat,xlim=c(.7,1),ylim=c(0,.3))
segments(x0=hdis[,1],x1=hdis[,2],y0=g1.hat)


par(mfrow=c(4,5),mar=c(1,1,1,1),
    oma=c(4,4,1,1))
for(i in 1:20){
  k = geo[i]
  
  hist(vec.list[[k]],breaks=30,main='',xlim=c(0,1))
  abline(v=optima[k],col='red')
  
  # plot(s.seq,vec.list[[k]],ylim=c(0,1),type='l',xaxt='n',yaxt='n')
  # # get closest value in vector
  # # https://stackoverflow.com/questions/43472234/fastest-way-to-find-nearest-value-in-vector  
  # index = which.min(abs(s.seq - s0.hat[k]))
  # 
  # segments(x0=s0.hat[k],y0=-1,y1=vec.list[[k]][index],lty='dotted')
  # 
  # # add estimates
  # points(x=s0.hat[k],y=0,pch=16,cex=1.2)
  # segments(x0=HPDI.s0[1,k],x1=HPDI.s0[2,k],y0=0)
  
  legend("topleft",siteNames[k],bty='n')
  
  # ifelse(i%in%c(1,6,11,16),axis(2L,las=1),NA)
  # ifelse(i%in%c(16:20),axis(1L),NA)
}

mtext("optimal g", side = 1, outer = TRUE, line = 2)

mcmc.draw = sample(1:length(s0[,1]),1000)

par(mfrow=c(4,5),mar=c(1,1,1,1),
    oma=c(4,4,1,1))
for(k in 1:20){

  # - ++get the population index  ----
  pop.index=index[,1]==k
  y_t_unc = rs[mcmc.draw,pop.index]
  
  # - ++draw the 15 years of reproductive success estimates  ----
   y_t_mode = rs.hat[pop.index]
  
  plot(y_t_mode,type='n',ylim=c(0,max(y_t_unc)))
  
  # - ++draw 10000 samples for reproductive success with replacement  ----
  # y_t.resample = sample(y_t,10000,replace=TRUE)
  
  for(j in 1:1000){
    
    # - ++draw the 15 years of reproductive success estimates  ----
    
    # - ++draw the 15 years of reproductive success estimates  ----
    y_t = y_t_unc[j,]
    
    # - ++draw 10000 samples for reproductive success with replacement  ----
   # y_t.resample = sample(y_t,10000,replace=TRUE)
    
    points(1:15+rnorm(15,mean=0,sd=.05),y_t,pch=16,col='black',cex=.5)
  }
  points(1:15,y_t_mode,col='red',pch=16)
}
