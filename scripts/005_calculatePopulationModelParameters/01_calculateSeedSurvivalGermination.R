####
####
# Script to calculate parameters for population model
# Seed survival, germination
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

mcmcDirectory = "outputs/005_calculatePopulationModelParameters/01_parameterPosteriorDistributions/"
outputDirectoryPopulation = "outputs/005_calculatePopulationModelParameters/02_populationModelParameters/"
outputDirectoryMatrix = "outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/"

# - +load libraries ----
library(MCMCvis)
library(tidyverse)
library(magrittr)

# - Read in what's needed ----

# - +Read in MCMC samples ----
mcmcSampleDirectory <- paste0(mcmcDirectory,list.files(mcmcDirectory))
mcmcSamplesSeedBags <- readRDS(mcmcSampleDirectory[[grep("seedBagExperimentParameters.RDS",mcmcSampleDirectory)]])
mcmcSamplesViability <- readRDS(mcmcSampleDirectory[[grep("viabilityTrialsParameters.RDS",mcmcSampleDirectory)]])

# - +Read in site names ----
siteAbiotic <- read.csv("data/siteAbioticData.csv",header=TRUE)
siteNames <- siteAbiotic$site

# ---
# - Seed bag burial posteriors ----
# ---

# - +Germination: population-level ----

mu0_g=MCMCchains(mcmcSamplesSeedBags,params="mu0_g")
mu_g=MCMCchains(mcmcSamplesSeedBags,params="mu_g")

# - ++index for age 1, 2, 3 ----
index.g1=grep(",1\\]",colnames(mu0_g))
index.g2=grep(",2\\]",colnames(mu0_g))
index.g3=grep(",3\\]",colnames(mu_g))

# - ++transform to probability scale and summarize ----
gamma1 = boot::inv.logit(mu0_g[,index.g1])
gamma1.sum=apply(gamma1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

gamma2 = boot::inv.logit(mu0_g[,index.g2])
gamma2.sum=apply(gamma2,2,quantile,probs=c(0.025,.25,.5,.75,.975))

gamma3 = boot::inv.logit(mu_g[,index.g3])
gamma3.sum=apply(gamma3,2,quantile,probs=c(0.025,.25,.5,.75,.975))

# - +Survival: population-level ----

mu0_s=MCMCchains(mcmcSamplesSeedBags,params="mu0_s")

index.s1=grep(",1\\]",colnames(mu0_s))
index.s2=grep(",2\\]",colnames(mu0_s))
index.s3=grep(",3\\]",colnames(mu0_s))

sigma1 = boot::inv.logit(mu0_s[,index.s1])
sigma1.sum=apply(sigma1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma2 = boot::inv.logit(mu0_s[,index.s2])
sigma2.sum=apply(sigma2,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma3 = boot::inv.logit(mu0_s[,index.s3])
sigma3.sum=apply(sigma3,2,quantile,probs=c(0.025,.25,.5,.75,.975))

# ---
# - Viability experiment posteriors ----
# ---

# - +Population estimates ----
# germination tests
mu0_g<-MCMCchains(mcmcSamplesViability, params = "mu0_g")
# viability assays
mu0_v<-MCMCchains(mcmcSamplesViability, params = "mu0_v")
# calculate total viability; calculation on latent scale
nu0=exp(mu0_g)/(1+exp(mu0_g))+(exp(-mu0_g+mu0_v)/(1+exp(-mu0_g)+exp(mu0_v)+exp(-mu0_g+mu0_v)))

# - +Interpolate estimates ----
# - ++January, age 0 seeds ----
nu0_ratio1=(nu0[,index.g1])^(1/3)
# - ++January, age 1 seeds ----
nu0_ratio2=nu0[,index.g1]*(nu0[,index.g2]/nu0[,index.g1])^(1/3)

# - +Relabel ----
# - ++January, age 0 seeds ----
nu_01.inter = nu0_ratio1
# - ++October, age 0 seeds ----
nu_1 = nu0[,index.g1]
# - ++January, age 1 seeds ----
nu_12.inter = nu0_ratio2
# - ++October, age 1 seeds ----
nu_2 = nu0[,index.g2]

# ---
# - Combine field and viability posteriors ----
# ---

# - +Germination, conditional on viability ----

g1.v = gamma1/(1-(1-nu_01.inter)*(1-gamma1))

# - +Survival, conditional on viability ----

# - ++s1 ----

s1.v = sigma1*(gamma1+(1-gamma1)*nu_01.inter)

# - ++s2 ----

s2.v = sigma2*(nu_1^(2/3))

# - ++s3 ----

s3.v = ((sigma1*(1-gamma1)*sigma2*sigma3)*(gamma2+(1-gamma2)*nu_12.inter))/(s1.v*(1-g1.v)*s2.v)

# function to resample values greater than 0 and redistribute probability mass 
resample=function(x){
  tmp = ifelse(x<=1,x,NA)
  len=sum(is.na(tmp))
  tmp2=x[x<=1]
  out = ifelse(is.na(tmp),sample(tmp2,len,replace=TRUE),tmp)
  return(out)
}

s3.v=apply(s3.v,2,resample)


# ---
# - Obtain estimates for s0 ----
# ---

# population level estimate
mu0_s0 = MCMCchains(mcmcSamplesSeedBags,params="mu0_s0")
s0 = boot::inv.logit(mu0_s0)

# population*year level estimate
mu_s0 = MCMCchains(mcmcSamplesSeedBags,params="mu_s0")
mu_s0 = boot::inv.logit(mu_s0)

# - Summarize parameters ----

g1.v.sum=apply(g1.v,2,quantile,probs=c(0.025,.25,.5,.75,.975))
s1.v.sum=apply(s1.v,2,quantile,probs=c(0.025,.25,.5,.75,.975))
s2.v.sum=apply(s2.v,2,quantile,probs=c(0.025,.25,.5,.75,.975))
s3.v.sum=apply(s3.v,2,quantile,probs=c(0.025,.25,.5,.75,.975))
s0.sum=apply(s0,2,quantile,probs=c(0.025,.25,.5,.75,.975))

# ---
# - Save structured model parameters ----
# ---

saveRDS(g1.v,paste0(outputDirectoryPopulation,"g1-population-level.RDS"))
saveRDS(s1.v,paste0(outputDirectoryPopulation,"s1-population-level.RDS"))
saveRDS(s2.v,paste0(outputDirectoryPopulation,"s2-population-level.RDS"))
saveRDS(s3.v,paste0(outputDirectoryPopulation,"s3-population-level.RDS"))
saveRDS(s0,paste0(outputDirectoryPopulation,"s0-population-level.RDS"))
saveRDS(mu_s0,paste0(outputDirectoryPopulation,"s0-population-year-level.RDS"))


# ---
# - Seed bag burial posteriors: population+year level for all parameters ----
# ---

# - +Germination: population-level ----

mu_g=MCMCchains(mcmcSamplesSeedBags,params="mu_g")

# - ++index for age 1, 2, 3 ----
index.g11=grep(",1\\]",colnames(mu_g))
index.g21=grep(",2\\]",colnames(mu_g))
index.g31=grep(",3\\]",colnames(mu_g))
index.g12=grep(",4\\]",colnames(mu_g))
index.g22=grep(",5\\]",colnames(mu_g))
index.g13=grep(",6\\]",colnames(mu_g))

# - ++transform to probability scale and summarize ----
gamma11 = boot::inv.logit(mu_g[,index.g11])
gamma21 = boot::inv.logit(mu_g[,index.g21])
gamma31 = boot::inv.logit(mu_g[,index.g31])
gamma12 = boot::inv.logit(mu_g[,index.g12])
gamma22 = boot::inv.logit(mu_g[,index.g22])
gamma13 = boot::inv.logit(mu_g[,index.g13])

# - +Survival: population-level ----

mu_s=MCMCchains(mcmcSamplesSeedBags,params="mu_s")

index.s11=grep(",1\\]",colnames(mu_s))
index.s21=grep(",2\\]",colnames(mu_s))
index.s31=grep(",3\\]",colnames(mu_s))
index.s41=grep(",4\\]",colnames(mu_s))
index.s51=grep(",5\\]",colnames(mu_s))
index.s61=grep(",6\\]",colnames(mu_s))

index.s12=grep(",7\\]",colnames(mu_s))
index.s22=grep(",8\\]",colnames(mu_s))
index.s32=grep(",9\\]",colnames(mu_s))
index.s42=grep(",10\\]",colnames(mu_s))

index.s13=grep(",11\\]",colnames(mu_s))
index.s23=grep(",12\\]",colnames(mu_s))

sigma11 = boot::inv.logit(mu_s[,index.s11])
sigma21 = boot::inv.logit(mu_s[,index.s21])
sigma31 = boot::inv.logit(mu_s[,index.s31])
sigma41 = boot::inv.logit(mu_s[,index.s41])
sigma51 = boot::inv.logit(mu_s[,index.s51])
sigma61 = boot::inv.logit(mu_s[,index.s61])

sigma12 = boot::inv.logit(mu_s[,index.s12])
sigma22 = boot::inv.logit(mu_s[,index.s22])
sigma32 = boot::inv.logit(mu_s[,index.s32])
sigma42 = boot::inv.logit(mu_s[,index.s42])

sigma13 = boot::inv.logit(mu_s[,index.s13])
sigma23 = boot::inv.logit(mu_s[,index.s23])

# ---
# - Viability experiment posteriors ----
# ---

# - +Population estimates ----
# germination tests
mu_g<-MCMCchains(mcmcSamplesViability, params = "mu_g")
# viability assays
mu_v<-MCMCchains(mcmcSamplesViability, params = "mu_v")
# calculate total viability; calculation on latent scale
nu=exp(mu_g)/(1+exp(mu_g))+(exp(-mu_g+mu_v)/(1+exp(-mu_g)+exp(mu_v)+exp(-mu_g+mu_v)))

# - +Interpolate estimates ----
# - ++January, age 0 seeds ----
nu_ratio11=(nu[,index.g11])^(1/3)
# - ++January, age 1 seeds ----
nu_ratio21=nu[,index.g11]*(nu[,index.g21]/nu[,index.g11])^(1/3)
# - ++January, age 2 seeds ----
nu_ratio31=nu[,index.g21]*(nu[,index.g31]/nu[,index.g21])^(1/3)

# - ++January, age 0 seeds ----
nu_ratio12=(nu[,index.g12])^(1/3)
# - ++January, age 1 seeds ----
nu_ratio22=nu[,index.g12]*(nu[,index.g22]/nu[,index.g12])^(1/3)

# - ++January, age 0 seeds ----
nu_ratio13=(nu[,index.g13])^(1/3)

# - +Relabel ----
# - ++October, age 0 seeds ----
nu_1 = nu[,index.g11]
# - ++October, age 1 seeds ----
nu_2 = nu[,index.g21]
# - ++October, age 1 seeds ----
nu_3 = nu[,index.g31]

# - ++October, age 0 seeds ----
nu_4 = nu[,index.g12]
# - ++October, age 1 seeds ----
nu_5 = nu[,index.g22]

# - ++October, age 1 seeds ----
nu_6 = nu[,index.g13]

# ---
# - Combine field and viability posteriors ----
# ---

# - ROUND 1 ----

# - +Germination, conditional on viability ----

g11.v = gamma11/(1-(1-nu_ratio11)*(1-gamma11))

# - +Survival, conditional on viability ----

# - ++s11 ----

s11.v = sigma11*(gamma11+(1-gamma11)*nu_ratio11)

# - ++s21 ----

s21.v = sigma21*(nu_1^(2/3))

# - ++s31 ----

s31.v = ((sigma11*(1-gamma11)*sigma21*sigma31)*(gamma21+(1-gamma21)*nu_ratio21))/(s11.v*(1-g11.v)*s21.v)

# - ++s41 ----

s41.v = sigma41*(nu_2^(2/3))

# - +Germination, conditional on viability ----

g21.v = gamma21/(1-(1-nu_ratio21)*(1-gamma21))

# - ++s51 ----

s51.v = ((sigma11*(1-gamma11)*sigma21*sigma31*(1-gamma21)*sigma41*sigma51)*(gamma31+(1-gamma31)*nu_ratio31))/(s11.v*(1-g11.v)*s21.v*s31.v*(1-g21.v)*s41.v)

# - ++s61 ----

s61.v = sigma61*(nu_3^(2/3))

# - +Germination, conditional on viability ----

g31.v = gamma31/(1-(1-nu_ratio31)*(1-gamma31))

# - ROUND 2 ----

# - +Germination, conditional on viability ----

g12.v = gamma12/(1-(1-nu_ratio12)*(1-gamma12))

# - +Survival, conditional on viability ----

# - ++s11 ----

s12.v = sigma12*(gamma12+(1-gamma12)*nu_ratio12)

# - ++s21 ----

s22.v = sigma22*(nu_4^(2/3))

# - ++s31 ----

s32.v = ((sigma12*(1-gamma12)*sigma22*sigma32)*(gamma22+(1-gamma22)*nu_ratio22))/(s12.v*(1-g12.v)*s22.v)

# - ++s41 ----

s42.v = sigma42*(nu_5^(2/3))

# - +Germination, conditional on viability ----

g22.v = gamma22/(1-(1-nu_ratio22)*(1-gamma22))

# - ROUND 3 ----

# - +Germination, conditional on viability ----

g13.v = gamma13/(1-(1-nu_ratio13)*(1-gamma13))

# - +Survival, conditional on viability ----

# - ++s11 ----

s13.v = sigma13*(gamma13+(1-gamma13)*nu_ratio13)

# - ++s21 ----

s23.v = sigma23*(nu_6^(2/3))

# lapply(list(g11.v,g21.v,g31.v,g12.v,g22.v,g13.v),max) %>% unlist
# lapply(list(s11.v,s21.v,s31.v,s41.v,s51.v,s61.v,s12.v,s22.v,s32.v,s42.v,s13.v,s23.v),max) %>% unlist

# function to resample values greater than 0 and redistribute probability mass 
resample=function(x){
  tmp = ifelse(x<=1,x,NA)
  len=sum(is.na(tmp))
  tmp2=x[x<=1]
  if(sum(tmp2)>0){
    out = ifelse(is.na(tmp),sample(tmp2,len,replace=TRUE),tmp)
  } else if(sum(tmp2)==0){
    out = rep(1,len)
  }
  return(out)
}

# the first 2 of these include estimate(s) on the boundary
s31.v=apply(s31.v,2,resample)
s51.v=apply(s51.v,2,resample)
s32.v=apply(s32.v,2,resample)


# ---
# - Obtain estimates for s0 ----
# ---

# population*year level estimate
mu_s0 = MCMCchains(mcmcSamplesSeedBags,params="mu_s0")
mu_s0 = boot::inv.logit(mu_s0)

# - +s0:convert to matrix ----

mu_s0.2006 = mu_s0[,grep(",1]",colnames(mu_s0))]
mu_s0.2007 = mu_s0[,grep(",2]",colnames(mu_s0))]

n.iter=dim(mu_s0.2006)[1]
dat.list <- list()
for(i in 1:20){
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  mat[,2] = mu_s0.2006[,i]
  mat[,3] = mu_s0.2007[,i]
  dat.list[[i]] <- mat
}

saveRDS(dat.list,paste0(outputDirectoryMatrix,"s0-ex1-population-year-level-mat.RDS"))

# - +s1:convert to matrix ----

dat.list <- list()
for(i in 1:20){
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  mat[,1] = s11.v[,i]
  mat[,2] = s12.v[,i]
  mat[,3] = s13.v[,i]
  dat.list[[i]] <- mat
}

saveRDS(dat.list,paste0(outputDirectoryMatrix,"s1-ex1-population-year-level-mat.RDS"))

# - +s2:convert to matrix ----

dat.list <- list()
for(i in 1:20){
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  mat[,1] = s21.v[,i]
  mat[,2] = s22.v[,i]
  mat[,3] = s23.v[,i]
  dat.list[[i]] <- mat
}

saveRDS(dat.list,paste0(outputDirectoryMatrix,"s2-ex1-population-year-level-mat.RDS"))

# - +s3:convert to matrix ----

dat.list <- list()
for(i in 1:20){
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  mat[,2] = s31.v[,i]
  mat[,3] = s32.v[,i]
  dat.list[[i]] <- mat
}

saveRDS(dat.list,paste0(outputDirectoryMatrix,"s3-ex1-population-year-level-mat.RDS"))

# - +s4:convert to matrix ----

dat.list <- list()
for(i in 1:20){
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  mat[,2] = s41.v[,i]
  mat[,3] = s42.v[,i]
  dat.list[[i]] <- mat
}

saveRDS(dat.list,paste0(outputDirectoryMatrix,"s4-ex1-population-year-level-mat.RDS"))

# - +s5:convert to matrix ----

dat.list <- list()
for(i in 1:20){
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  mat[,3] = s51.v[,i]
  dat.list[[i]] <- mat
}

saveRDS(dat.list,paste0(outputDirectoryMatrix,"s5-ex1-population-year-level-mat.RDS"))

# - +s6:convert to matrix ----

dat.list <- list()
for(i in 1:20){
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  mat[,3] = s61.v[,i]
  dat.list[[i]] <- mat
}

saveRDS(dat.list,paste0(outputDirectoryMatrix,"s6-ex1-population-year-level-mat.RDS"))


# - +g1:convert to matrix ----

dat.list <- list()
for(i in 1:20){
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  mat[,1] = g11.v[,i]
  mat[,2] = g12.v[,i]
  mat[,3] = g13.v[,i]
  dat.list[[i]] <- mat
}

saveRDS(dat.list,paste0(outputDirectoryMatrix,"g1-ex1-population-year-level-mat.RDS"))

# - +g2:convert to matrix ----

dat.list <- list()
for(i in 1:20){
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  mat[,2] = g21.v[,i]
  mat[,3] = g22.v[,i]
  dat.list[[i]] <- mat
}

saveRDS(dat.list,paste0(outputDirectoryMatrix,"g2-ex1-population-year-level-mat.RDS"))

# - +g3:convert to matrix ----

dat.list <- list()
for(i in 1:20){
  mat = matrix(NA,nrow = n.iter,ncol = 15)
  mat[,3] = g31.v[,i]
  dat.list[[i]] <- mat
}

saveRDS(dat.list,paste0(outputDirectoryMatrix,"g3-ex1-population-year-level-mat.RDS"))
