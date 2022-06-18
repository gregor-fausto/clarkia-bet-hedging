
# - Environment ----

fileDirectory = "/Users/Gregor/Dropbox/dataLibrary/seed-banks-reanalysis-2/"
outDataDirectory = "/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/002_statisticalModelFitting/data/"

rm(list=setdiff(ls(all=TRUE),c("fileDirectory","outputDirectory"))) # if using in source(script)
options(stringsAsFactors = FALSE)

# - Libraries ----

library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)

# - Read in data ----

# - +Read seed bag trials & posterior ----

mcmcSampleDirectory <- paste0(fileDirectory,list.files(fileDirectory))
mcmcSamples <- readRDS(mcmcSampleDirectory[[grep("seedSamplesNonparametric-twoYear.rds",mcmcSampleDirectory)]])

# - +Read viability trials & posterior ----

mcmcSamplesViability <- readRDS(mcmcSampleDirectory[[grep("viabilityTrialSamples.rds",mcmcSampleDirectory)]])

# - +Site names ----

position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
  dplyr::select(site,easting,dominant.surface.rock.type) %>%
  dplyr::mutate(easting=easting/1000)

siteNames <- unique(position$site)


# ---
# - Seed bag burial posteriors ----
# ---

# - +Germination: population- and year-level ----

mu_g=MCMCchains(mcmcSamples,params="mu_g")

# - ++index for age 1, 2, 3 ----
index.g1=grep(c(",1\\]|4\\]|6\\]"),colnames(mu_g))
index.g2=grep(",2\\]|5\\]",colnames(mu_g))
index.g3=grep(",3\\]",colnames(mu_g))

# - ++transform to probability scale and summarize ----
gamma1 = boot::inv.logit(mu_g[,index.g1])
gamma1.sum=apply(gamma1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

gamma2 = boot::inv.logit(mu_g[,index.g2])
gamma2.sum=apply(gamma2,2,quantile,probs=c(0.025,.25,.5,.75,.975))

gamma3 = boot::inv.logit(mu_g[,index.g3])
gamma3.sum=apply(gamma3,2,quantile,probs=c(0.025,.25,.5,.75,.975))

# - +Survival: population-level ----

mu_s=MCMCchains(mcmcSamples,params="mu_s")

index.s1=grep(",1\\]|,7\\]|,11\\]",colnames(mu_s))
index.s2=grep(",2\\]|,8\\]|,12\\]",colnames(mu_s))
index.s3=grep(",3\\]|,9\\]",colnames(mu_s))
index.s4=grep(",4\\]|,10\\]",colnames(mu_s))
index.s5=grep(",5\\]",colnames(mu_s))
index.s6=grep(",6\\]",colnames(mu_s))

sigma1 = boot::inv.logit(mu_s[,index.s1])
sigma1.sum=apply(sigma1,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma2 = boot::inv.logit(mu_s[,index.s2])
sigma2.sum=apply(sigma2,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma3 = boot::inv.logit(mu_s[,index.s3])
sigma3.sum=apply(sigma3,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma4 = boot::inv.logit(mu_s[,index.s4])
sigma4.sum=apply(sigma4,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma5 = boot::inv.logit(mu_s[,index.s5])
sigma5.sum=apply(sigma5,2,quantile,probs=c(0.025,.25,.5,.75,.975))

sigma6 = boot::inv.logit(mu_s[,index.s6])
sigma6.sum=apply(sigma6,2,quantile,probs=c(0.025,.25,.5,.75,.975))

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
nu_ratio1=(nu[,index.g1])^(1/3)
# - ++January, age 1 seeds ----
nu_ratio2=nu[,index.g1[1:40]]*(nu[,index.g2[1:40]]/nu[,index.g1[1:40]])^(1/3)
# - ++January, age 2 seeds ----
nu_ratio3=nu[,index.g2[1:20]]*(nu[,index.g3[1:20]]/nu[,index.g2[1:20]])^(1/3)

# - +Relabel ----
# - ++January, age 0 seeds ----
nu_01.inter = nu_ratio1
# - ++October, age 0 seeds ----
nu_1 = nu[,index.g1]
# - ++January, age 1 seeds ----
nu_12.inter = nu_ratio2
# - ++October, age 1 seeds ----
nu_2 = nu[,index.g2]
# - ++January, age 2 seeds ----
nu_23.inter = nu_ratio3
# - ++October, age 2 seeds ----
nu_3 = nu[,index.g3]

# ---
# - Combine field and viability posteriors ----
# ---

# - +Germination, conditional on viability ----

g1.v = gamma1/(1-(1-nu_01.inter)*(1-gamma1))
g2.v = gamma2/(1-(1-nu_12.inter)*(1-gamma2))
g3.v = gamma3/(1-(1-nu_23.inter)*(1-gamma3))

# - ++ save germination ----

germinationList <- list ( g1.v, g2.v, g3.v )

# - +Survival, conditional on viability ----

# - ++s1 ----

s1.v = sigma1*(gamma1+(1-gamma1)*nu_01.inter)

# - ++s2 ----

s2.v = sigma2*(nu_1^(2/3))

# - ++s3 ----

s3.v = ((sigma1[1:40]*(1-gamma1[1:40])*sigma2[1:40]*sigma3)*(gamma2+(1-gamma2)*nu_12.inter))/(s1.v[1:40]*(1-g1.v[1:40])*s2.v[1:40])

# - ++s4 ----

s4.v = sigma4*(nu_2^(2/3))

# - ++s5 ----

s5.v = ((sigma1[1:20]*(1-gamma1[1:20])*sigma2[1:20]*sigma3[1:20]*(1-gamma2[1:20])*sigma4[1:20]*sigma5[1:20])*(gamma3+(1-gamma3)*nu_23.inter))/(s1.v[1:20]*(1-g1.v[1:20])*s2.v[1:20]*s3.v[1:20]*(1-g2.v[1:20]*s4.v[1:20]))

# - ++s6 ----

s6.v = sigma6*(nu_3^(2/3))

# function to resample values greater than 0 and redistribute probability mass 
resample=function(x){
  tmp = ifelse(x<=1,x,NA)
  len=sum(is.na(tmp))
  tmp2=x[x<=1]
  out = ifelse(is.na(tmp),sample(tmp2,len,replace=TRUE),tmp)
  return(out)
}

s3.v=apply(s3.v,2,resample)

survivalList <- list ( s1.v, s2.v, s3.v , s4.v, s5.v, s6.v)

bayesianEstimates <- list(germinationList,survivalList)

saveRDS(bayesianEstimates,paste0("/Users/Gregor/Dropbox/dataLibrary/chapter-3/","seedRateEstimatesBayesianCombined.rds"))

