####
####
# Script to fit model for fruits per plant
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

tmpDataDirectory = "outputs/001_prepareDataForModels/"
modelDataDirectory = "outputs/002_fitStatisticalModels/data/"
modelDirectory = "models/"
outMCMCDirectory = "outputs/002_fitStatisticalModels/mcmcSamples/"

# - Libraries ----
library(rjags) 
library(parallel)
library(stringr)

# - Read in and process data ----
data <- readRDS(file=paste0(tmpDataDirectory,"fruitsPerPlantAllPlants.RDS"))

# - +Filter out unused variables ----
# none

# - +Save data for model fitting ----
saveRDS(data,file=paste0(modelDataDirectory,"fruitsData.RDS"))

# - Set up model ----

# - +MCMC conditions ----
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 10000/1
n.update = 15000/1
n.iterations = 15000/1
n.thin = 1

# - +Initial values ----
# - ++Functions to draw from prior ----

initsMu0 <- function(samps = data$n_site){
  rgamma(n = samps, shape = 1, rate = 1)
}

initsSigma0 <- function(samps = data$n_site){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

# initsSigma <- function(samps = data$n_siteYearIndex_tot){
#   extraDistr::rhnorm(n = samps, sigma = 1)
# }

# priors for proportion damaged
initsMu0_dam <- function(samps = data$n_site2){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0_dam <- function(samps = data$n_site2){
  extraDistr::rhnorm(n = samps, sigma = 2)
}

# initsSigma_dam <- function(samps = data$n_siteYearIndex_tot){
#   extraDistr::rhnorm(n = samps, sigma = 1)
# }

# - +Set inits ----
# function to initialize model in parallel

initsFun = function(){
  temp = list(initsMu0(), initsMu0(), 
              initsSigma0(), initsSigma0(), 
          #    initsSigma(samps = data$n_siteYearIndex_tfe),
          #    initsSigma(samps = data$n_siteYearIndex_tot),
              initsMu0_dam(), initsSigma0_dam(), #initsSigma_dam(samps=data$n_siteYearIndex_tot),
              
              "base::Mersenne-Twister", runif(1, 1, 2000))
  
  names(temp) = c(paste(rep("nu",2),c("tfe","tot"),sep="_"),
                  paste(rep("sigma0",2),c("tfe","tot"),sep="_"),
               #   paste(rep("sigma",2),c("tfe","tot"),sep="_"),
                  "mu0",
                  "sigma0",
                  ".RNG.name",".RNG.seed")
  return(temp)
}

inits = list(initsFun(),initsFun(),initsFun())

# - +Variables to monitor ----
# variables in models
parsToMonitor = c(paste(rep("nu",2),c("tfe","tot"),sep="_"),
                  paste(rep("sigma0",2),c("tfe","tot"),sep="_"),
                  paste(rep("mu",2),c("tfe","tot"),sep="_"),
                #  paste(rep("sigma",2),c("tfe","tot"),sep="_"),
                  "mu0","sigma0","mu")

# variables for model checking
parsToCheck = c("y_tfe.sim","chi2.tfe.obs","chi2.tfe.sim",
                       "y_tot.sim","chi2.tot.obs","chi2.tot.sim",
                       "y_dam.sim","chi2.dam.obs","chi2.dam.sim")

# - +Parallelize computation ----

cl <- makeCluster(3)
myWorkers <- NA
for(i in 1:3) myWorkers[i] <- stringr::word(capture.output(cl[[i]]), -1)

parallel::clusterExport(cl, c("myWorkers","data", "inits", "n.adapt", "n.update",
                              "n.iterations", "parsToMonitor", "parsToCheck", "modelDirectory"))

# - Run model ----
# - +Call JAGS in parallel ----

out <- clusterEvalQ(cl, {
  library(rjags)
  # create object for sampling from posterior
  jm = jags.model(file=paste0(modelDirectory,"jags-fruitspropSimple.R"),
                  data = data, n.chains = 1,
                  n.adapt = n.adapt, 
                  inits = inits[[which(myWorkers==Sys.getpid())]])
  update(jm, n.iter = n.update)
  # sample from the posterior
  zm = coda.samples(jm, variable.names = c(parsToMonitor,parsToCheck),
                    n.iter = n.iterations, thin = 1)
  return(as.mcmc(zm))
})

# end parallel computation
stopCluster(cl)

# - +Write MCMC samples to list ----
samples.rjags = mcmc.list(out)

# - Save MCMC samples ----

dir.create(file.path(outMCMCDirectory), showWarnings = FALSE)
saveRDS(samples.rjags,file=paste0(outMCMCDirectory,"fruitsPerPlantSamplesSimple.RDS"))
