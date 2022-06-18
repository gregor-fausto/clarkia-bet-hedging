####
####
# Script to fit model for seedling survival to fruiting
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
data <- readRDS(file=paste0(tmpDataDirectory,"seedlingFruitingPlantCountsPermanentPlots.RDS"))

# - +Filter out unused variables ----

# variables not used in the model
notModelVariables <- c("transect","n_transect","position","n_position",
                    "year",'n_year',"year_observed","transect_observed",
                    "siteYearIndex_observed","siteTransectIndex_observed")

data <- data[!(names(data) %in% notModelVariables)]  

# - +Save data for model fitting ----
saveRDS(data,file=paste0(modelDataDirectory,"seedlingSurvivalData.RDS"))

# - Set up model ----

# - +MCMC conditions ----
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 10000
n.update = 15000
n.iterations = 15000
n.thin = 1

# - +Initial values ----
# - ++Functions to draw from prior ----

initsMu0 <- function(samps = data$n_site){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0 <- function(samps = data$n_site){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

# - +Set inits ----
# function to initialize model in parallel
initsFun = function(){
  temp = list(initsMu0(), initsSigma0(), 
              "base::Mersenne-Twister", runif(1, 1, 2000))
  
  names(temp) = c("mu0",
                  "sigma0",
                  ".RNG.name",".RNG.seed")
  return(temp)
}

# set initial values
inits = list(initsFun(),initsFun(),initsFun())

# - +Variables to monitor ----

parsToMonitor = c("mu0","sigma0","mu")
parsToCheck = c("fruitplNumber_sim","chi2.obs","chi2.sim")

# - +Parallelize computation ----

cl <- makeCluster(3)
myWorkers <- NA
for(i in 1:3) myWorkers[i] <- stringr::word(capture.output(cl[[i]]), -1)

parallel::clusterExport(cl, c("myWorkers","data", "inits", "n.adapt", "n.update",
                              "n.iterations", "parsToMonitor", "parsToCheck","modelDirectory"))

# - Run model ----

# - +Call JAGS in parallel ----

out <- clusterEvalQ(cl, {
  library(rjags)
  jm = jags.model(file=paste0(modelDirectory,"jags-seedlingSurvival.R"),
                  data = data, n.chains = 1,
                  n.adapt = n.adapt, 
                  inits = inits[[which(myWorkers==Sys.getpid())]])
  update(jm, n.iter = n.update)
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
saveRDS(samples.rjags,file=paste0(outMCMCDirectory,"seedlingSurvivalPosteriorSamples.RDS"))
