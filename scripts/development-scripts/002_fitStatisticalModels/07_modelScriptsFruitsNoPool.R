####
####
# Script to fit model for fruits per plant
# *this model does not pool years
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
saveRDS(data,file=paste0(modelDataDirectory,"fruitsNoPoolData.RDS"))

# - Set up model ----

# - +MCMC conditions ----
# scalars that specify the 
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
n.adapt = 10000/1
n.update = 15000/1
n.iterations = 15000/1
n.thin = 10

# - +Initial values ----
# - ++Functions to draw from prior ----

initsMu <- function(rows = data$n_site, cols = data$n_year){
  matrix(rgamma(n = rows*cols, shape = 1, rate = 1), rows, cols)
}

initsSigma <- function(rows = data$n_site, cols = data$n_year){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 2), rows, cols)
}

# - +Set inits ----

initsFun = function(){
  temp = list(initsMu(rows = data$n_site,cols=data$n_year), 
              initsMu(rows = data$n_site2,cols=data$n_year2), 
              initsMu(rows = data$n_site2,cols=data$n_year2), 
              initsSigma(rows = data$n_site,cols=data$n_year),
              initsSigma(rows = data$n_site2,cols=data$n_year2),
              initsSigma(rows = data$n_site2,cols=data$n_year2),
              "base::Mersenne-Twister", runif(1, 1, 2000))
  
  names(temp) = c(paste(rep("mu.log",3),c("tfe","und","dam"),sep="_"),
                  paste(rep("sigma",3),c("tfe","und","dam"),sep="_"),
                  ".RNG.name",".RNG.seed")
  return(temp)
}

# set initial values
inits = list(initsFun(),initsFun(),initsFun())

# - +Variables to monitor ----

parsToMonitor = c(paste(rep("mu.log",3),c("tfe","und","dam"),sep="_"),
                  paste(rep("sigma",3),c("tfe","und","dam"),sep="_"))
parsToCheck = c("y_tfe.sim","chi2.tfe.obs","chi2.tfe.sim",
                "y_und.sim","chi2.und.obs","chi2.und.sim",
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
  jm = jags.model(file=paste0(modelDirectory,"jags-fruits-noPool.R"),
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
saveRDS(samples.rjags,file=paste0(outMCMCDirectory,"fruitsSamples-noPool.rds"))
