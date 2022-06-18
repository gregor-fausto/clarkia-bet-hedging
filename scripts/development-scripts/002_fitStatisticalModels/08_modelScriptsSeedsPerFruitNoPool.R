####
####
# Script to fit model for seeds per fruit
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
data <- readRDS(file=paste0(tmpDataDirectory,"seedsPerFruit.RDS"))

# - +Filter out variables not used in fitting the model ----

# variables not used in the model
modelVariables <- c("n")

data <- data[!(names(data) %in% modelVariables)]  

# - +Save data used for model ----

saveRDS(data, file=paste0(modelDataDirectory,"seedsPerFruitDataNoPool.rds"))

# - Set up model ----

# - +Initial values ----
# - ++Functions to draw from prior ----

initsMu <- function(rows = data$n_site3, cols = data$n_year3){
  matrix(rgamma(n = rows*cols, shape = 1, rate = 1), rows, cols)
}

initsSigma <- function(rows = data$n_site3, cols = data$n_year3){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 2), rows, cols)
}
# - +Set inits ----

initsFun = function(){temp = list(initsMu(rows = data$n_site3,cols=data$n_year3),
                                  initsMu(rows = data$n_site4,cols=data$n_year4),
                                  initsSigma(rows = data$n_site3,cols=data$n_year3),
                                  initsSigma(rows = data$n_site4,cols=data$n_year4),
                                  "base::Mersenne-Twister", runif(1, 1, 2000))

names(temp) = c(paste(rep("mu.log",2),c("seeds","dam_seeds"),sep="_"),
                paste(rep("sigma",2),c("seeds","dam_seeds"),sep="_"),
                ".RNG.name",".RNG.seed")
return(temp)
}

# set initial values
inits = list(initsFun(),initsFun(),initsFun())

# - +Variables to monitor ----

parsToMonitor = c(paste(rep("mu.log",2),c("seeds","dam_seeds"),sep="_"),
                  paste(rep("sigma",2),c("seeds","dam_seeds"),sep="_"))
parsToCheck = c("y_sd.sim","chi2.sd.obs","chi2.sd.sim",
                "y_sd_dam.sim","chi2.sd_dam.obs","chi2.sd_dam.sim")

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
  jm = jags.model(file=paste0(modelDirectory,"jags-seedsperfruit-noPool.R"),
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
saveRDS(samples.rjags,file=paste0(outMCMCDirectory,"seedsPerFruitSamples-noPool.rds"))