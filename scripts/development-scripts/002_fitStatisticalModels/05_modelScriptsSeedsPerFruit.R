####
####
# Script to fit model for seeds per fruit
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

saveRDS(data, file=paste0(modelDataDirectory,"seedsPerFruitData.RDS"))

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

initsNu <- function(samps = data$n_site3){
  rgamma(n = samps, shape = 1, rate = 1)
}

initsSigma0 <- function(samps = data$n_site3){
  extraDistr::rhnorm(n = samps, sigma = 2)
}

initsSigma <- function(samps = data$n_siteYearIndex_und){
  extraDistr::rhnorm(n = samps, sigma = 2)
}

# priors for proportion damaged
initsMu0_dam <- function(samps = data$n_site3){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0_dam <- function(samps = data$n_site3){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsSigma_dam <- function(samps = data$n_siteYearIndex_dam){
  extraDistr::rhnorm(n = samps, sigma = 1)
}



# - +Set inits ----
# function to initialize model in parallel
# initsFun = function(){temp = list(initsMu0(), initsMu0(), 
#                                   initsSigma0(), initsSigma0(), 
#                                   initsSigma(rows = data$n_site3,cols=data$n_year3),
#                                   initsSigma(rows = data$n_site3,cols=data$n_year4),
#                                   "base::Mersenne-Twister", runif(1, 1, 2000))
# 
# names(temp) = c(paste(rep("nu",2),c("seeds","dam_seeds"),sep="_"),
#                 paste(rep("sigma0",2),c("seeds","dam_seeds"),sep="_"),
#                 paste(rep("sigma",2),c("seeds","dam_seeds"),sep="_"),
#                 ".RNG.name",".RNG.seed")

initsFun = function(){temp = list(initsNu(), initsSigma0(), #initsSigma(samps=data$n_siteYearIndex_und),
                                  initsMu0_dam(), initsSigma0_dam(),# initsSigma_dam(samps=data$n_siteYearIndex_dam),
                                  "base::Mersenne-Twister", runif(1, 1, 2000))

names(temp) = c(paste(rep("nu",1),c("seeds"),sep="_"),
                paste(rep("sigma0",1),c("seeds"),sep="_"),
            #    paste(rep("sigma",1),c("seeds"),sep="_"),
                "mu0",
                "sigma0",
              #  "sigma",
                ".RNG.name",".RNG.seed")

return(temp)
}

# set initial values
inits = list(initsFun(),initsFun(),initsFun())


# - +Variables to monitor ----
# variables in model
parsToMonitor = c(paste(rep("nu",1),c("seeds"),sep="_"),
                  paste(rep("sigma0",1),c("seeds"),sep="_"),
                  "mu_seeds",
                 #paste(rep("sigma",1),c("seeds"),sep="_"),
                  "mu0",
                  "sigma0",
                  "mu")
# variables in for model checking
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
  jm = jags.model(file=paste0(modelDirectory,"jags-seedsperfruitpropSimple.R"),
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
saveRDS(samples.rjags,file=paste0(outMCMCDirectory,"seedsPerFruitSamplesSimple.RDS"))
