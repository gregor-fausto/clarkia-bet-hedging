####
####
# Script to fit model for viability trials
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
data <- readRDS(file=paste0(tmpDataDirectory,"viabilityData.RDS"))

# - +Filter out variables not used in fitting the model ----

# variables not used in the model
notModelVariables <- c("yearViab","n_yearViab","ageViab","bag","n_bag",
                    "yearViab_v","n_yearViab_v","ageViab_v","bag_v","n_bag_v")

data <- data[!(names(data) %in% notModelVariables)]  

# - +Save data used for model ----

saveRDS(data, file=paste0(modelDataDirectory,"viabilityData.RDS"))

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
# - ++Functions to sample from prior ----

initsMu0 <- function(rows = data$n_siteViab, cols = data$n_yearViab){
  matrix( rnorm(n = rows*cols, mean = 0, sd = 1), rows, cols)
}

initsSigma0 <- function(rows = data$n_siteViab, cols = data$n_yearViab){
  matrix( extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
}

# - +Set inits ----
# function to initialize model in parallel
initsFun = function(){
  temp = list(initsMu0(cols=2),initsMu0(cols=2), 
              initsSigma0(cols=2), initsSigma0(cols=2), 
              "base::Mersenne-Twister", runif(1, 1, 2000))

  names(temp) = c("mu0_g","mu0_v",
                  "sigma0_g","sigma0_v",
                  ".RNG.name",".RNG.seed")
  return(temp)
}

# set initial values
inits = list(initsFun(),initsFun(),initsFun())

# - +Variables to monitor ----

parsToMonitor_g = c("mu0_g","sigma0_g","mu_g","germCount_sim","chi2.germCount.obs","chi2.germCount.sim")
parsToMonitor_v = c("mu0_v","sigma0_v","mu_v","viabStain_sim","chi2.viabStain.obs","chi2.viabStain.sim")

# - +Parallelize computation ----

cl <- makeCluster(3)
myWorkers <- NA
for(i in 1:3) myWorkers[i] <- stringr::word(capture.output(cl[[i]]), -1)

parallel::clusterExport(cl, c("myWorkers","data", "inits", "n.adapt", "n.update",
                              "n.iterations", "parsToMonitor_g", "parsToMonitor_v","modelDirectory"))

# - Run model ----
# - +Call JAGS in parallel ----
out <- clusterEvalQ(cl, {
  library(rjags)
  # create object for sampling from posterior
  jm = jags.model(file=paste0(modelDirectory,"jags-viabilityTrials.R"),
                  data = data, n.chains = 1,
                  n.adapt = n.adapt,
                  inits = inits[[which(myWorkers==Sys.getpid())]])
  update(jm, n.iter = n.update)
  # sample from the posterior
  zm = coda.samples(jm, variable.names = c(parsToMonitor_g,parsToMonitor_v),
                    n.iter = n.iterations, thin = 1)
  return(as.mcmc(zm))
})

# end parallel computation
stopCluster(cl)

# - +Write MCMC samples to list ----
samples.rjags = mcmc.list(out)

# - Save MCMC samples ----
dir.create(file.path(outMCMCDirectory), showWarnings = FALSE)
saveRDS(samples.rjags,file=paste0(outMCMCDirectory,"viabilityPosteriorSamples.RDS"))
