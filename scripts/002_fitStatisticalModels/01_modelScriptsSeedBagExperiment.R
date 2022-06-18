####
####
# Script to fit model for seed survival and germination
# to observations from seed bag experiment
# and seed rain + seedlings
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
data <- readRDS(file=paste0(tmpDataDirectory,"seedData.RDS"))

# read in matrix used for initialization
matrixForInitialization=readRDS(file=paste0(tmpDataDirectory,"seedBagsMatrixForInitialization.RDS"))

# - +Filter out variables not used in fitting the model ----

# variables used in the model
modelVariables <- c("siteSurvival","yearSurvival","y","seedStart","siteGermination",
                    "germinationIndex","seedlingJan","totalJan","sitePlot","yearPlot",
                    "fecIndex","plotSeedlings","fec","n1","n2","n3","fecCompIndex",
                    "n_t2","n_t1",
                    "n_siteGermination","n_yearGermination",
                    "n_siteSurvival","n_yearSurvival","compIndex","months")

data <- data[names(data) %in% modelVariables]

# - +Save data used for model ----

saveRDS(data, file=paste0(modelDataDirectory,"seedData.RDS"))

# - +Filter out unused variables ----

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
# - ++Functions to set initial conditions from prior ----

initsMu0Mat <- function(rows = data$n_siteGermination, cols = data$n_yearGermination){
  matrix(rnorm(n = rows*cols, mean = 0, sd = 1), rows, cols)
}

initsSigma0LessWeakMat <- function(rows = data$n_siteSurvival, cols = data$n_yearSurvival){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 2), rows, cols)
}

initsMu0 <- function(samps = data$n_siteSurvival){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0LessWeak <- function(samps = data$n_siteSurvival){
  extraDistr::rhnorm(n = samps, sigma = 2)
}

f=function(x){(sample(1:x,1))}

f2=function(x){
  x=as.matrix(x)
  if(x[1]==0){
    c(0,x[3])
  } else if(x[2]==0){
    c(x[3],0)
  } else if(x[1]!=0&x[2]!=0&x[3]!=0){
    tmp=f(x[3])
    c(tmp,x[3]-tmp)
  } else if(x[1]!=0&x[2]!=0&x[3]==0){
    c(0,0)
  }
}

# - +Set inits ----
# function to initialize model in parallel
initsFun = function(){

  tmp=matrixForInitialization[,3:5]
  initsSeedlings<-apply(tmp,1,f2)

  temp = list(initsMu0Mat(rows=data$n_siteSurvival, cols=4), 
              initsSigma0LessWeakMat(rows=data$n_siteSurvival,cols=4),
              initsMu0Mat(cols=2), initsSigma0LessWeakMat(cols=2), 
              initsMu0(), initsSigma0LessWeak(),
              initsSeedlings[2,],  initsSeedlings[1,],
              "base::Mersenne-Twister", runif(1, 1, 2000))

  names(temp) = c("mu0_s","sigma0_s",
                  "mu0_g", "sigma0_g",
                  "mu0_s0", "sigma0_s0",
                  "y_t1", "y_t2",
                  ".RNG.name",".RNG.seed")
  return(temp)
}

# set initial values
inits = list(initsFun(),initsFun(),initsFun())

# - +Variables to monitor ----
# model parameters
parsToMonitor = c("mu0_s","sigma0_s","mu_s",
                  "mu0_g", "sigma0_g","mu_g",
                  "mu0_s0", "sigma0_s0", "mu_s0",
                  "s0_1","s0_2")
# model checking
parsToCheck = c("y_sim","chi2.yobs","chi2.ysim",
                "seedlingJan_sim","chi2.obs","chi2.sim",
                "plotSeedlings_sim","chi2.plot.obs","chi2.plot.sim",
                "y_t1","y_t2", "s0_1", "s0_2")

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
  jm = jags.model(file=paste0(modelDirectory,"jags-seedBagExperiment.R"),
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
saveRDS(samples.rjags,file=paste0(outMCMCDirectory,"seedBagExperimentPosteriorSamples.RDS"))
