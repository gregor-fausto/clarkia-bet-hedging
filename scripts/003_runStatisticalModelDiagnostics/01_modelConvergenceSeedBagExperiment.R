####
####
# Script to evaluate convergence
# Seed survival and germination model
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list = setdiff(ls(all = TRUE), 
                  c("scriptConvergenceDirectory", "outMCMCDirectory", "outputDirectory")))  # if using in source(script)
options(stringsAsFactors = FALSE)

# - Libraries ----
library(MCMCvis)
library(tidybayes)
library(tidyverse)
library(magrittr)
library(bayesplot)
library(rethinking)

# - Read in MCMC samples ----
mcmcSampleFiles <- paste0(outMCMCDirectory, list.files(outMCMCDirectory))
mcmcSamples <- readRDS(mcmcSampleFiles[[grep("seedBagExperimentPosteriorSamples.RDS", mcmcSampleFiles)]])

# - Source function ----
# function produces jpegs of trace plots
source("scripts/003_runStatisticalModelDiagnostics/00_tracePlotFunction.R")


# - Convergence diagnostics: seed survival ----
outputDirectory = "outputs/003_runStatisticalModelDiagnostics/01_tracePlots/1_seedBagExperiment/"

# - +Graphical checks on trace plots: survival ----
f(x = "mu0_s", model = "seedBagExperiment")
f(x = "sigma0_s", model = "seedBagExperiment")

# - +Graphical checks on trace plots: germination ----
f(x = "mu0_g", model = "seedBagExperiment")
f(x = "sigma0_g", model = "seedBagExperiment")

# - +Graphical checks on trace plots: s0 ----
f(x = "mu0_s0", model = "seedBagExperiment")
f(x = "sigma0_s0", model = "seedBagExperiment")


# - Convergence diagnostics ----
outputDirectory = "outputs/003_runStatisticalModelDiagnostics/02_rhatDistribution/"

# - +Recover chains ----
all.chains = MCMCchains(mcmcSamples, params = c("mu0_s", "sigma0_s"), 
                        mcmc.list = TRUE)

# R-hat for all parameters
rhat = coda::gelman.diag(all.chains, confidence = 0.95)$psrf[, 1]

# - +Heidelberger and Welch convergence diagnostic ----
hd = coda::heidel.diag(all.chains)
p = c()
for (i in 1:3) {
  p[i] = sum(hd[[i]][, 1])/length(hd[[i]][, 1])
}
p = signif(p, 2)
dt = hist(coda::gelman.diag(all.chains, confidence = 0.95)$psrf[, 1], breaks = 25, 
          plot = FALSE)

# - +Plot Rhat and HW ----
par(mfrow = c(1, 1))
jpeg(filename = paste0(outputDirectory, "rhat-hw-seedBagExperiment-seedSurvival", ".jpeg"), 
     quality = 75)
# plot distribution of Rhat values
hist(rhat, col = "black", border = "white", breaks = 25,
     main = NULL,
     xlab ="R-hat values")
title(main = "A. Distribution of R-hat for seed survival in seed bank",adj=0)

plot.hist = hist(rhat, breaks = 25, plot = FALSE)
prob = sum(rhat < 1.05)/length(rhat)
text(max(plot.hist$mids), max(plot.hist$counts) * 0.95, 
     paste0("Summary of R-hat diagnostic"), pos = 2)
text(max(plot.hist$mids), max(plot.hist$counts) * 0.9, paste0("Percent of R-hat < 1.05: ", signif(prob, 2)), pos = 2)

# add portion of each chain that passes stationarity test
text(max(plot.hist$mids), max(plot.hist$counts) * 0.8, 
     paste0("Heidelberg-Welch diagnostic"), pos = 2)
text(1 * max(dt$mids), 0.75 * max(dt$counts), paste0("% passing (chain 1): ", p[1]), pos = 2)
text(1 * max(dt$mids), 0.7 * max(dt$counts), paste0("% passing (chain 2): ", p[2]), pos = 2)
text(1 * max(dt$mids), 0.65 * max(dt$counts), paste0("% passing (chain 3): ", p[3]), pos = 2)
dev.off()

# - Convergence diagnostics: germination ----


# - +Recover chains ----
all.chains = MCMCchains(mcmcSamples, params = c("mu0_g", "sigma0_g"), mcmc.list = TRUE)

# - +Rhat for all parameters ----
rhat = coda::gelman.diag(all.chains, confidence = 0.95)$psrf[, 1]

# - +Heidelberger and Welch convergence diagnostic ----
hd = coda::heidel.diag(all.chains)
p = c()
for (i in 1:3) {
  p[i] = sum(hd[[i]][, 1])/length(hd[[i]][, 1])
}
p = signif(p, 2)
dt = hist(coda::gelman.diag(all.chains, confidence = 0.95)$psrf[, 1], breaks = 25, 
  plot = FALSE)

# - +Plot Rhat and HW ----

par(mfrow = c(1, 1))
jpeg(filename = paste0(outputDirectory, "rhat-hw-seedBagExperiment-germination", ".jpeg"), 
  quality = 75)

# plot distribution of Rhat values
hist(rhat, col = "black", border = "white", breaks = 25, 
     main = NULL,
     xlab ="R-hat values")
title(main = "B. Distribution of R-hat for seed germination",adj=0)

plot.hist = hist(rhat, breaks = 25, plot = FALSE)
prob = sum(rhat < 1.05)/length(rhat)
text(max(plot.hist$mids), max(plot.hist$counts) * 0.95, 
     paste0("Summary of R-hat diagnostic"), pos = 2)
text(max(plot.hist$mids), max(plot.hist$counts) * 0.9, paste0("Percent of R-hat < 1.05: ", 
  signif(prob, 2)), pos = 2)

# add portion of each chain that passes stationarity test
text(max(plot.hist$mids), max(plot.hist$counts) * 0.8, 
     paste0("Heidelberg-Welch diagnostic"), pos = 2)
text(1 * max(dt$mids), 0.75 * max(dt$counts), paste0("% passing (chain 1): ", p[1]), pos = 2)
text(1 * max(dt$mids), 0.7 * max(dt$counts), paste0("% passing (chain 2): ", p[2]), pos = 2)
text(1 * max(dt$mids), 0.65 * max(dt$counts), paste0("% passing (chain 3): ", p[3]), pos = 2)
dev.off()



# - Convergence diagnostics: s0 ----



# - +R-hat ----
all.chains = MCMCchains(mcmcSamples, params = c("mu0_s0", "sigma0_s0"), 
  mcmc.list = TRUE)
# R-hat for all parameters
rhat = coda::gelman.diag(all.chains, confidence = 0.95)$psrf[, 1]

# - +Heidelberger and Welch convergence diagnostic ----
hd = coda::heidel.diag(all.chains)
p = c()
for (i in 1:3) {
  p[i] = sum(hd[[i]][, 1])/length(hd[[i]][, 1])
}
p = signif(p, 2)
dt = hist(coda::gelman.diag(all.chains, confidence = 0.95)$psrf[, 1], breaks = 25, 
  plot = FALSE)

# - +Plot Rhat and HW ----
par(mfrow = c(1, 1))
jpeg(filename = paste0(outputDirectory, "rhat-hw-seedBagExperiment-s0", ".jpeg"), quality = 75)
# plot distribution of Rhat values
hist(rhat, col = "black", border = "white", breaks = 25, 
     main = NULL,
     xlab ="R-hat values")
title(main = "C. Distribution of R-hat for seed survival, s0",adj=0)

plot.hist = hist(rhat, breaks = 25, plot = FALSE)
prob = sum(rhat < 1.05)/length(rhat)
text(max(plot.hist$mids), max(plot.hist$counts) * 0.95, 
     paste0("Summary of R-hat diagnostic"), pos = 2)
text(max(plot.hist$mids), max(plot.hist$counts) * 0.9, paste0("Percent of R-hat < 1.05: ", 
  signif(prob, 2)), pos = 2)
# add portion of each chain that passes stationarity test
text(max(plot.hist$mids), max(plot.hist$counts) * 0.8, 
     paste0("Heidelberg-Welch diagnostic"), pos = 2)
text(1 * max(dt$mids), 0.75 * max(dt$counts), paste0("% passing (chain 1): ", p[1]), pos = 2)
text(1 * max(dt$mids), 0.7 * max(dt$counts), paste0("% passing (chain 2): ", p[2]), pos = 2)
text(1 * max(dt$mids), 0.65 * max(dt$counts), paste0("% passing (chain 3): ", p[3]), pos = 2)
dev.off()
