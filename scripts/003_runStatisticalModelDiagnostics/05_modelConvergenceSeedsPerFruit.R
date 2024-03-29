####
####
# Script to evaluate convergence
# Seeds per fruit model
####
####

# - Environment ----
# clear environment but keep directories for data, models, and output files
rm(list=(ls())) # if using in source(script), include variables to keep
options(stringsAsFactors = FALSE)

outMCMCDirectory = "outputs/002_fitStatisticalModels/mcmcSamples/"

# - Libraries ----
library(MCMCvis)
library(coda)
library(magrittr)

# - Read in MCMC samples ----
mcmcSampleFiles <- paste0(outMCMCDirectory, list.files(outMCMCDirectory))
mcmcSamples <- readRDS(mcmcSampleFiles[[grep("seedsPerFruitPosteriorSamples.RDS", mcmcSampleFiles)]])

# - Source function ----
# function produces jpegs of trace plots
source("scripts/003_runStatisticalModelDiagnostics/00_tracePlotFunction.R")


# - Convergence diagnostics: seed survival ----

# - +Graphical checks on trace plots ----
outputDirectory = "outputs/003_runStatisticalModelDiagnostics/01_tracePlots/5_seedsPerFruit/"

f(x = "nu_seeds", model = "seedsUndamaged")
f(x = "sigma0_seeds", model = "seedsUndamaged")
f(x = "sigma_overdisp_seeds", model = "seedsUndamaged")

# print trace plots for each parameter
f(x = "mu0", model = "propDamage")
f(x = "sigma0", model = "propDamage")

# - +Recover chains for undamaged seeds ----
all.chains = MCMCchains(mcmcSamples, params = c("nu_seeds", "sigma0_seeds", "sigma_overdisp_seeds"), mcmc.list = TRUE)

# - +R-hat ----
par(mfrow = c(1, 1))
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

# plot distribution of Rhat values
outputDirectory = "outputs/003_runStatisticalModelDiagnostics/02_rhatDistribution/"

jpeg(filename = paste0(outputDirectory, "rhat-heidelberg-seedsUndamaged", ".jpeg"), 
  quality = 75)
hist(rhat, col = "black", border = "white", breaks = 25, 
     main = NULL,
     xlab ="R-hat values",cex.axis=1.2,cex.lab=1.4)
title(main = "E. Distribution of R-hat for seeds/undamaged fruit",adj=0,cex.main=1.4)
plot.hist = hist(rhat, breaks = 25, plot = FALSE)
prob = sum(rhat < 1.05)/length(rhat)
text(max(plot.hist$mids), max(plot.hist$counts) * 0.95, 
     paste0("Summary of R-hat diagnostic"), pos = 2,cex=1.4)
text(max(plot.hist$mids), max(plot.hist$counts) * 0.9, paste0("Percent of R-hat < 1.05: ", 
  signif(prob, 2)), pos = 2,cex=1.4)
# add portion of each chain that passes stationarity test
text(max(plot.hist$mids), max(plot.hist$counts) * 0.8, 
     paste0("Heidelberg-Welch diagnostic"), pos = 2,cex=1.4)
text(1 * max(dt$mids), 0.75 * max(dt$counts), paste0("% passing (chain 1): ", 
  p[1]), pos = 2,cex=1.4)
text(1 * max(dt$mids), 0.7 * max(dt$counts), paste0("% passing (chain 2): ", 
  p[2]), pos = 2,cex=1.4)
text(1 * max(dt$mids), 0.65 * max(dt$counts), paste0("% passing (chain 3): ", 
  p[3]), pos = 2,cex=1.4)
dev.off()


# - +Recover chains for proportion of seeds damaged ----
all.chains = MCMCchains(mcmcSamples, params = c("mu0", "sigma0"), mcmc.list = TRUE)

# - +R-hat ----
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

jpeg(filename = paste0(outputDirectory, "rhat-heidelberg-seedsPropDamaged", ".jpeg"), 
     quality = 75)
# plot distribution of Rhat values
hist(rhat, col = "black", border = "white", breaks = 25, 
     main = NULL,
     xlab ="R-hat values",cex.axis=1.2,cex.lab=1.4)
title(main = "F. Distribution of R-hat for proportion\ndamaged seeds",adj=0,cex.main=1.4)

plot.hist = hist(rhat, breaks = 25, plot = FALSE)
prob = sum(rhat < 1.05)/length(rhat)
text(max(plot.hist$mids), max(plot.hist$counts) * 0.95, 
     paste0("Summary of R-hat diagnostic"), pos = 2, cex = 1.4)
text(max(plot.hist$mids), max(plot.hist$counts) * 0.9, paste0("Percent of R-hat < 1.05: ", 
                                                              signif(prob, 2)), pos = 2, cex = 1.4)
# add portion of each chain that passes stationarity test
text(max(plot.hist$mids), max(plot.hist$counts) * 0.8, 
     paste0("Heidelberg-Welch diagnostic"), pos = 2, cex = 1.4)
text(1 * max(dt$mids), 0.75 * max(dt$counts), paste0("% passing (chain 1): ", 
                                                         p[1]), pos = 2, cex = 1.4)
text(1 * max(dt$mids), 0.7 * max(dt$counts), paste0("% passing (chain 2): ", 
                                                        p[2]), pos = 2, cex = 1.4)
text(1 * max(dt$mids), 0.65 * max(dt$counts), paste0("% passing (chain 3): ", 
                                                         p[3]), pos = 2, cex = 1.4)
dev.off()
