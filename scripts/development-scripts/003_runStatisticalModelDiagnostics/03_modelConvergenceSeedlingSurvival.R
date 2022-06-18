####
####
# Script to evaluate convergence
# Seedling survival to fruiting model
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
mcmcSamples <- readRDS(mcmcSampleFiles[[grep("seedlingSurvivalSamplesSimplePlots.RDS", mcmcSampleFiles)]])

# - Source function ----
# function produces jpegs of trace plots
source("scripts/003_runStatisticalModelDiagnostics/00_tracePlotFunction.R")

# - Convergence diagnostics: seed survival ----

# - +Graphical checks on trace plots ----
outputDirectory = "outputs/003_runStatisticalModelDiagnostics/01_tracePlots/3_seedlingSurvival/"

f(x = "mu0", model = "seedlingSurvivalPlots")
f(x = "sigma0", model = "seedlingSurvivalPlots")
f(x = "sigma_transect", model = "seedlingSurvivalPlots")

# - +Recover chains ----
all.chains = MCMCchains(mcmcSamples, params = c("mu0", "sigma0", "sigma_transect"), mcmc.list = TRUE)

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
outputDirectory = "outputs/003_runStatisticalModelDiagnostics/02_rhatDistribution/"

par(mfrow=c(1,1))
jpeg(filename = paste0(outputDirectory, "rhat-hw-seedlingSurvivalPlots", ".jpeg"), 
  quality = 75)
# plot distribution of Rhat values
hist(rhat, col = "black", border = "white", breaks = 25, main = "Distribution of R-hat; seedling survival")
plot.hist = hist(rhat, breaks = 25, plot = FALSE)
prob = sum(rhat < 1.05)/length(rhat)
text(max(plot.hist$mids), max(plot.hist$counts) * 0.9, paste0("Percent of R-hat < 1.05: ", 
  signif(prob, 2)), pos = 2)
# add portion of each chain that passes stationarity test
text(0.999 * max(dt$mids), 0.85 * max(dt$counts), paste0("% passing (chain 1): ", 
  p[1]), pos = 2)
text(0.999 * max(dt$mids), 0.8 * max(dt$counts), paste0("% passing (chain 2): ", 
  p[2]), pos = 2)
text(0.999 * max(dt$mids), 0.75 * max(dt$counts), paste0("% passing (chain 3): ", 
  p[3]), pos = 2)
dev.off()
