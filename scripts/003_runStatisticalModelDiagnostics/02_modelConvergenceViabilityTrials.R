####
####
# Script to evaluate convergence
# Viability model
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
mcmcSamples <- readRDS(mcmcSampleFiles[[grep("viabilityPosteriorSamples.RDS", mcmcSampleFiles)]])

# - Source function ----
# function produces jpegs of trace plots
source("scripts/003_runStatisticalModelDiagnostics/00_tracePlotFunction.R")

# - Convergence diagnostics ----

# - +Graphical checks on trace plots ----
outputDirectory = "outputs/003_runStatisticalModelDiagnostics/01_tracePlots/2_viability/"

f(x = "mu0_g", model = "germinationTrials")
f(x = "sigma0_g", model = "germinationTrials")
f(x = "mu0_v", model = "viabilityTrials")
f(x = "sigma0_v", model = "viabilityTrials")

# - +Recover chains ----
all.chains = MCMCchains(mcmcSamples, params = c("mu0_g", "sigma0_g", "mu0_v","sigma0_v"), 
                        mcmc.list = TRUE)

# - +R-hat for all parameters ----
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

par(mfrow = c(1, 1))
jpeg(filename = paste0(outputDirectory, "rhat-hw-viabilityTrials", ".jpeg"), quality = 75)
# plot distribution of Rhat values
hist(rhat, col = "black", border = "white", breaks = 25, 
     main = NULL,
     xlab ="R-hat values",cex.axis=1.2,cex.lab=1.4)
title(main = "D. Distribution of R-hat for seed viability",adj=0,cex.main=1.4)

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
