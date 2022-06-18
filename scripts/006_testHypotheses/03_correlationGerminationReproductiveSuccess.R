# -------------------------------------------------------------------
# Analysis of correlation between germination and RS
# -------------------------------------------------------------------
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(rjags) # jags interface
library(MCMCvis)
library(tidyverse)
library(reshape2)
library(HDInterval)
library(bayesplot)
library(rethinking)

# -------------------------------------------------------------------
# Functions for use when analyzing data
# -------------------------------------------------------------------
temporal_variance <- function(x,fun=var){
  apply(x,1,fun)
}

cols_fun <- function(x,fun=var){
  apply(x,2,fun)
}

# geometric var
gsd <- function(x){
  y <- exp(sd(log(x)))
  return(y)
}

# sample based calculation of gsd
gsd.am <- function(x){
  x=x+.5
  n = length(x[!is.na(x)])
  mu = exp(mean(log(x),na.rm=TRUE))
  y <- exp(sqrt(sum((log(x/mu))^2,na.rm=TRUE)/(n-1)))
  return(y)
}

f<-function(x="param"){
  chain<-MCMCchains(zc,params = x)
  BCI <- t(apply(p,2,FUN = function(x) quantile(x, c(.025, .5, .975))))
  return(BCI)
}

posterior.mode = function(x){
 if(!is.na(x[1])){ x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]}else if(is.na(x[1])){
    modeParam <- NA
  }
  return(modeParam)
}

# -------------------------------------------------------------------
# Read in samples from posterior distributions
# -------------------------------------------------------------------

g1 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/g1-population-level.RDS")


g1 <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/g1-population-level.RDS")
sigma <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/sigma-population-year-level-mat.RDS")
fec <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/combinedF-population-year-level-mat.RDS")
phi <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/phi-population-year-level-mat.RDS")
rs <- readRDS("/Users/Gregor/Dropbox/clarkia-bet-hedging/outputs/005_calculatePopulationModelParameters/reproductiveSuccess-population-year-level-mat.RDS")

# get mode from full posterior and calculate var in RS based on modes
sigma.mode <- lapply(sigma,apply,2,posterior.mode)
fec.mode <- lapply(fec,apply,2,posterior.mode)
phi.mode <- lapply(phi,apply,2,posterior.mode)
rs.mode <- lapply(rs,apply,2,posterior.mode)



# ---
# - Site names by position ----
# ---
siteAbiotic <- read.csv("data/siteAbiotic.csv",header=TRUE)

position<-siteAbiotic %>% 
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteIndex <- order(position$easting,decreasing=FALSE)
siteNames = unique(position$site)[siteIndex]

siteNames <- unique(position$site)

# -------------------------------------------------------------------
# Compare calculation of posterior by mode of components vs. mode of RS
# -------------------------------------------------------------------

# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# for(i in 1:20){
#   obj = sigma.mode
#   index=grep(paste0("\\[",i,","),names(obj))
#   plot(sigma.mode[index]*fec.mode[index]*phi.mode[index],rs.mode[index],pch=19);
#   abline(a=0,b=1,col='red');text(sigma.mode[index]*fec.mode[index]*phi.mode[index],rs.mode[index],index)
# }

# -------------------------------------------------------------------
# Compare calculation of geometric SD in reproductive success
# -------------------------------------------------------------------
par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)
plot(NA,NA,xlim=c(0,10),ylim=c(0,10),pch=19);

for(i in 1:20){
  obj = sigma.mode
  index=grep(paste0("\\[",i,","),names(obj))
  text(gsd.am(sigma.mode[index]*fec.mode[index]*phi.mode[index]),gsd.am(rs.mode[index]),i);
  abline(a=0,b=1,col='red')
}

# -------------------------------------------------------------------
# Calculate geometric SD in reproductive success
# -------------------------------------------------------------------

df.list = list()
for(i in 1:20){
  obj = rs[[i]]
  tmp.df=apply(obj,1,gsd.am)
  df.list[[i]] = tmp.df
}

gsdSummary=do.call(cbind,df.list)

# -------------------------------------------------------------------
# Analyze correlation of germination and GSD per-capita RS
# -------------------------------------------------------------------

probability.g1 = g1

# - +Summary statistics on germination ----
# - ++Calculate the 95% credible interval ----
CI.g1 <- apply(probability.g1,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# - ++Calculate the 95% highest posterior density interval ----
HPDI.g1 <- apply(probability.g1,2,FUN = function(x) hdi(x, .68))
# - ++Calculate the posterior mode ----
mode.g1 <- apply(probability.g1,2, FUN = posterior.mode)

# - ++Construct data frame ----
# use highest posterior density interval and mode
g1PosteriorSummary <- data.frame(cbind(t(HPDI.g1),mode.g1))
names(g1PosteriorSummary) <- c("lo.g1","hi.g1","mode.g1")

# - +Summary statistics on variability in RS ----
# - ++Calculate the 95% credible interval ----
CI.rs <- apply(gsdSummary,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# - ++Calculate the 95% highest posterior density interval ----
HPDI.rs <- apply(gsdSummary,2,FUN = function(x) hdi(x, .68))
# - ++Calculate the posterior mode ----
mode.rs <- apply(gsdSummary,2, FUN = posterior.mode)

# - ++Construct data frame ----
rsPosteriorSummary<-data.frame(cbind(t(HPDI.rs),mode.rs))
names(rsPosteriorSummary) <- c("lo.rs","hi.rs","mode.rs")

# create empty vector for the correlation
posterior.correlation<-c()

n.iter = dim(gsdSummary)[1]
# calculate correlation for each draw from the posterior
for(i in 1:n.iter){
  posterior.correlation[i]<-cor(probability.g1[i,],gsdSummary[i,])
}

# - +Summary statistics on correlation ----
# - ++Calculate the 95% credible interval ----
CI.correlation <- quantile(posterior.correlation, c(.025, .5, .975))
# - ++Calculate the 95% highest posterior density interval ----
HPDI.correlation <- hdi(posterior.correlation, .95)
# - ++Calculate the 50% highest posterior density interval ----
HPDI.correlation.lo <- hdi(posterior.correlation, c(.5))
# - ++Calculate the posterior mode ----
mode.correlation <- posterior.mode(posterior.correlation)

# - ++Construct data frame ----
correlationPosteriorSummary<-data.frame(cbind(t(HPDI.correlation),mode.correlation))
names(correlationPosteriorSummary) <- c("lo.corr","hi.corr","mode.corr")

signif(correlationPosteriorSummary,3)

# - Manuscript figure ----

# - +set font sizes ----
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12
pt5 = 5/12

dev.off()

# tiff(filename=paste0("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/products/figures/correlationGerminationVariabilityRS.tif"),
#      height=3.2,width=3.2,units="in",res=800,compression="lzw",pointsize=12)
tiff(filename=paste0("~/Dropbox/clarkia-bet-hedging/outputs/correlationGerminationVariabilityRS.tif"),
     height=3.2,width=3.2,units="in",res=800,compression="lzw",pointsize=12)

par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(2,2,.7,0)+.1,mgp=c(3,.45,0))
plot(x = NA,
     y = NA,
     xlim=c(0,10),ylim=c(0,.4),
     pch=16, cex = 0.5,
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt8)


d.plot=data.frame(s=rsPosteriorSummary$mode.rs,
                  g=g1PosteriorSummary$mode.g1,
                  site=siteNames)

# d.plot[12,1:2] = d.plot[12,1:2]*c(1,1.04)
# d.plot[3,1:2] = d.plot[3,1:2]*c(1,1.05)
# d.plot[6,1:2] = c(4.1,.07)
# d.plot[11,1:2] = c(3.25,.05)
# d.plot[8,1:2] = c(4.6,.129)
# d.plot[5,1:2] = c(3.875,.13)

# segments(y0=.055,x0=3.25,y1=.07,x1=3.175,col='gray90')
# segments(y0=.1295,x0=4.45,y1=.1275,x1=4.35,col='gray90')
# segments(x0=3.98,x1=4,y0=.125,y1=.12,col='gray90')
# segments(x0=4.1,x1=3.8675,y0=.075,y1=.11,col='gray90')

segments(x0=rsPosteriorSummary$lo.rs,x1=rsPosteriorSummary$hi.rs,
         y0=g1PosteriorSummary$mode.g1, lwd=1)
segments(x0=rsPosteriorSummary$mode.rs,
         y0=g1PosteriorSummary$lo.g1, y1=g1PosteriorSummary$hi.g1,
         lwd=1)
points(rsPosteriorSummary$mode.rs,g1PosteriorSummary$mode.g1,
       pch=21,col='black',bg='white',cex=2.5)
points(rsPosteriorSummary$mode.rs[15],g1PosteriorSummary$mode.g1[15],
       pch=21,col='black',bg='white',cex=2.5)
# text(d.plot[,1:2],siteNames,cex=4/12)




axis(1, seq(0,10,by=2), padj = -.5,
     labels = seq(0,10,by=2), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(1, seq(0,10,by=1),labels=FALSE)
axis(2, seq(0,1,by=.1),
     labels = seq(0,1,by=.1), las = 1, line = 0, hadj= 1.2,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(2, seq(.05,1,by=.1),labels=FALSE)

mtext("Germination probability",
      side=2,line=1.5,adj=.5,col='black',cex=pt8)
mtext("Geometric standard deviation of per-capita reproductive success",
      side=1,line=1,adj=0,col='black',cex=pt8,at=-2)

box()
text(x=2.8,y=.4,
     paste0("Pearson's r=",round(mode.correlation,3)),
     cex=pt8)

mtext("A.", adj = 0, cex=pt10)

dev.off()

# 
# 
# pdf(
#   "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/correlation-germ-rs.pdf",
#   height = 8, width = 6)
# par(mar=c(4,4,2,1))
# par(fig=c(0,10,4,10)/10)
# # plot median of g1 vs. median of RS with CIs
# plot(x = NA,
#      y = NA,
#      xlim=c(0,12),ylim=c(0,.6),
#      pch=16, cex = 0.5,
#      xlab = "",
#      ylab = "",
#      xaxt= "n", yaxt="n",
#      cex.lab = 1, cex.axis = 1)
# 
# segments(x0=rsPosteriorSummary$lo.rs,x1=rsPosteriorSummary$hi.rs,
#          y0=g1PosteriorSummary$mode.g1, lwd=1)
# segments(x0=rsPosteriorSummary$mode.rs,
#          y0=g1PosteriorSummary$lo.g1, y1=g1PosteriorSummary$hi.g1,
#          lwd=1)
# points(rsPosteriorSummary$mode.rs,g1PosteriorSummary$mode.g1,
#        pch=21,col='black',bg='white',cex=1.25)
# 
# axis(1, seq(0,12,by=2),
#      labels = seq(0,12,by=2), las = 1, line = 0,
#      col = NA, col.ticks = 1, cex.axis = 1)
# axis(2, seq(0,1,by=.2),
#      labels = seq(0,1,by=.2), las = 1, line = 0,
#      col = NA, col.ticks = 1, cex.axis = 1)
# mtext("Germination probability",
#       side=2,line=2.5,adj=.5,col='black',cex=1)
# mtext("Geometric SD RS",
#       side=1,line=2,adj=.5,col='black',cex=1)
# 
# text(x=0,y=.58,
#      paste0("Pearson's r=",round(CI.correlation[2],2)),
#      cex=1,adj=c(0,0))
# #abline(a=0,b=1)
# 
# par(fig=c(0,10,0,4.5)/10)
# par(new=T)
# # plot posterior of correlation coefficient
# hist(posterior.correlation,breaks = 50, 
#      main = "", xlab = "", ylab='', xaxt='n',yaxt='n', 
#      xlim = c(-1, 1), ylim=c(0,2.5), border= "white",
#      freq = FALSE, col = "gray75", 
#      cex.lab = 1.25,cex.axis=1.5)
# 
# # as in Duskey dissertation
# segments(x0=HPDI.correlation[1],x1=HPDI.correlation[2],y0=2.4,lwd=1.25)
# segments(x0=HPDI.correlation.lo[1],x1=HPDI.correlation.lo[2],y0=2.4,lwd=2.5)
# points(x=median(posterior.correlation),y=2.4,pch=21,col='black',bg='white',cex=1.25)
# 
# axis(1, seq(-2,2,by=.2),
#      labels = seq(-2,2,by=.2), las = 1, line = 0,
#      col = NA, col.ticks = 1, cex.axis = 1)
# axis(2, seq(0,3,by=.4),
#      labels = seq(0,3,by=.4), las = 1, line = -.5,
#      col = NA, col.ticks = 1, cex.axis = 1)
# segments(x0=-1,x1=1,y0=-.1,lwd=2)
# segments(x0=-1.035,y0=0,y1=2.4,lwd=1.5)
# mtext("Density",
#       side=2,line=2.5,adj=.5,col='black',cex=1)
# mtext("Correlation of germination and  GSD RS",
#       side=1,line=2,adj=.5,col='black',cex=1)
# dev.off()

# -------------------------------------------------------------------
# Uncomment below to label points
# -------------------------------------------------------------------
# climate <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/climateData-2021.RDS")
# siteNames=climate %>% dplyr::filter(intenseDemography==1) %>% 
#   dplyr::select(site) %>% unique
# 
# df <- g1PosteriorSummary %>%
#   dplyr::bind_cols(site=siteNames) %>%
#   dplyr::left_join(rsPosteriorSummary %>%
#                      dplyr::bind_cols(site=siteNames),by="site")
# 
# library(ggrepel)
# 
# g1.plot <- ggplot(df,aes(x=mode.rs,y=mode.g1,label=site)) +
#   geom_point() +
#   geom_text_repel(size=3,color="black") +
#   # annotate("text", label =  paste0("Pearson's r=",round(CI.correlation[1],2)), x = 2.5, y = .29, size = 4) +
#   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
#   # scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
#   # scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
#   xlab("Geometric SD of reproductive success") +
#   ylab("Mean germination probability [P(G)]") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/correlation-germ-rs-labeled.pdf",
#        plot=g1.plot,width=4,height=4)


# -------------------------------------------------------------------
# Repeat with low fitness years set to 0
# -------------------------------------------------------------------
# 
# names<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>%
#   dplyr::select(site)
# siteNames = unique(names$site)

siteIndex <- data.frame(site=siteNames,siteIndex=1:20)
yearIndex <- data.frame(year=2006:2020,yearIndex=1:15)

lowFitnessYears <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/outputs/007_hypothesisTesting/objects/lowFitnessYears.RDS")

lowFitnessYears=lowFitnessYears %>% 
  dplyr::left_join(yearIndex) %>% 
  dplyr::left_join(siteIndex)

# set years without any plants at the population-level to zero fitness
df.list = list()
for(i in 1:20){
  obj = rs
  zeroYears=lowFitnessYears[lowFitnessYears$siteIndex==i,]$yearIndex
  index=grep(paste0("\\[",i,","),colnames(obj))
  obj=obj[,index]
  obj[,zeroYears] = 0
  tmp.df=apply(obj,1,gsd.am)
  df.list[[i]] = tmp.df
}

gsdSummary=do.call(cbind,df.list)

n.iter=dim(gsdSummary)[1]

probability.g1 = g1

# - +Summary statistics on germination ----
# - ++Calculate the 95% credible interval ----
CI.g1 <- apply(probability.g1,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# - ++Calculate the 95% highest posterior density interval ----
HPDI.g1 <- apply(probability.g1,2,FUN = function(x) hdi(x, .68))
# - ++Calculate the posterior mode ----
mode.g1 <- apply(probability.g1,2, FUN = posterior.mode)

# - ++Construct data frame ----
# use highest posterior density interval and mode
g1PosteriorSummary <- data.frame(cbind(t(HPDI.g1),mode.g1))
names(g1PosteriorSummary) <- c("lo.g1","hi.g1","mode.g1")

# - +Summary statistics on variability in RS ----
# - ++Calculate the 95% credible interval ----
CI.rs <- apply(gsdSummary,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# - ++Calculate the 95% highest posterior density interval ----
HPDI.rs <- apply(gsdSummary,2,FUN = function(x) hdi(x, .68))
# - ++Calculate the posterior mode ----
mode.rs <- apply(gsdSummary,2, FUN = posterior.mode)

# - ++Construct data frame ----
rsPosteriorSummary<-data.frame(cbind(t(HPDI.rs),mode.rs))
names(rsPosteriorSummary) <- c("lo.rs","hi.rs","mode.rs")

# create empty vector for the correlation
posterior.correlation<-c()

n.iter = dim(rs)[1]
# calculate correlation for each draw from the posterior
for(i in 1:n.iter){
  posterior.correlation[i]<-cor(probability.g1[i,],gsdSummary[i,])
}

# - +Summary statistics on correlation ----
# - ++Calculate the 95% credible interval ----
CI.correlation <- quantile(posterior.correlation, c(.025, .5, .975))
# - ++Calculate the 95% highest posterior density interval ----
HPDI.correlation <- hdi(posterior.correlation, .95)
# - ++Calculate the 50% highest posterior density interval ----
HPDI.correlation.lo <- hdi(posterior.correlation, c(.5))
# - ++Calculate the posterior mode ----
mode.correlation <- posterior.mode(posterior.correlation)

# - ++Construct data frame ----
correlationPosteriorSummary<-data.frame(cbind(t(HPDI.correlation),mode.correlation))
names(correlationPosteriorSummary) <- c("lo.corr","hi.corr","mode.corr")

signif(correlationPosteriorSummary,3)


# - Manuscript figure ----

# - +set font sizes ----
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12
pt5 = 5/12

dev.off()

# tiff(filename=paste0("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/products/figures/correlationGerminationVariabilityRS-lowFitness.tif"),
#      height=3.2,width=3.2,units="in",res=800,compression="lzw",pointsize=12)
tiff(filename=paste0("~/Desktop/figures/correlationGerminationVariabilityRS-lowFitness.tif"),
     height=3.2,width=3.2,units="in",res=800,compression="lzw",pointsize=12)

par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(2,2,.7,0)+.1,mgp=c(3,.45,0))
plot(x = NA,
     y = NA,
     xlim=c(2,10),ylim=c(0,.4),
     pch=16, cex = 0.5,
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt8)

segments(y0=.075,y1=.1,x0=4.1,x1=4.12,col='gray90')
segments(y0=.075,y1=.08,x0=2.85,x1=3.2,col='gray90')
segments(y0=.065,y1=.09,x0=3.9,x1=3.8,col='gray90')
segments(y0=.0725,y1=.09,x0=7.45,x1=7.3,col='gray90')

segments(x0=rsPosteriorSummary$lo.rs,x1=rsPosteriorSummary$hi.rs,
         y0=g1PosteriorSummary$mode.g1, lwd=1)
segments(x0=rsPosteriorSummary$mode.rs,
         y0=g1PosteriorSummary$lo.g1, y1=g1PosteriorSummary$hi.g1,
         lwd=1)
points(rsPosteriorSummary$mode.rs,g1PosteriorSummary$mode.g1,
       pch=21,col='black',bg='white',cex=2.5)
points(rsPosteriorSummary$mode.rs[15],g1PosteriorSummary$mode.g1[15],
       pch=21,col='black',bg='white',cex=2.5)
points(rsPosteriorSummary$mode.rs[6],g1PosteriorSummary$mode.g1[6],
       pch=21,col='black',bg='white',cex=2.5)

d.plot=data.frame(s=rsPosteriorSummary$mode.rs,
                  g=g1PosteriorSummary$mode.g1,
                  site=siteNames)

d.plot[3,1:2] = d.plot[3,1:2]*c(.9875,1.02)
d.plot[5,1:2] = c(4.169,.07)
d.plot[6,1:2] = d.plot[6,1:2]*c(1,1)
d.plot[9,1:2] = d.plot[9,1:2]*c(1,1.02)
d.plot[12,1:2] = d.plot[12,1:2]*c(1.05,.65)
d.plot[13,1:2] = d.plot[13,1:2]*c(1.05,.65)
d.plot[15,1:2] = d.plot[15,1:2]*c(.85,.8)
d.plot[17,1:2] = d.plot[17,1:2]*c(.9875,.9775)
text(d.plot[,1:2],siteNames,cex=4/12)

axis(1, seq(0,10,by=2), padj = -.5,
     labels = seq(0,10,by=2), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(1, seq(0,10,by=1),labels=FALSE)
axis(2, seq(0,1,by=.1),
     labels = seq(0,1,by=.1), las = 1, line = 0, hadj= 1.2,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(2, seq(.05,1,by=.1),labels=FALSE)

mtext("Germination probability",
      side=2,line=1.5,adj=.5,col='black',cex=pt8)
mtext("Geometric standard deviation of per-capita reproductive success",
      side=1,line=1,adj=-.39,col='black',cex=pt8,at=-3.5)

box()
text(x=3.25,y=.4,
     paste0("Pearson's r=",round(mode.correlation,3)),
     cex=pt8)

mtext("B.", adj = 0, cex=pt10)


dev.off()

# 
# pdf(
#   "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/correlation-germ-rs-lowfitness.pdf",
#   height = 8, width = 6)
# par(mar=c(4,4,2,1))
# par(fig=c(0,10,4,10)/10)
# # plot median of g1 vs. median of RS with CIs
# plot(x = NA,
#      y = NA,
#      xlim=c(0,12),ylim=c(0,.6),
#      pch=16, cex = 0.5,
#      xlab = "",
#      ylab = "",
#      xaxt= "n", yaxt="n",
#      cex.lab = 1, cex.axis = 1)
# 
# segments(x0=rsPosteriorSummary$lo.rs,x1=rsPosteriorSummary$hi.rs,
#          y0=g1PosteriorSummary$mode.g1, lwd=1)
# segments(x0=rsPosteriorSummary$mode.rs,
#          y0=g1PosteriorSummary$lo.g1, y1=g1PosteriorSummary$hi.g1,
#          lwd=1)
# points(rsPosteriorSummary$mode.rs,g1PosteriorSummary$mode.g1,
#        pch=21,col='black',bg='white',cex=1.25)
# 
# axis(1, seq(0,12,by=2),
#      labels = seq(0,12,by=2), las = 1, line = 0,
#      col = NA, col.ticks = 1, cex.axis = 1)
# axis(2, seq(0,1,by=.2),
#      labels = seq(0,1,by=.2), las = 1, line = 0,
#      col = NA, col.ticks = 1, cex.axis = 1)
# mtext("Germination probability",
#       side=2,line=2.5,adj=.5,col='black',cex=1)
# mtext("Geometric SD RS",
#       side=1,line=2,adj=.5,col='black',cex=1)
# 
# text(x=0,y=.58,
#      paste0("Pearson's r=",round(CI.correlation[2],2)),
#      cex=1,adj=c(0,0))
# #abline(a=0,b=1)
# 
# par(fig=c(0,10,0,4.5)/10)
# par(new=T)
# # plot posterior of correlation coefficient
# hist(posterior.correlation,breaks = 50, 
#      main = "", xlab = "", ylab='', xaxt='n',yaxt='n', 
#      xlim = c(-1, 1), ylim=c(0,2.5), border='white',
#      freq = FALSE, col = "gray75", 
#      cex.lab = 1.25,cex.axis=1.5)
# 
# # as in Duskey dissertation
# segments(x0=HPDI.correlation[1],x1=HPDI.correlation[2],y0=2.4,lwd=1.25)
# #segments(x0=CI.correlation[2],y0=2.3,y1=2.5,lwd=1.5)
# segments(x0=HPDI.correlation.lo[1],x1=HPDI.correlation.lo[2],y0=2.4,lwd=2.5)
# points(x=median(posterior.correlation),y=2.4,pch=21,col='black',bg='white',cex=1.25)
# #segments(x0=CI.correlation[2],y0=2.4,y1=0,lwd=2,lty='dotted')
# 
# axis(1, seq(-2,2,by=.2),
#      labels = seq(-2,2,by=.2), las = 1, line = 0,
#      col = NA, col.ticks = 1, cex.axis = 1)
# axis(2, seq(0,3,by=.4),
#      labels = seq(0,3,by=.4), las = 1, line = -.5,
#      col = NA, col.ticks = 1, cex.axis = 1)
# segments(x0=-1,x1=1,y0=-.1,lwd=2)
# segments(x0=-1.035,y0=0,y1=2.4,lwd=1.5)
# mtext("Density",
#       side=2,line=2.5,adj=.5,col='black',cex=1)
# mtext("Correlation of germination and  GSD RS",
#       side=1,line=2,adj=.5,col='black',cex=1)
# dev.off()

# -------------------------------------------------------------------
# Uncomment below to label points
# -------------------------------------------------------------------
# 
# df <- g1PosteriorSummary %>%
#   dplyr::bind_cols(site=siteNames) %>%
#   dplyr::left_join(rsPosteriorSummary %>%
#                      dplyr::bind_cols(site=siteNames),by="site")
# 
# library(ggrepel)
# 
# g1 <- ggplot(df,aes(x=med.rs,y=med.g1,label=site)) +
#   geom_point() +
#   geom_text_repel(size=3,color="black") +
#  # annotate("text", label =  paste0("Pearson's r=",round(CI.correlation[1],2)), x = 2.5, y = .29, size = 4) +
#  theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
#  # scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
#  # scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
#   xlab("Geometric SD of reproductive success") +
#   ylab("Mean germination probability [P(G)]") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/correlation-germ-rs-lowfitness-labeled.pdf",
#        plot=g1,width=4,height=4)

# compare RS calculations
# 
# rs <- readRDS("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/modelAnalysis/rsPosterior.RDS")
# lowFitnessYears <- readRDS("/Users/Gregor/Dropbox/clarkiaSeedBanks/scriptsAnalysis/output/lowFitnessYearsPlots.RDS")
# lowFitnessYears=lowFitnessYears %>% dplyr::left_join(yearIndex) %>% dplyr::left_join(siteIndex)
# 
# df.list = list()
# for(i in 1:20){
#   obj = rs
#   zeroYears=lowFitnessYears[lowFitnessYears$siteIndex==i,]$yearIndex
#   index=grep(paste0("\\[",i,","),colnames(obj))
#   obj=obj[,index]
#   obj[,zeroYears] = 0
#   tmp.df=apply(obj,1,gsd.am)
#   df.list[[i]] = tmp.df
# }
# 
# gsdSummaryLowFitness=do.call(cbind,df.list)
# 
# # calculate the 95% credible interval and HPDI for probability of RS
# CI.rsLowFitness <- apply(gsdSummaryLowFitness,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# 
# 
# df.list = list()
# for(i in 1:20){
#   obj = rs
#   index=grep(paste0("\\[",i,","),colnames(obj))
#   tmp.df=apply(obj[,index],1,gsd.am)
#   df.list[[i]] = tmp.df
# }
# 
# gsdSummary=do.call(cbind,df.list)
# 
# # calculate the 95% credible interval and HPDI for probability of RS
# CI.rs <- apply(gsdSummary,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# 
# 
# plot(CI.rsLowFitness[2,],CI.rs[2,],
#      xlim=c(2,12),ylim=c(2,12),
#      pch=19)
# abline(a=0,b=1)
# 
# delta = CI.rsLowFitness[2,]-CI.rs[2,]
# index=order(delta)
# delta=delta[index]
# plot(NA,NA,type='n',ylim=c(0,20),xlim=c(0,8),
#      axes=FALSE,frame=FALSE,
#      xlab="",ylab="")
# 
# for(i in 1:20){
#   tmp<-delta[i]
#  # segments(y0=tmp[1],y1=tmp[5],x0=y.pt[i])
#   #segments(y0=tmp[2],y1=tmp[4],x0=y.pt[i],lwd=3)
#   points(x=tmp,y=i,pch=21,bg='white')
# }
# axis(1,  seq(0,8,by=1), col.ticks = 1)
# axis(2, 1:20,
#      labels = (siteNames[index]), las = 2, 
#      col = NA, col.ticks = 1, cex.axis = 1)
# 
# position<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>% 
#   dplyr::select(site,easting,dominant.surface.rock.type) %>%
#   dplyr::mutate(easting=easting/1000)
# 
# delta = CI.rsLowFitness[2,]-CI.rs[2,]
# df=data.frame(delta,siteIndex) %>%
#   dplyr::left_join(position,by='site')
# 
# plot(df$easting,df$delta,type='n')
# text(x=df$easting,y=df$delta,labels=df$site)
# 
# rs.var.lowFitness = CI.rsLowFitness[2,]
# df=data.frame(rs.var.lowFitness,siteIndex) %>%
#   dplyr::left_join(position,by='site')
# 
# plot(df$easting,df$rs.var.lowFitness,type='n',ylim=c(0,20))
# text(x=df$easting,y=df$rs.var.lowFitness,labels=df$site)
# 
# rs.var = CI.rs[2,]
# df=data.frame(rs.var,siteIndex) %>%
#   dplyr::left_join(position,by='site')
# 
# plot(df$easting,df$rs.var,type='n',ylim=c(0,10))
# text(x=df$easting,y=df$rs.var,labels=df$site)
# # 
# # 
# # g1.blank <- ggplot(df,aes(x=var.rs,y=med,label=site)) +
# #  # geom_point() +
# #   # geom_text_repel(size=3,color="black") +
# #  # annotate("text", label =  paste0("Pearson's r=",round(CI.correlation[1],2)), x = 1.5, y = .29, size = 4) +
# #   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
# #   scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
# #   scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
# #   xlab("Geometric SD of reproductive success") +
# #   ylab("Mean germination probability [P(G)]") +
# #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# # 
# # ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-blank.pdf",
# #        plot=g1.blank,width=4,height=4)
# # 
# # 
# # # g1.point <- ggplot(df %>% dplyr::filter(site=="BG"),aes(x=var.rs,y=med,label=site)) +
# # #    geom_point() +
# # #   # geom_text_repel(size=3,color="black") +
# # #   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
# # #   scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
# # #   scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
# # #   xlab("Geometric SD of reproductive success") +
# # #   ylab("Mean germination probability [P(G)]") +
# # #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# # # 
# # # ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-point.pdf",
# # #        plot=g1.point,width=4,height=4)
# # 
# # matplot(c0[,2],c0[,3:12],
# #         xlab="Mean site germination probability",
# #         ylab="Annual mean fitness",
# #         main=expression(paste("cor(",sigma,",P(G))=-.75")),
# #         pch=16,col='gray')
# # plot(c0[1:nsites,1],c0[1:nsites,2],xlab="Variance in Fitness",ylab="Probability of Germination")
# # 
# # set.seed(11)
# # 
# # # COPULA
# # 
# # library(MASS)
# # 
# # nsites = 20
# # nyears = 10
# # 
# # mu = mu = rnorm(n=nsites,mean=25,sd=0)
# # 
# # sim <- function(rho=0,nsites=20){
# #   # Defines the  sequence for stochastic trials.
# #   reps=1000
# #   
# #   # a and b are hyperparameters of the gamma distribution 
# #   # that define both the expected value and variance.   
# #   a = 6
# #   b = 1.75
# #   
# #   # alpha and beta are hyperparameters of the beta distribution that define both the expected value and
# #   # variance.  
# #   alpha =  1
# #   beta =  1
# #   
# #   # Defines the temporal correlation between the two parameters.
# #   rho = rho
# #   
# #   # Generates standard multivariate normal data with correlation structure defined by rho.
# #   Z <- mvrnorm(n=reps,mu=c(0,0), 
# #                matrix(data=c(1,rho,rho,1),
# #                       nrow=2,ncol=2))
# #   
# #   # Apply the Normal CDF function to Z to obtain data that is uniform on the interval [0,1], but still correlated.
# #   U = pnorm(Z)
# #   
# #   # x is gamma distributed
# #   X = qgamma(U[,1],shape=a,rate=b) 
# #  # X = qunif(U[,1],0,.75) 
# #   
# #   # y is beta distributed
# #   #X <- cbind(X,qbeta(U[,2],shape1=alpha,shape2=beta) )
# #   X <- cbind(X,qunif(U[,2],0,.3) )
# #   
# #   # gamma marginal of multivariate X.
# #   # hist(X[,1])
# #   # beta marginal of multivariate X
# #   # hist(X[,2])
# #   
# #   # plot(X[,1],X[,2])
# #   
# #   nsites = nsites
# #   nyears = 10
# #   X[1:nsites,1:2]
# #   
# #   # generate samples for each site with a different variance
# #   d<-lapply(X[1:nsites,1],rnorm,n=nyears,mean=rnorm(n=nsites,mean=mu,sd=0))
# #   # check that the variances are appropriately sampled
# #   cbind(unlist(lapply(d, sd)),X[1:nsites,1])
# #   
# #   d<-data.frame(d)
# #   names(d) <- 1:20
# #   
# #   # check variances again
# #   cbind(apply(d,2, sd),X[1:nsites,1])
# #   
# #   dim(d)
# #   
# #   dt<-cbind(X[1:nsites,1:2],t(d))
# #   return(dt)
# # }
# # 
# # 
# # 
# # c0<-sim(rho=-.75,nsites=20)
# # 
# # df.sim <- data.frame(df$site,c0)
# # g1.hypothesis<-ggplot(df.sim,aes(x=X,y=V2)) +
# #   geom_point(color='gray') +
# #   # geom_text_repel(size=3,color="black") +
# #    annotate("text", label =  paste0("Pearson's r=",-.75), x = 6, y = .29, size = 4,color='gray') +
# #   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
# #   scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
# #   scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
# #   xlab("Geometric SD of reproductive success") +
# #   ylab("Mean germination probability [P(G)]") +
# #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# # 
# # ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-hypothesis.pdf",
# #        plot=g1.hypothesis,width=4,height=4)
# # 
# # 
# # g1.point<-ggplot(df.sim[1,],aes(x=X,y=V2)) +
# #   geom_segment(aes(x=X,y=0,xend=X,yend=V2),color='#E69F00') +
# #   geom_segment(aes(x=1,y=V2,xend=X,yend=V2),color='#CC79A7') +
# #   
# #   geom_point(color='gray',size=5) +
# #   # geom_text_repel(size=3,color="black") +
# #  # annotate("text", label =  paste0("Pearson's r=",-.75), x = 6, y = .29, size = 4,color='gray') +
# #   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
# #   scale_x_continuous(limits = c(.99,8), expand = c(0, 0), breaks = c(1, 3, 5, 7)) +
# #   scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
# #   xlab("Geometric SD of reproductive success") +
# #   ylab("Mean germination probability [P(G)]") +
# #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# # 
# # ggsave(filename=  "~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-point.pdf",
# #        plot=g1.point,width=4,height=4)
# # 
# # 
# # # # Germination
# # # # extract parameters for analysis
# # #  posterior.g1<-MCMCchains(belowground,params = "g1")
# # # 
# # # # calculate the 95% credible interval and HPDI for g1
# # # CI.g1 <- apply(posterior.g1,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# # # HPDI.g1 <- apply(posterior.g1,2,FUN = function(x) rethinking::HPDI(x, .95))
# # # 
# # # ## Reproductive success
# # # lapply(aboveground,)
# # # 
# # # d<-lapply(aboveground,temporal_variance,fun=gsd)
# # # reproductiveSuccess<-matrix(unlist(d), ncol = 20, byrow = FALSE)
# # # # # exclude problem sites
# # # # reproductiveSuccess<-reproductiveSuccess[,-probs]
# # # 
# # # # don't currently have the full distribution?
# # # CI.reproductiveSuccess <- apply(reproductiveSuccess,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# # # HPDI.reproductiveSuccess <- apply(reproductiveSuccess,2,FUN = function(x) hdi(x, .95))
# # # 
# # # g1PosteriorSummary<-data.frame(t(CI.g1))
# # # names(g1PosteriorSummary) <- c("lo.g1","med.g1","hi.g1")
# # # 
# # # ## need to finish this section
# # # reproductiveSuccessPosteriorSummary<-data.frame(t(CI.reproductiveSuccess))
# # # names(reproductiveSuccessPosteriorSummary) <- c("lo.rs","med.rs","hi.rs")
# # # 
# # # # calculate correlation for each draw from the posterior
# # # n.iter=3000
# # # posterior.correlation<-c()
# # # 
# # # for(i in 1:n.iter){
# # #   posterior.correlation[i]<-cor(posterior.g1[i,],reproductiveSuccess[i,])
# # # }
# # # 
# # # # calculate the 95% credible interval and HPDI for the correlation
# # # CI.correlation <- quantile(posterior.correlation, c(.025, .5, .975))
# # # HPDI.correlation <- hdi(posterior.correlation, .95)
# # # 
# # # 
# # #   
# # # 
# # # ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/germ_rs_correlation-labeled.pdf",
# # #        plot=g1,width=6,height=6)
# # # 
# # # 
# # # 
# # # hist(posterior.correlation,breaks = 50, main = "", xlab = "", xlim = c(-1, 1),
# # #      freq = FALSE, col = "azure1", cex.lab = 1.5,cex.axis=1.5)
# # # 
# # # title(xlab="Correlation of germination and geometric SD of fitness \n (Pearson's r)", line=4, cex.lab=1.5)
# # # 
# # # abline(v=CI.correlation[c(1,3)],lty='dashed',lwd='2')
# # # abline(v=CI.correlation[2],lty='solid',lwd='2',col='red')
# 
# 
# # -------------------------------------------------------------------
# # Code below was comparing calculation of variance in reproductive success
# # -------------------------------------------------------------------
# 
# f.rs = function(s,f,p){
#   tmp = s*f*p
#   df.list = list()
#   for(i in 1:20){
#     obj = tmp
#     index=grep(paste0("\\[",i,","),names(obj))
#     tmp.df=gsd.am(obj[index])
#     df.list[[i]] = tmp.df
#   }
#   unlist(df.list)
# }
# 
# # calculate rs. mode based on modes of components
# rs.mode = sigma.mode*fec.mode*phi.mode
# 
# df.list = list()
# for(i in 1:20){
#   obj = rs.mode
#   index=grep(paste0("\\[",i,","),names(obj))
#   tmp.df=gsd.am(obj[index])
#   df.list[[i]] = tmp.df
# }
# 
# gsdSummary=unlist(df.list)
# 
# 
# for(i in 1:20){
#   obj = sigma.mode*fec.mode*phi.mode
#   index=grep(paste0("\\[",i,","),colnames(obj))
#   tmp.df=apply(obj[,index],1,gsd.am)
#   df.list[[i]] = tmp.df
# }
# 
# n.iter = dim(sigma)[1]
# gsdSummaryComponent=matrix(NA,ncol=20,nrow=n.iter)
# 
# for(i in 1:20){gsdSummaryComponent[,i]=df.list[[i]]}
# 
# 
# par(mfrow = c(4,5),
#     oma = c(5,4,0,0) + 0.1,
#     mar = c(0,0,1,1) + 0.1)
# for(i in 1:20){
#   hist(df.list[[i]],breaks=100);abline(v=gsdSummary[i],col='red',lwd=2)
# }
# 
# # get mode from resampling posterior and calculate var in RS based on modes
# n.iter=dim(sigma)[1]
# n.sub = 1000
# n.rep=100
# 
# # sig.samples=sample(n.iter,n.rep)
# # fec.samples=sample(n.iter,n.rep)
# # phi.samples=sample(n.iter,n.rep)
# 
# # repeat calculation from above using subsamples from the posterior
# rs.mat=matrix(NA,nrow=n.rep,ncol=20)
# for(i in 1:n.rep){
#   sig.samples=sample(1:n.iter,n.sub)
#   fec.samples=sample(1:n.iter,n.sub)
#   phi.samples=sample(1:n.iter,n.sub)
#   
#   sigma.mode.tmp = apply(sigma[sig.samples,],2,posterior.mode);
#   fec.mode.tmp = apply(fec[fec.samples,],2,posterior.mode);
#   phi.mode.tmp = apply(phi[phi.samples,],2,posterior.mode);
#   
#   rs.mat[i,]=f.rs(sigma.mode.tmp,fec.mode.tmp,phi.mode.tmp)
# }
# 
# # variance in reproductive success based on modes alone is higher than
# # variance in reproductive success based on full distribution
# # this is likely because when using the full distribution there is some chance
# # of capturing high probability years
# # issue in site 10, 12, 20, 17, 18
# plot(rs.mat)
# 
# # par(mfrow=c(4,5),oma=c(0,0,0,0),mar=c(0,0,0,0))
# # for(i in 1:20){
# # plot(rs.mat[,i],ylim=c(0,10));abline(h=gsdSummary[i],col="red",lwd=2)
# # }
# # 
# # for(i in 1:20){hist(rs.mat[,i],breaks=25,main="");abline(v=gsdSummary[i],col="red",lwd=2)}
# # n.iter=dim(gsdSummary)[1]