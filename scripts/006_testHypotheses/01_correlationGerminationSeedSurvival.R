# --
#  Script to analyze correlation between g1 and s2s3
# --

# ---
# - Set up environment ----
# ---

# - +remove unused objects ----
rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE,max.print=100000)

# - +source scripts ----
source("scripts/006_testHypotheses/00_utilityFunctions.R")

# - +load packages ----
library(MCMCvis)
library(tidyverse)
library(HDInterval)
library(bayesplot)

# ---
# - Site names by position ----
# ---
siteAbiotic <- read.csv("data/siteAbioticData.csv",header=TRUE)

position<-siteAbiotic %>%
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteIndex <- order(position$easting,decreasing=FALSE)
siteNames = unique(position$site)[siteIndex]

# ---
# - Read in germination and survival estimates ----
# ---

g1 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/g1-population-level.RDS")
s2 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s2-population-level.RDS")
s3 <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s3-population-level.RDS")

# ---
# - Histograms of posteriors ----
# ---
par(mfrow = c(4,5),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1    )

# in histograms, the following quantities are noted
# Median: red solid line
# Mode: purple solid line
# 95% percentile interval (.025, .975): orange dotted lines
# Highest posterior density interval (95%): blue dotted lines; calculated with HDinterval

# - +Germination ----

for(i in 1:20){
  tmp.i = siteIndex[i]
  hist(g1[,tmp.i],breaks=100,main='');
  abline(v=median(g1[,tmp.i]),lwd=2,col='red');
  abline(v=posterior.mode(g1[,tmp.i]),lwd=2,col='purple');
  abline(v=quantile(g1[,tmp.i],c(.025, .5, .975)),lwd=1,col='orange');
  abline(v=HDInterval::hdi(g1[,tmp.i],c(.95)),lwd=1,col='blue');
  mtext(siteNames[i],adj=0)
}

# - +Seed survival (s2*s3) ----

for(i in 1:20){
  tmp.i = siteIndex[i]
  hist(s2[,tmp.i]*s3[,tmp.i],breaks=100,main='');
  abline(v=median(s2[,tmp.i]*s3[,tmp.i]),lwd=2,col='red');
  abline(v=posterior.mode(s2[,tmp.i]*s3[,tmp.i]),lwd=2,col='purple');
  abline(v=quantile(s2[,tmp.i]*s3[,tmp.i],c(.025, .975)),lwd=1,col='orange');
  abline(v=HDInterval::hdi(s2[,tmp.i]*s3[,tmp.i],c(.95)),lwd=1,col='blue');
  mtext(siteNames[i],adj=0)
}

# in both cases, the mode & HPDI seem more appropriate as a summary,
# since the posteriors are bounded [0,1] and not necessarily symmetric

# ---
# - Analyze correlations ----
# ---

# - +Dimensions of matrix ----
n.iter = dim(g1)[1]

# - +Rename objects ----
probability.survival = s2*s3
probability.g1 = g1

# - +Summary statistics on germination ----
# - ++Calculate the 95% percentile interval ----
CI.g1 <- apply(probability.g1,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# - ++Calculate the 68% highest posterior density interval ----
HPDI.g1.68 <- apply(probability.g1,2,FUN = function(x) hdi(x, .68))
HPDI.g1.95 <- apply(probability.g1,2,FUN = function(x) hdi(x, .95))
# - ++Calculate the posterior mode ----
mode.g1 <- apply(probability.g1,2, FUN = posterior.mode)

# - ++Construct data frame ----
# use highest posterior density interval and mode
g1PosteriorSummary <- data.frame(cbind(t(HPDI.g1.68),t(HPDI.g1.95),mode.g1))
names(g1PosteriorSummary) <- c("lo.g1","hi.g1","lolo.g1","hihi.g1","mode.g1")

# - +Summary statistics on seed survival (s2*s3) ----
# - ++Calculate the 95% percentile interval ----
CI.survival <- apply(probability.survival,2,FUN = function(x) quantile(x, c(.025, .5, .975)))
# - ++Calculate the 68% highest posterior density interval ----
HPDI.survival.68 <- apply(probability.survival,2,FUN = function(x) hdi(x, .68))
HPDI.survival.95 <- apply(probability.survival,2,FUN = function(x) hdi(x, .95))
# - ++Calculate the posterior mode ----
mode.survival <- apply(probability.survival,2, FUN = posterior.mode)

# - ++Construct data frame ----
survivalPosteriorSummary<-data.frame(cbind(t(HPDI.survival.68),t(HPDI.survival.95),mode.survival))
names(survivalPosteriorSummary) <- c("lo.surv","hi.surv","lolo.surv","hihi.surv","mode.surv")

# - +Empty vector for the correlation ----
posterior.correlation<-c()

# - +Calculate correlation for each draw from the posterior ----
# calculation is done rowwise
for(i in 1:n.iter){
  posterior.correlation[i]<-cor(probability.g1[i,],probability.survival[i,])
}

# - +Summary statistics on correlation ----
# - ++Calculate the 95% credible interval ----
CI.correlation <- quantile(posterior.correlation, c(.025, .5, .975))
# - ++Calculate the 95% highest posterior density interval ----
HPDI.correlation <- hdi(posterior.correlation, .95)
# - ++Calculate the posterior mode ----
mode.correlation <- posterior.mode(posterior.correlation)

# - ++Construct data frame ----
correlationPosteriorSummary<-data.frame(cbind(t(HPDI.correlation),mode.correlation))
names(correlationPosteriorSummary) <- c("lo.corr","hi.corr","mode.corr")

signif(correlationPosteriorSummary,3)

# - Write out summaries ----


# df.summary<-data.frame(t(CI.survival))
# names(df.summary) <- c("lo.s2s3","med.s2s3","hi.s2s3")
# names<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>%
#   dplyr::select(site)
# siteNames = unique(names$site)
# df.summary=df.summary %>%
#   dplyr::bind_cols(site=siteNames)
# df.summary=cbind(site=df.summary[,4],df.summary[,1:3])
#
# write.csv(df.summary,"~/Downloads/s2s3-estimates.csv")


# - +set font sizes ----
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12
pt5 = 5/12

dev.off()

tiff(filename=paste0("products/figures/correlationGerminationSeedSurvival.tif"),
     height=3.2,width=3.2,units="in",res=800,compression="lzw",pointsize=12)

par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(2,2.2,.7,0)+.1,mgp=c(3,.45,0))
# plot mode of g1 vs. mode of s2*s3 with CIs
plot(x = NA,
     y = NA,
     xlim=c(.2,.8),ylim=c(0,.425),
     pch=16, cex = 0.5,
     xlab = "",
     ylab = "",
     xaxt= "n", yaxt="n",
     cex.lab = pt10, cex.axis = pt8)

# commented code would be to label points
# d.plot=data.frame(s=survivalPosteriorSummary$mode.surv,
#                   g=g1PosteriorSummary$mode.g1,
#                   site=siteNames)

# d.plot[13,1:2] = c(.61, .07)
# d.plot[4,1:2] = c(.65,.07)
# d.plot[3,1:2] = c(.689,.145)
# d.plot[8,1:2] = c(.78,.125)
# d.plot[14,1:2] = c(.73,.07)
# d.plot[10,1:2] = c(.68,.06)


segments(x0=survivalPosteriorSummary$lo.surv,x1=survivalPosteriorSummary$hi.surv,
         y0=g1PosteriorSummary$mode.g1, lwd=1)
segments(x0=survivalPosteriorSummary$mode.surv,
         y0=g1PosteriorSummary$lo.g1, y1=g1PosteriorSummary$hi.g1,
         lwd=1)

points(survivalPosteriorSummary$mode.surv,g1PosteriorSummary$mode.g1,
       pch=21,cex=1.5,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 1),lwd=0)
points(survivalPosteriorSummary$mode.surv,g1PosteriorSummary$mode.g1,
       pch=21,col='black',cex=1.5,
       bg=rgb(red = 1, green = 1, blue = 1, alpha = 0.5))

# commented code would be to label points
# text(d.plot[,1:2],siteNames,cex=4/12)


axis(1, seq(0,1,by=.1), padj = -.5,
     labels = seq(0,1,by=.1), line = 0,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(1, seq(.25,1,by=.1),labels=FALSE)
axis(2, seq(0,1,by=.1),
     labels = seq(0,1,by=.1), las = 1, line = 0, hadj= 1.2,
     col = NA, col.ticks = 1, cex.axis = pt8)
axis(2, seq(.05,1,by=.1),labels=FALSE)

mtext("Germination probability",
      side=2,line=1.5,adj=.5,col='black',cex=pt10)
mtext("Seed survival probability",
      side=1,line=1,adj=.5,col='black',cex=pt10)

box()
legend("topleft",bty='n',inset=c(-.05,0),
     paste0("Pearson's r=",round(correlationPosteriorSummary[3],3),
            " (",round(correlationPosteriorSummary[1],3),", ",
            round(correlationPosteriorSummary[2],3),")"),
     cex=pt7)

mtext("A.", adj = 0, cex=pt10)


dev.off()

# -------------------------------------------------------------------
# Commented code below is for labeling populations
# -------------------------------------------------------------------
#
# names<-read.csv(file="~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE) %>%
#   dplyr::select(site)
# df <- data.frame(names,survivalPosteriorSummary,g1PosteriorSummary)
# ciPosteriorSummary <- data.frame(t(CI.correlation))
#
#
# library(ggrepel)
#
# g1 <- ggplot(df,aes(x=med.surv,y=med.g1,label=site)) +
#   geom_point() +
#   geom_text_repel(size=3,color="black") +
#   theme_bw()+
#   xlab("Probability of seed survival [P(S)]") +
#   ylab("Mean germination probability [P(G)]") +
#   #annotate("text", label =  paste0("Pearson's r=",round(CI.correlation[2],2)), x = .15, y = .29, size = 4) +
#   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
#  # scale_x_continuous(limits = c(0,.75), expand = c(0, 0)) +
#  # scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#
#
# tiff(filename=paste0("/Users/Gregor/Dropbox/clarkiaSeedBanks/analysis/products/figures/correlationGerminationSeedSurvival-2.tif"),
#      height=4,width=4,units="in",res=800,compression="lzw",pointsize=12)
#
# # pdf(
# #  "~/Dropbox/clarkiaSeedBanks/products/figures/analysis/correlation-germ-surv.pdf",
# #  height = 8, width = 6)
#
# par(mar=c(3,2.5,0,1)+.1,oma=c(0,0,0,0))
# par(fig=c(0,10,4,10)/10)
# # plot median of g1 vs. median of s2*s3 with CIs
# plot(x = NA,
#      y = NA,
#      xlim=c(.2,.9),ylim=c(0,.55),
#      pch=16, cex = 0.5,
#      xlab = "",
#      ylab = "",
#      xaxt= "n", yaxt="n",
#      cex.lab = pt10, cex.axis = pt8)
#
# segments(x0=survivalPosteriorSummary$lo.surv,x1=survivalPosteriorSummary$hi.surv,
#          y0=g1PosteriorSummary$mode.g1, lwd=1)
# segments(x0=survivalPosteriorSummary$mode.surv,
#          y0=g1PosteriorSummary$lo.g1, y1=g1PosteriorSummary$hi.g1,
#          lwd=1)
# points(survivalPosteriorSummary$mode.surv,g1PosteriorSummary$mode.g1,
#        pch=21,col='black',bg='white',cex=1.25)
#
# axis(1, seq(0,1,by=.1), padj = -2,
#      labels = seq(0,1,by=.1), line = 0,
#      col = NA, col.ticks = 1, cex.axis = pt8)
# axis(2, seq(0,1,by=.1),
#      labels = seq(0,1,by=.1), las = 1, line = 0, hadj= .5,
#      col = NA, col.ticks = 1, cex.axis = pt8)
# mtext("Germination probability",
#       side=2,line=1.5,adj=.5,col='black',cex=pt10)
# mtext("Seed survival probability",
#       side=1,line=1,adj=.5,col='black',cex=pt10)
#
# text(x=.3,y=.5,
#      paste0("Pearson's r=",round(CI.correlation[2],2)),
#      cex=pt10)
# #abline(a=0,b=1)
#
# par(fig=c(0,10,0,4)/10)
# par(new=T)
# # plot posterior of correlation coefficient
# hist(posterior.correlation,breaks = 50,
#      main = "", xlab = "", ylab='', xaxt='n',yaxt='n',
#      xlim = c(-1, 1), ylim=c(0,2),
#      freq = FALSE, col = "gray75", border='white',
#      cex.lab = 1.25,cex.axis=1.5)
#
# # as in Duskey dissertation
# segments(x0=HPDI.correlation[1],x1=HPDI.correlation[2],y0=2,lwd=1.25)
# #segments(x0=CI.correlation[2],y0=2.3,y1=2.5,lwd=1.5)
# segments(x0=HPDI.correlation.lo[1],x1=HPDI.correlation.lo[2],y0=2,lwd=2.5)
# points(x=median(posterior.correlation),y=2,pch=21,col='black',bg='white',cex=1.25)
# #segments(x0=CI.correlation[2],y0=2.4,y1=0,lwd=2,col='gray90')
#
# ylab=formatC(seq(0,2,by=.5), format="f", digits=1)
# axis(1, seq(-2,2,by=.2), padj = -2,
#      labels = seq(-2,2,by=.2), las = 1, line = 0,
#      col = NA, col.ticks = 1, cex.axis = pt8)
# axis(2, seq(0,2,by=.5), hadj = .5,
#      labels = ylab, las = 2, line = -.15,
#      col = NA, col.ticks = 1, cex.axis = pt8)
# segments(x0=-1,x1=1,y0=-.08,lwd=2)
# segments(x0=-1.05,y0=0,y1=2,lwd=1.5)
# mtext("Density",
#       side=2,line=1.5,adj=.5,col='black',cex=pt10)
# mtext("Correlation of germination and survival",
#       side=1,line=2,adj=.5,col='black',cex=pt10)
#
# dev.off()





# ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/germ_surv_correlation.pdf",
#        plot=g1,width=4,height=4)


# c0<-sim(rho=-.75,nsites=20)
#
# df.sim <- data.frame(df$site,c0)
# g1.hypothesis<-ggplot(df.sim,aes(x=X,y=V2)) +
#   geom_point(color='gray') +
#   #geom_text_repel(size=3,color="black") +
#   theme_bw()+
#   xlab("Probability of seed survival [P(S)]") +
#   ylab("Mean germination probability [P(G)]") +
#   annotate("text", label =  paste0("Pearson's r=-.75"), x = .6, y = .29, size = 4,color='gray') +
#   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
#   scale_x_continuous(limits = c(0,.75), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#
# # ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/germ_surv_correlation-hypothesis.pdf",
# #        plot=g1.hypothesis,width=4,height=4)
#
# g1.blank<-ggplot(df.sim,aes(x=X,y=V2)) +
#   #geom_point(color='gray') +
#   #geom_text_repel(size=3,color="black") +
#   theme_bw()+
#   xlab("Probability of seed survival [P(S)]") +
#   ylab("Mean germination probability [P(G)]") +
#   annotate("text", label =  paste0("Pearson's r=-.75"), x = .6, y = .29, size = 4,color='gray') +
#   theme_bw() + #xlim(c(0,8)) + ylim(c(0,.3)) +
#   scale_x_continuous(limits = c(0,.75), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(0,.31), expand = c(0, 0)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggsave(filename="~/Dropbox/clarkiaSeedBanks/products/figures/germ_surv_correlation-blank.pdf",
#        plot=g1.blank,width=4,height=4)
