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
siteAbiotic <- read.csv("~/Dropbox/clarkia-demography-projects/data/siteAbioticData.csv",header=TRUE)

position<-siteAbiotic %>%
  dplyr::select(site,easting) %>%
  dplyr::mutate(easting=easting/1000)

siteIndex <- order(position$easting,decreasing=FALSE)
siteNames = unique(position$site)[siteIndex]

# ---
# - Read in germination and survival estimates ----
# ---

source("scripts/006_testHypotheses/00_utilityFunctions.R")

# ---
# - Colors ----
# ---

library(khroma)

################################################################################
# Germination
#################################################################################

g1 <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/g1-ex1-population-year-level-mat.RDS")
g1.pop <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/g1-population-level.RDS")

vr.popYear <- g1
vr.pop <- g1.pop

df.list <- list()
for(i in 1:20){
  tmp<-vr.popYear[[i]][,1:3]
  df.list[[i]]<-data.frame(site=unique(position$site)[i],
                           year = c(2006, 2007, 2008),
                           median=apply(tmp,2,median),
                           mode = apply(tmp,2,posterior.mode),
                           mean = apply(tmp,2,mean),
                           ci.lo = apply(tmp,2,function(x) hdi(x,credMass=0.95))[1,], 
                           ci.hi = apply(tmp,2,function(x) hdi(x,credMass=0.95))[2,],
                           ci.lo2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[1,], 
                           ci.hi2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[2,])
}
df.summary<-do.call(rbind,df.list)

vr.site=data.frame(site=unique(position$site), 
                   med= apply(vr.pop,2,posterior.mode),
                   ord= unique(position$easting),
                   ci.lo.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.025)), 
                   ci.hi.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.975)),
                   ci.lo2.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.25)), 
                   ci.hi2.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.75)))

vr.site=cbind(vr.site,t(hdi(vr.pop, credMass = 0.95) ),
              t(hdi(vr.pop, credMass = 0.5) ))
names(vr.site)[10:11] = c("lower2","upper2")

interannualG1 <- df.summary %>%
  dplyr::left_join(vr.site, by="site") %>%
  dplyr::arrange(ord) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(mode),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  # mutate(id = row_number()) %>% 
  ggplot(aes(x = year , y = mode)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper),
            alpha = 1/5,
            fill = "gray90") +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower2, ymax = upper2),
            alpha = 1/5,
            fill = "gray80") +
  geom_hline(aes(yintercept=med),linetype='dotted') +
  
  
  geom_linerange(aes(x=year,ymin=ci.lo,ymax=ci.hi),size=.25) +
  geom_linerange(aes(x=year,ymin=ci.lo2,ymax=ci.hi2),size=.5) +
  geom_point(aes(color=year)) +
  
  # coord_flip() +
  facet_grid(. ~ site, scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(1,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(
    axis.text.x=element_blank()) +
  xlab("Year") + ylab("Probability of germination (g1)") + labs(color="Year") +
  theme(legend.position="bottom",
        legend.background = element_rect(fill="gray95"),
        legend.key = element_rect(fill="gray95")) +
  scale_color_highcontrast()


ggsave(filename=paste0("outputs/005_calculatePopulationModelParameters/05_graphicalSummaries/interannual-G1.pdf"),
       plot=interannualG1,width=10,height=4)

################################################################################
# Seed survival s0
#################################################################################

s0 <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/s0-ex1-population-year-level-mat.RDS")
s0.pop <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s0-population-level.RDS")

vr.popYear <- s0
vr.pop <- s0.pop

df.list <- list()
for(i in 1:20){
  tmp<-vr.popYear[[i]][,2:3]
  df.list[[i]]<-data.frame(site=unique(position$site)[i],
                           year = c(2006, 2007),
                           median=apply(tmp,2,median),
                           mode = apply(tmp,2,posterior.mode),
                           mean = apply(tmp,2,mean),
                           ci.lo = apply(tmp,2,function(x) hdi(x,credMass=0.95))[1,], 
                           ci.hi = apply(tmp,2,function(x) hdi(x,credMass=0.95))[2,],
                           ci.lo2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[1,], 
                           ci.hi2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[2,])
}
df.summary<-do.call(rbind,df.list)

vr.site=data.frame(site=unique(position$site), 
                   med= apply(vr.pop,2,posterior.mode),
                   ord= unique(position$easting),
                   ci.lo.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.025)), 
                   ci.hi.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.975)),
                   ci.lo2.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.25)), 
                   ci.hi2.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.75)))

vr.site=cbind(vr.site,t(hdi(vr.pop, credMass = 0.95) ),
              t(hdi(vr.pop, credMass = 0.5) ))
names(vr.site)[10:11] = c("lower2","upper2")

interannualS0 <- df.summary %>%
  dplyr::left_join(vr.site, by="site") %>%
  dplyr::arrange(ord) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(mode),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  # mutate(id = row_number()) %>% 
  ggplot(aes(x = year , y = mode)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper),
            alpha = 1/5,
            fill = "gray90") +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower2, ymax = upper2),
            alpha = 1/5,
            fill = "gray80") +
  geom_hline(aes(yintercept=med),linetype='dotted') +
  
  geom_linerange(aes(x=year,ymin=ci.lo,ymax=ci.hi),size=.25) +
  geom_linerange(aes(x=year,ymin=ci.lo2,ymax=ci.hi2),size=.5) +
  geom_point(aes(color=year)) +
  
  # coord_flip() +
  facet_grid(. ~ site, scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(1,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(
    axis.text.x=element_blank()) +
  xlab("Year") + ylab("Probability of seed survival (s0)") + labs(color="Year") +
  theme(legend.position="bottom",
        legend.background = element_rect(fill="gray95"),
        legend.key = element_rect(fill="gray95"))+
  scale_color_highcontrast()

ggsave(filename=paste0("outputs/005_calculatePopulationModelParameters/05_graphicalSummaries/interannual-S0.pdf"),
       plot=interannualS0,width=10,height=4)


################################################################################
# Seed survival s1
#################################################################################

s1 <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/s1-ex1-population-year-level-mat.RDS")
s1.pop <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s1-population-level.RDS")

vr.popYear <- s1
vr.pop <- s1.pop

df.list <- list()
for(i in 1:20){
  tmp<-vr.popYear[[i]][,1:3]
  df.list[[i]]<-data.frame(site=unique(position$site)[i],
                           year = c(2006, 2007, 2008),
                           median=apply(tmp,2,median),
                           mode = apply(tmp,2,posterior.mode),
                           mean = apply(tmp,2,mean),
                           ci.lo = apply(tmp,2,function(x) hdi(x,credMass=0.95))[1,], 
                           ci.hi = apply(tmp,2,function(x) hdi(x,credMass=0.95))[2,],
                           ci.lo2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[1,], 
                           ci.hi2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[2,])
}
df.summary<-do.call(rbind,df.list)

vr.site=data.frame(site=unique(position$site), 
                   med= apply(vr.pop,2,posterior.mode),
                   ord= unique(position$easting),
                   ci.lo.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.025)), 
                   ci.hi.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.975)),
                   ci.lo2.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.25)), 
                   ci.hi2.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.75)))

vr.site=cbind(vr.site,t(hdi(vr.pop, credMass = 0.95) ),
              t(hdi(vr.pop, credMass = 0.5) ))
names(vr.site)[10:11] = c("lower2","upper2")

interannualS1 <- df.summary %>%
  dplyr::left_join(vr.site, by="site") %>%
  dplyr::arrange(ord) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(mode),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  # mutate(id = row_number()) %>% 
  ggplot(aes(x = year , y = mode)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper),
            alpha = 1/5,
            fill = "gray90") +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower2, ymax = upper2),
            alpha = 1/5,
            fill = "gray80") +
  geom_hline(aes(yintercept=med),linetype='dotted') +
  
  geom_linerange(aes(x=year,ymin=ci.lo,ymax=ci.hi),size=.25) +
  geom_linerange(aes(x=year,ymin=ci.lo2,ymax=ci.hi2),size=.5) +
  geom_point(aes(color=year)) +
  
  # coord_flip() +
  facet_grid(. ~ site, scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(1,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(
    axis.text.x=element_blank()) +
  xlab("Year") + ylab("Probability of seed survival (s1)") + labs(color="Year") +
  theme(legend.position="bottom",
        legend.background = element_rect(fill="gray95"),
        legend.key = element_rect(fill="gray95"))+
  scale_color_highcontrast()

ggsave(filename=paste0("outputs/005_calculatePopulationModelParameters/05_graphicalSummaries/interannual-S1.pdf"),
       plot=interannualS1,width=10,height=4)

################################################################################
# Seed survival s2
#################################################################################

s2 <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/s2-ex1-population-year-level-mat.RDS")
s2.pop <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s2-population-level.RDS")

vr.popYear <- s2
vr.pop <- s2.pop

df.list <- list()
for(i in 1:20){
  tmp<-vr.popYear[[i]][,1:3]
  df.list[[i]]<-data.frame(site=unique(position$site)[i],
                           year = c(2006, 2007, 2008),
                           median=apply(tmp,2,median),
                           mode = apply(tmp,2,posterior.mode),
                           mean = apply(tmp,2,mean),
                           ci.lo = apply(tmp,2,function(x) hdi(x,credMass=0.95))[1,], 
                           ci.hi = apply(tmp,2,function(x) hdi(x,credMass=0.95))[2,],
                           ci.lo2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[1,], 
                           ci.hi2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[2,])
}
df.summary<-do.call(rbind,df.list)

vr.site=data.frame(site=unique(position$site), 
                   med= apply(vr.pop,2,posterior.mode),
                   ord= unique(position$easting),
                   ci.lo.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.025)), 
                   ci.hi.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.975)),
                   ci.lo2.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.25)), 
                   ci.hi2.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.75)))

vr.site=cbind(vr.site,t(hdi(vr.pop, credMass = 0.95) ),
              t(hdi(vr.pop, credMass = 0.5) ))
names(vr.site)[10:11] = c("lower2","upper2")

interannualS2 <- df.summary %>%
  dplyr::left_join(vr.site, by="site") %>%
  dplyr::arrange(ord) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(mode),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  # mutate(id = row_number()) %>% 
  ggplot(aes(x = year , y = mode)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper),
            alpha = 1/5,
            fill = "gray90") +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower2, ymax = upper2),
            alpha = 1/5,
            fill = "gray80") +
  geom_hline(aes(yintercept=med),linetype='dotted') +
  
  geom_linerange(aes(x=year,ymin=ci.lo,ymax=ci.hi),size=.25) +
  geom_linerange(aes(x=year,ymin=ci.lo2,ymax=ci.hi2),size=.5) +
  geom_point(aes(color=year)) +
  
  # coord_flip() +
  facet_grid(. ~ site, scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(1,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(
    axis.text.x=element_blank()) +
  xlab("Year") + ylab("Probability of seed survival (s2)") + labs(color="Year") +
  theme(legend.position="bottom",
        legend.background = element_rect(fill="gray95"),
        legend.key = element_rect(fill="gray95"))+
  scale_color_highcontrast()

ggsave(filename=paste0("outputs/005_calculatePopulationModelParameters/05_graphicalSummaries/interannual-S2.pdf"),
       plot=interannualS2,width=10,height=4)

################################################################################
# Seed survival s3
#################################################################################

s3 <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/s3-ex1-population-year-level-mat.RDS")
s3.pop <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/s3-population-level.RDS")

vr.popYear <- s3
vr.pop <- s3.pop

df.list <- list()
for(i in 1:20){
  tmp<-vr.popYear[[i]][,2:3]
  df.list[[i]]<-data.frame(site=unique(position$site)[i],
                           year = c(2007, 2008),
                           median=apply(tmp,2,median),
                           mode = apply(tmp,2,posterior.mode),
                           mean = apply(tmp,2,mean),
                           ci.lo = apply(tmp,2,function(x) hdi(x,credMass=0.95))[1,], 
                           ci.hi = apply(tmp,2,function(x) hdi(x,credMass=0.95))[2,],
                           ci.lo2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[1,], 
                           ci.hi2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[2,])
}
df.summary<-do.call(rbind,df.list)

vr.site=data.frame(site=unique(position$site), 
                   med= apply(vr.pop,2,posterior.mode),
                   ord= unique(position$easting),
                   ci.lo.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.025)), 
                   ci.hi.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.975)),
                   ci.lo2.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.25)), 
                   ci.hi2.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.75)))

vr.site=cbind(vr.site,t(hdi(vr.pop, credMass = 0.95) ),
              t(hdi(vr.pop, credMass = 0.5) ))
names(vr.site)[10:11] = c("lower2","upper2")

interannualS3 <- df.summary %>%
  dplyr::left_join(vr.site, by="site") %>%
  dplyr::arrange(ord) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(desc(mode),.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  # mutate(id = row_number()) %>% 
  ggplot(aes(x = year , y = mode)) + 
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper),
            alpha = 1/5,
            fill = "gray90") +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower2, ymax = upper2),
            alpha = 1/5,
            fill = "gray80") +
  geom_hline(aes(yintercept=med),linetype='dotted') +
  
  geom_linerange(aes(x=year,ymin=ci.lo,ymax=ci.hi),size=.25) +
  geom_linerange(aes(x=year,ymin=ci.lo2,ymax=ci.hi2),size=.5) +
  geom_point(aes(color=year)) +
  
  # coord_flip() +
  facet_grid(. ~ site, scales="free_x", space="free_x") +
  theme_bw() +
  theme(panel.spacing=unit(1,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(
    axis.text.x=element_blank()) +
  xlab("Year") + ylab("Probability of seed survival (s3)") + labs(color="Year") +
  theme(legend.position="bottom",
        legend.background = element_rect(fill="gray95"),
        legend.key = element_rect(fill="gray95"))+
  scale_color_manual(values=c("#DDAA33", "#BB5566"))

ggsave(filename=paste0("outputs/005_calculatePopulationModelParameters/05_graphicalSummaries/interannual-S3.pdf"),
       plot=interannualS3,width=10,height=4)

################################################################################
# sigma
#################################################################################

sigma <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/sigma-population-year-level-mat.RDS")
sigma.pop <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/sigma-population-level.RDS")

vr.popYear <- sigma
vr.pop <- sigma.pop

df.list <- list()
for(i in 1:20){
  tmp<-vr.popYear[[i]][,1:15]
  df.list[[i]]<-data.frame(site=unique(position$site)[i],
                           year = c(2006:2020),
                           median=apply(tmp,2,median),
                           mode = apply(tmp,2,posterior.mode),
                           mean = apply(tmp,2,mean),
                           ci.lo = apply(tmp,2,function(x) hdi(x,credMass=0.95))[1,], 
                           ci.hi = apply(tmp,2,function(x) hdi(x,credMass=0.95))[2,],
                           ci.lo2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[1,], 
                           ci.hi2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[2,])
}
df.summary<-do.call(rbind,df.list)

vr.site=data.frame(site=unique(position$site), 
                   med= apply(vr.pop,2,posterior.mode),
                   ord= unique(position$easting),
                   ci.lo.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.025)), 
                   ci.hi.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.975)),
                   ci.lo2.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.25)), 
                   ci.hi2.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.75)))

vr.site=cbind(vr.site,t(hdi(vr.pop, credMass = 0.95) ),
              t(hdi(vr.pop, credMass = 0.5) ))
names(vr.site)[10:11] = c("lower2","upper2")

interannualSigma <- df.summary %>%
  dplyr::left_join(vr.site, by="site") %>%
  dplyr::arrange(ord) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(mode,.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(x = year , y = mode)) + 
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper),
  #           alpha = 1/5,
  #           fill = "gray90") +
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower2, ymax = upper2),
  #           alpha = 1/5,
  #           fill = "gray80") +
  # geom_hline(aes(yintercept=med),linetype='dotted') +
  
  geom_linerange(aes(x=year,ymin=ci.lo,ymax=ci.hi),size=.25) +
  geom_linerange(aes(x=year,ymin=ci.lo2,ymax=ci.hi2),size=.5) +
  geom_point(aes(color=year)) +
  
  # coord_flip() +
  facet_wrap(. ~ site, nrow=2) +
  theme_bw() +
  theme(panel.spacing=unit(1,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(
    axis.text.x=element_blank()) +
  xlab("Year") + ylab("Probability of seedling survival to fruiting") + labs(color="Year") +
  theme(legend.position="bottom",
        legend.background = element_rect(fill="gray95"),
        legend.key = element_rect(fill="gray95")) +
  scale_color_batlow(discrete=TRUE,reverse=TRUE)

ggsave(filename=paste0("outputs/005_calculatePopulationModelParameters/05_graphicalSummaries/interannual-Sigma.pdf"),
       plot=interannualSigma,width=10,height=6)
# 
# annualSigma<-df.summary %>%
#   dplyr::left_join(vr.site, by="site") %>%
#   dplyr::arrange(med) %>% 
#   dplyr::mutate(site=factor(site,levels=unique(site))) %>%
#   dplyr::group_by(site) %>%
#   dplyr::arrange(mode,.by_group = TRUE) %>%
#   dplyr::mutate(year=factor(year)) %>%
#   mutate(id = row_number()) %>% 
#   ggplot(aes(x = year , y = mode)) + 
#   # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper),
#   #           alpha = 1/5,
#   #           fill = "gray90") +
#   # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower2, ymax = upper2),
#   #           alpha = 1/5,
#   #           fill = "gray80") +
#   #geom_hline(aes(yintercept=med),linetype='dotted') +
#   
#   geom_linerange(aes(x=year,ymin=ci.lo,ymax=ci.hi),size=.25) +
#   geom_linerange(aes(x=year,ymin=ci.lo2,ymax=ci.hi2),size=.5) +
#   geom_point() +
#   geom_line(aes(x=year,y=mode,group=site)) +
#   
#   #coord_flip() +
#   facet_wrap( ~ site) +
#   theme_bw() +
#   theme(panel.spacing=unit(2,"pt"), 
#         panel.border=element_rect(colour="grey50", fill=NA)) +
#   theme(
#     axis.text.x=element_blank()) +
#   xlab("Year") + ylab("Probability of seedling survival to fruiting") + labs(color="Year") 
# 
# ggsave(filename=paste0("outputs/005_calculatePopulationModelParameters/05_graphicalSummaries/annual-Sigma.pdf"),
#        plot=annualSigma,width=7,height=7)


################################################################################
# fecundity
#################################################################################

fec <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/combinedF-population-year-level-mat.RDS")

vr.popYear <- fec

df.list <- list()
for(i in 1:20){
  tmp<-vr.popYear[[i]][,1:15]
  df.list[[i]]<-data.frame(site=unique(position$site)[i],
                           year = c(2006:2020),
                           median=apply(tmp,2,median),
                           mode = apply(tmp,2,posterior.mode),
                           mean = apply(tmp,2,mean),
                           ci.lo = apply(tmp,2,function(x) hdi(x,credMass=0.95))[1,], 
                           ci.hi = apply(tmp,2,function(x) hdi(x,credMass=0.95))[2,],
                           ci.lo2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[1,], 
                           ci.hi2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[2,])
}
df.summary<-do.call(rbind,df.list)

vr.site=data.frame(site=unique(position$site), 
                   med= apply(vr.pop,2,posterior.mode),
                   ord= unique(position$easting))

# vr.site=cbind(vr.site,t(hdi(vr.pop, credMass = 0.95) ),
#               t(hdi(vr.pop, credMass = 0.5) ))
# names(vr.site)[10:11] = c("lower2","upper2")

interannualFec<- df.summary %>%
  dplyr::left_join(vr.site, by="site") %>%
  dplyr::arrange(ord) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(mode,.by_group = TRUE) %>%
  dplyr::mutate(year.numeric=(year)) %>%
  dplyr::mutate(year=as.factor(year)) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(x = year.numeric , y = mode)) + 
  geom_rect(aes(xmin = 2005, xmax = 2012, ymin = -Inf, ymax = Inf),
            alpha = 1/20,
            fill = "gray95") +
  geom_linerange(aes(x=year.numeric,ymin=ci.lo,ymax=ci.hi),size=.25) +
  geom_linerange(aes(x=year.numeric,ymin=ci.lo2,ymax=ci.hi2),size=.5) +
  geom_point(aes(color=year)) +
  
  # coord_flip() +
  facet_wrap(. ~ site, nrow=2) +
  theme_bw() +
  theme(panel.spacing=unit(1,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(
    axis.text.x=element_blank()) +
  xlab("Year") + ylab("Fruits per plant") + labs(color="Year") +
  theme(legend.position="bottom",
        legend.background = element_rect(fill="gray95"),
        legend.key = element_rect(fill="gray95")) +
  scale_color_batlow(discrete=TRUE,reverse=TRUE)

ggsave(filename=paste0("outputs/005_calculatePopulationModelParameters/05_graphicalSummaries/interannual-fec.pdf"),
       plot=interannualFec,width=10,height=6)

################################################################################
# phi
#################################################################################

phi <- readRDS("outputs/005_calculatePopulationModelParameters/03_populationModelParametersMatrix/phi-population-year-level-mat.RDS")
phi.pop <- readRDS("outputs/005_calculatePopulationModelParameters/02_populationModelParameters/phi-population-level.RDS")

vr.popYear <- phi
vr.pop <- phi.pop

df.list <- list()
for(i in 1:20){
  tmp<-vr.popYear[[i]][,1:15]
  df.list[[i]]<-data.frame(site=unique(position$site)[i],
                           year = c(2006:2020),
                           median=apply(tmp,2,median),
                           mode = apply(tmp,2,posterior.mode),
                           mean = apply(tmp,2,mean),
                           ci.lo = apply(tmp,2,function(x) hdi(x,credMass=0.95))[1,], 
                           ci.hi = apply(tmp,2,function(x) hdi(x,credMass=0.95))[2,],
                           ci.lo2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[1,], 
                           ci.hi2 = apply(tmp,2,function(x) hdi(x,credMass=0.5))[2,])
}
df.summary<-do.call(rbind,df.list)

vr.site=data.frame(site=unique(position$site), 
                   med= apply(vr.pop,2,posterior.mode),
                   ord= unique(position$easting),
                   ci.lo.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.025)), 
                   ci.hi.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.975)),
                   ci.lo2.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.25)), 
                   ci.hi2.pop = apply(vr.pop,2,function(x) quantile(x,probs=0.75)))

vr.site=cbind(vr.site,t(hdi(vr.pop, credMass = 0.95) ),
              t(hdi(vr.pop, credMass = 0.5) ))
names(vr.site)[10:11] = c("lower2","upper2")

interannualPhi<-df.summary %>%
  dplyr::left_join(vr.site, by="site") %>%
  dplyr::arrange(ord) %>% 
  dplyr::mutate(site=factor(site,levels=unique(site))) %>%
  dplyr::group_by(site) %>%
  dplyr::arrange(mode,.by_group = TRUE) %>%
  dplyr::mutate(year=factor(year)) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(x = year , y = mode)) + 
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper),
  #           alpha = 1/5,
  #           fill = "gray90") +
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = lower2, ymax = upper2),
  #           alpha = 1/5,
  #           fill = "gray80") +
  # geom_hline(aes(yintercept=med),linetype='dotted') +
   
  geom_linerange(aes(x=year,ymin=ci.lo,ymax=ci.hi),size=.25) +
  geom_linerange(aes(x=year,ymin=ci.lo2,ymax=ci.hi2),size=.5) +
  geom_point(aes(color=year)) +
  
  # coord_flip() +
  facet_wrap(. ~ site, nrow=2) +
  theme_bw() +
  theme(panel.spacing=unit(1,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA)) +
  theme(
    axis.text.x=element_blank()) +
  xlab("Year") + ylab("Seeds per undamaged fruit") + labs(color="Year") +
  theme(legend.position="bottom",
        legend.background = element_rect(fill="gray95"),
        legend.key = element_rect(fill="gray95")) +
  scale_color_batlow(discrete=TRUE,reverse=TRUE)

ggsave(filename=paste0("outputs/005_calculatePopulationModelParameters/05_graphicalSummaries/interannual-phi.pdf"),
       plot=interannualPhi,width=10,height=6)
