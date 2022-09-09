####
####
# Script to create figures that illustrate
# the structure of the models for aboveground observations
# using seedling survival to fruiting as the example dataset
# the figures link observations -> estimated parameters
# Produces graphs in Figure 2A focused on aboveground parts of demography
####
####

# - Load packages ----
library(rjags)
library(tidybayes)
library(tidyverse)
library(MCMCvis)

# - Simulate dataset ----
set.seed(17)

# 30 observations
n.obs = 30
# per 5 years
n.years = 5

# Population mean and standard deviation
mu0 = -1
sigma0 = 1

# Simulate the annual mean for 5 years of observations
mu = rnorm(n.years,mean=mu0,sd=sigma0)
# for each of those 5 years, simulate 5 values for the standard deviation
# sigma = runif(n.years,min=0,max=1)

# pre allocate lists to hold output
# theta is each individual's intercept
# theta_p is each individual's probability of success
theta = theta_p = list()

# n is the number of seedlings at the start of the observation period, drawn from 1:1000
# y is the number of seedlings surviving to fruiting, sampled with probability theta_p
n = y = list()

# simulate individual observations
for(i in 1:n.years){
  theta[[i]] = mu[i] # rnorm(n.obs,mean=mu[i],sd=sigma[i])
  theta_p[[i]] = boot::inv.logit(theta[[i]])
  n[[i]] = sample(50,n.obs,replace=TRUE)
  y[[i]] = rbinom(n.obs, n[[i]], theta_p[[i]])
}

# turn observations into vectors
y = unlist(y)
n = unlist(n)
# assign "site" and "year" to the observations
site = rep(as.factor(1),length(n))
year = rep(as.factor(1:n.years),each=n.obs)

# construct a data frame 
df=data.frame(fruitplNumber=y,seedlingNumber=n,site=site,year=year,siteYearIndex=year,site_observed=site)

# - Prepare data for analysis with JAGS ----
# compose_data turns the data frame into a list for JAGS
data <- tidybayes::compose_data(df)

# detach tidyverse to avoid package conflicts
detach("package:tidyverse", unload=TRUE)

# - Set JAGS parameters and random seed ----
# number of iterations in the chain for adaptation
# number of iterations for burn-in
# number of samples in the final chain
# no thinning the chain here

n.adapt = 3000
n.update = 5000
n.iterations = 2000
n.thin = 1

# - Define functions that randomly set the initial conditions ----

initsMu0 <- function(samps = data$n_site){
  rnorm(n = samps, mean = 0, sd = 1)
}

initsSigma0 <- function(samps = data$n_site){
  extraDistr::rhnorm(n = samps, sigma = 1)
}

initsSigma <- function(rows = data$n_site, cols = data$n_year){
  matrix(extraDistr::rhnorm(n = rows*cols, sigma = 1), rows, cols)
}

# - Set the initial conditions for JAGS ----

inits <- list()

for(i in 1:3){
  inits[[i]] <- list(initsMu0(), initsSigma0(), initsSigma() )
  
  names(inits[[i]]) = c("mu0","sigma0","sigma")
  
}

# - Call JAGS and fit model ----

# tuning step (n.adapt)
jm = jags.model("models/jags-seedlingSurvival.R",
                data = data, inits = inits,n.chains = length(inits), n.adapt = n.adapt)

# burn-in step (n.update)
update(jm, n.iter = n.update)

# monitor these parameters in draws from posterior
parsToMonitor = c("mu0","sigma0","mu","sigma","theta")

# sample chain (n.iter)
samples.rjags = coda.samples(jm, 
                             variable.names = c(parsToMonitor), 
                             n.iter = n.iterations, thin = n.thin)

# - Summarize the model fit ----

# view summary of parameters
MCMCsummary(samples.rjags,params=c("mu0","sigma0","mu","sigma"))

# extract the posterior chains for all parameters
mu0.post = MCMCchains(samples.rjags,params="mu0")
sigma0.post = MCMCchains(samples.rjags,params="sigma0")
mu.post = MCMCchains(samples.rjags,params="mu")
sigma.post = MCMCchains(samples.rjags,params="sigma")

# - Function to create object for plotting ----
# extracts a kernel density with a smoothed resolution

f=function(x){
  tmp=density(x,from=min(x),to=max(x),adjust=2.5)
  df=cbind(tmp$x,tmp$y)
  df=rbind(c(tmp$x[1],0),df)
  df=rbind(df,c(tmp$x[length(tmp$x)],0))
  return(df)
}



# - Plot figures ----

# set font sizes
pt12 = 1
pt10 = 10/12
pt9 = 9/12
pt8 = 8/12
pt7 = 7/12
pt6 = 6/12

# - + Plot marginal posterior ----

tiff("products/figures/marginal-posterior.tif",
     units='px',height = 1600, width = 1600,res=800,compression='lzw')

par(mfrow=c(1,1),mar=c(.5,.5,.75,0),oma=c(.75,.75,0,0)+.1,mgp=c(3,.1,0))

plot(x=NA,y=NA,
     type='n',
     ylim=c(0,12),xlim=c(0,1),
     frame=FALSE,xaxt='n',yaxt='n')

year=1:5
sigma = boot::inv.logit(mu.post)

# for each year, plot the marginal posterior and then add the observations below
for(i in 1:5){
  df.tmp=sigma[,i]
  full.post=mu.post[,i] #rnorm(length(mu0.post[,i]),mean=mu0.post[,i],sd=sigma0.post[,i])
  full.post = boot::inv.logit(full.post)
  dat.tmp = df[df$year==i,] %>%
    dplyr::mutate(p = fruitplNumber/seedlingNumber) %>%
    dplyr::filter(!is.na(p))
  
  upper.limit=max(f(full.post)[,2])*.8
  polygon(y=.5+2*year[i]+f(full.post)[,2]/upper.limit,x=f(full.post)[,1],col='gray80',border='gray80')
  
  n.obs = length(dat.tmp$seedlingNumber)
  size = dat.tmp$seedlingNumber/max(df$seedlingNumber,na.rm=TRUE)
  points(y=.5+2*rep(i-0.1,n.obs)+rnorm(n.obs,0,.025),x=dat.tmp$p,
         pch = 19, cex = .45,col=rgb(0,0,0,.5))
}

# add axes
axis(1, at = seq(0,1,by=.2), labels = seq(0,1,by=.2), 
     las = 1,  tck=-0.01, padj=-1.1,
     col = 'black', col.ticks = 1, 
     cex.axis = pt6)
axis(2, at = 2*(c(2:6)-.5), labels = c(1:5), 
     las = 1,  tck=-0.01,
     col = 'black', col.ticks = 1, 
     cex.axis = pt6)

# add axis labels
mtext("Year",side=2,adj=.6,line=.5,cex=pt7,cex.lab=pt7)
mtext("Pr(seedling survival to fruiting)",
      side=1,adj=.4,padj=-.5,
      line=.5,cex=pt7,cex.lab=pt7)

# add the population-level marginal posterior distribution
df.tmp=boot::inv.logit(mu0.post)
full.post=rnorm(length(mu0.post),mean=mu0.post,sd=sigma0.post)
full.post = boot::inv.logit(full.post)

upper.limit=max(f(full.post)[,2])*.8
polygon(y=f(full.post)[,2]/upper.limit,x=f(full.post)[,1],col='gray80',border='gray80')

# add some text to the plot to emphasize the population and year level vs. population level parts
text(.375,11.9,
     labels = "Population and year level" ,cex = pt6, pos = 4)
text(.585,.9,
     labels = "Population level" ,cex = pt6, pos = 4)

abline(h=1.5,lty='dotted')

box()
mtext("Marginalized probabilitities",adj=0,cex=pt9)

dev.off()


# - + Plot full posterior distribution of parameter estimates ----

tiff("products/figures/parameter.tif",
     units='px',height = 1600, width = 1600,res=800,compression='lzw')

par(mfrow=c(1,2),mar=c(.5,.25,.75,0),oma=c(.75,1,0,.1)+.1,mgp=c(3,.1,0))

plot(x=NA,y=NA,
     type='n',
     ylim=c(0,12),xlim=c(-4,2),
     frame=FALSE,xaxt='n',yaxt='n')

year=1:5
sigma = (mu.post)

# for each year, plot the posterior distribution of the annual estimate for the
# mean seedling survival to fruiting, on the untransformed scale
for(i in 1:5){
  
  df.tmp=sigma[,i]
  dat.tmp = df[df$year==i,] %>%
    dplyr::mutate(p = fruitplNumber/seedlingNumber) %>%
    dplyr::filter(!is.na(p))
  
  upper.limit=max(f(df.tmp)[,2])*1
  polygon(y=.5+2*year[i]+f(df.tmp)[,2]/upper.limit,x=f(df.tmp)[,1],col='#d95f02',border='#d95f02')
  
  n.obs = length(dat.tmp$seedlingNumber)
  size = dat.tmp$seedlingNumber/max(df$seedlingNumber,na.rm=TRUE)
  
}

# add axes
axis(1, at= seq(-6,2,by=2), labels = seq(-6,2,by=2), 
     las = 1, col = 'black', padj=-1.1,
     col.ticks = 1,  tck=-0.01, 
     cex.axis = pt6)
axis(2, at=  2*(c(2:6)-.5), labels = c(1:5),
     las = 1,  tck=-0.01, 
     col = 'black', col.ticks = 1, 
     cex.axis = pt6)

# add population level marginal posterior distribution for mean seedling survival to fruiting 
df.tmp=(mu0.post)

upper.limit=max(f(df.tmp)[,2])*1
polygon(y=f(df.tmp)[,2]/upper.limit,x=f(df.tmp)[,1],col='#7570b3',border='#7570b3')

abline(h=1.5,lty='dotted')
box()

# add axis labels
mtext("Year",side=2,adj=.6,line=.5,cex=pt7,cex.lab=pt7)
mtext("Model parameters",
      side=3,adj=0,cex=pt9)
mtext("Means",
      side=1,adj=.5,padj=-.5,
      line=.5,cex=pt7,cex.lab=pt7)

plot(x=NA,NA,
     type='n',
     ylim=c(0,12),xlim=c(0,2),
     #axes=FALSE,
     frame=FALSE,xaxt='n',yaxt='n')

# sigma = (sigma.post)

# add population and year marginal posterior distribution for SD of seedling survival to fruiting 
# for(i in 1:5){
#   df.tmp=sigma[,i]
#   dat.tmp = df[df$year==i,] %>%
#     dplyr::mutate(p = fruitplNumber/seedlingNumber) %>%
#     dplyr::filter(!is.na(p))
#   upper.limit=max(f(df.tmp)[,2])*1
#   polygon(y=.5+2*year[i]+f(df.tmp)[,2]/upper.limit,x=f(df.tmp)[,1],col='#d95f02',border='#d95f02')
#   
#   n.obs = length(dat.tmp$seedlingNumber)
#   size = dat.tmp$seedlingNumber/max(df$seedlingNumber,na.rm=TRUE)
# }

rect(xleft=-1,xright=3,ybottom=1.5,ytop=13,border=FALSE,col='gray90')

# add axis
axis(1, at= seq(0,4,by=1), labels = seq(0,4,by=1), 
     las = 1, col = 'black', padj=-1.1,
     col.ticks = 1,  tck=-0.01, 
     cex.axis = pt6)

# add population level marginal posterior distribution for SD of seedling survival to fruiting 

df.tmp=boot::inv.logit(sigma0.post)

upper.limit=max(f(df.tmp)[,2])*1
polygon(y=f(df.tmp)[,2]/upper.limit,x=f(df.tmp)[,1],col='#7570b3',border='#7570b3')

abline(h=1.5,lty='dotted')

# add axis
mtext("Standard deviations",
      side=1,adj=.4,padj=-.5,
      line=.5,cex=pt7,cex.lab=pt7)
box()

dev.off()