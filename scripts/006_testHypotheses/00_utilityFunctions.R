# --
#  Functions used to analyze results
# --

#' Apply a function to columns
#' 
#' @param x matrix object, generally a matrix of MCMC samples from a posterior 
#' @param toDo function to call 
#' 
#' @return numeric output of toDo, with dimensions according to call to toDo

cols_fun <- function(x,fun=toDo){
 return( apply(x,2,toDo) )
}

#' Sample-based estimate of the mode of an continuous density estimate
#' NOTE: this is based on rethinking::chainmode in the rethinking package
#' with modifications to account for missing values and lower bounds. 
#' We need this here because not all years have observations.
#' 
#' @param x samples from a posterior, using MCMC 
#' 
#' @return numeric estimate for the mode of the posterior
posterior.mode = function(x){
  if(!is.na(x[1])){ x.max=max(x)
  x.min=min(x)
  dres <- density( x ,from = x.min, to = x.max)
  modeParam <- dres$x[which.max(dres$y)]}else if(is.na(x[1])){
    modeParam <- NA
  }
  return(modeParam)
}


#' COMMENTED OUT BECAUSE THIS WAS NOT SUED
#' Sample-based estimate for the geometric mean
#' NOTE: from https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
#' NOTE: from https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
#' 
#' @param x values for which to calculate the geometric mean
#' @param na.rm logical for whether to remove NAs from the calculation
#' 
#' @return numeric estimate for the geometric mean
# gm_mean = function(x, na.rm=TRUE){
#   exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
# }

#' Sample-based estimate for the geometric mean
#' 
#' @param x values for which to calculate the geometric mean
#' @param na.rm logical for whether to remove NAs from the calculation
#' 
#' @return numeric estimate for the geometric mean
gm_mean = function(x, na.rm=TRUE){
  gm_est = exp(sum(log(x), na.rm=na.rm) / length(x))
  return(gm_est)
}

#' Sample-based estimate for the geometric standard deviation.
#' The implementation accounts for 0s by adding 0.5 to all values.
#' The geometric mean is then calculated by exponentiating the 
#' arithmetic mean of the natural log of the values. 
#' The geometric standard deviation is then calculated using 
#' this geometric mean.
#' 
#' @param x values for which to calculate the geometric standard deviation
#' 
#' @return numeric value for the geometric standard deviation
gsd.am <- function(x){
  x=x+.5
  n = length(x[!is.na(x)])
  geom_mu = exp(mean(log(x),na.rm=TRUE))
  geom_sd <- exp(sqrt(sum((log(x/geom_mu))^2,na.rm=TRUE)/(n-1)))
  return(geom_sd)
}

#' Calculates the geometric standard deviation.
#' This does not account for 0s.
#' 
#' @param x values for which to calculate the geometric standard deviation
#' 
#' @return numeric value for the geometric standard deviation
gsd <- function(x){
  y <- exp(sd(log(x)))
  return(y)
}

#' Estimates percentile intervals for samples from a posterior.
#' The percentile intervals assign equal mass to both tails of the distribution.
#' Follows the methods for PI/PCI in the rethinking package
#' 
#' @param zc MCMC object
#' @param x names of parameters to be summarized
#' @param prob by default, calculates the 50 and 95% percentile intervals and the median
#' 
#' @return numeric value for the geometric standard deviation
posteriorPercentileInterval<-function(postSamples = zc, paramNames="param", prob = c(.025, .25, .5, .75, .975)){
  chain<-MCMCchains(postSamples, params = paramNames)
  percentileInterval <- t(apply(p, 2, FUN = function(x) quantile(x, prob)))
  return(percentileInterval)
}
