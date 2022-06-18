model {

  # priors
  for(k in 1:n_site){

    # assign a population-level intercept to each population
    mu0[k] ~  dnorm(0, 1)
    # place prior on standard deviation of random year effects
    sigma0[k] ~ dnorm(0, 1) T(0,)
    # calculate the precision, which is what JAGS uses
    tau0[k] <- 1/(sigma0[k]*sigma0[k])

  }

  # indexed by site*year combinations
  for(i in 1:n_siteYearIndex){

    # for each site*year combination, use the site (site_observed)
    # to get the population-level intercept and standard deviation of random year effects
    # draw the random year effects
    mu[i] ~ dnorm(mu0[site_observed[i]], tau0[site_observed[i]])

  }

  # Likelihood ------------------------------------------------------------------

  for(i in 1:n){

    logit(theta[i]) <- mu[siteYearIndex[i]]

    # LIKELIHOOD
    fruitplNumber[i] ~ dbinom(theta[i], seedlingNumber[i])

    # POSTERIOR PREDICTIVE
    fruitplNumber_sim[i] ~ dbinom(theta[i], seedlingNumber[i])

    chi2.obs[i] <- pow((fruitplNumber[i]- theta[i]*seedlingNumber[i]),2) / (theta[i]*seedlingNumber[i]+.001)         # obs.
    chi2.sim[i] <- pow((fruitplNumber_sim[i]- theta[i]*seedlingNumber[i]),2) / (theta[i]*seedlingNumber[i]+.001)         # sim.

  }

}
