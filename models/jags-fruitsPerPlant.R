
model {

  # PRIORS ------------------------------------------------

  for(i in 1:n_site){

    # TOTAL FRUIT EQUIVALENTS -----------------------------------------

    # population-level intercept for each population
    nu_tfe[i] ~ dgamma(1, 1)
  #  nu_tfe[i] ~ dnorm(0, 1)
    # place prior on standard deviation on random year effects
    # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
    sigma0_tfe[i] ~ dnorm(0, 1) T(0,)
    # calculate precision, which is what JAGS uses
    tau0_tfe[i] <- 1/(sigma0_tfe[i]*sigma0_tfe[i])

    sigma_overdisp_tfe[i] ~ dnorm(0, 1) T(0,)
    tau_overdisp_tfe[i]<-1/( sigma_overdisp_tfe[i]*sigma_overdisp_tfe[i])

    # TOTAL/DAMAGED FRUITS ------------------------------------------------

    # population-level intercept for each population
    nu_tot[i] ~ dgamma(1, 1)
    # place prior on standard deviation on random year effects
    # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
    sigma0_tot[i] ~ dnorm(0, 1) T(0,)
    # calculate precision, which is what JAGS uses
    tau0_tot[i] <- 1/(sigma0_tot[i]*sigma0_tot[i])

    sigma_overdisp_tot[i] ~ dnorm(0, 1) T(0,)
    tau_overdisp_tot[i]<-1/( sigma_overdisp_tot[i]*sigma_overdisp_tot[i])

    # Proportion of fruits damaged  -----------------------------------------
    # population-level intercept for each population
    mu0[i] ~  dnorm(0, 1)
    # place prior on standard deviation on random year effects
    sigma0[i] ~ dnorm(0, 1) T(0,)
    # calculate precision, which is what JAGS uses
    tau0[i] <- 1/(sigma0[i]*sigma0[i])

  }

  # indexed by site*year combinations
  for(i in 1:n_siteYearIndex_tfe){

    # for each site*year combination, use the site (site_observed)
    # to get the population-level intercept and standard deviation of random year effects
    # draw the random year effects
    mu_tfe[i] ~ dlnorm(nu_tfe[site_tfe_observed[i]], tau0_tfe[site_tfe_observed[i]])
    mu_log_tfe[i] <- log(mu_tfe[i])

  }

  # indexed by site*year combinations
  for(i in 1:n_siteYearIndex_tot){
    # for each site*year combination, use the site (site_observed)
    # to get the population-level intercept and standard deviation of random year effects
    # draw the random year effects
    mu_tot[i] ~ dlnorm(nu_tot[site_tot_observed[i]], tau0_tot[site_tot_observed[i]])
    mu_log_tot[i] <- log(mu_tot[i])

    mu[i] ~ dnorm(mu0[site_tot_observed[i]], tau0[site_tot_observed[i]])

  }


  # LIKELIHOODS -------------------------------------------------------------

  for (i in 1:n){

    eps_tfe_overdisp[i] ~ dnorm(0,tau_overdisp_tfe[site[i]])

    #  TOTAL FRUIT EQUIVALENTS -------------------------------------------------
    log(z_tfe[i]) <- mu_log_tfe[siteYearIndex_tfe[i]] + eps_tfe_overdisp[i]
    y_tfe[i] ~ dpois(z_tfe[i])

    # Posterior predictive and chi-squared
    # y_tfe.sim[i] ~ dpois(z_tfe[i])
    # chi2.tfe.obs[i] <- pow((y_tfe[i] - z_tfe[i]),2) / (z_tfe[i])
    # chi2.tfe.sim[i] <- pow((y_tfe.sim[i]- z_tfe[i]),2) / (z_tfe[i])

    # rss.tfe.obs[i] <- pow((y_tfe[i] - z_tfe[i]),2)
    # rss.tfe.sim[i] <- pow((y_tfe.sim[i]- z_tfe[i]),2)

  }

  for (i in 1:n2){

    eps_tot_overdisp[i] ~ dnorm(0,tau_overdisp_tot[site2[i]])

    #  UNDAMAGED FRUITS -------------------------------------------------
    log(z_tot[i]) <- mu_log_tot[siteYearIndex_tot[i]] + eps_tot_overdisp[i]
    y_tot[i] ~ dpois(z_tot[i])

    # Posterior predictive and chi-squared
    # y_tot.sim[i] ~ dpois(z_tot[i])
    # chi2.tot.obs[i] <- pow((y_tot[i] - z_tot[i]),2) / (z_tot[i])
    # chi2.tot.sim[i] <- pow((y_tot.sim[i]- z_tot[i]),2) / (z_tot[i])

    # rss.tot.obs[i] <- pow((y_tot[i] - z_tot[i]),2)
    # rss.tot.sim[i] <- pow((y_tot.sim[i]- z_tot[i]),2)

    #  DAMAGED FRUITS -------------------------------------------------
    # proportion fruits damaged
    logit(prop_dam[i]) <- mu[siteYearIndex_tot[i]]

    # get number of seeds per damaged fruit
    y_dam[i] ~ dbinom(prop_dam[i],y_tot[i])

    # Posterior predictive and chi-squared
    # y_dam.sim[i] ~ dbinom(prop_dam[i],y_tot[i])
    # chi2.dam.obs[i] <- pow((y_dam[i] - prop_dam[i]*y_tot[i]),2) / (prop_dam[i]*y_tot[i]+0.001)
    # chi2.dam.sim[i] <- pow((y_dam.sim[i] - prop_dam[i]*y_tot[i]),2) / (prop_dam[i]*y_tot[i]+0.001)

    # rss.dam.obs[i] <- pow((y_dam[i] - prop_dam[i]*y_tot[i]),2)
    # rss.dam.sim[i] <- pow((y_dam.sim[i] - prop_dam[i]*y_tot[i]),2)
  }

}
