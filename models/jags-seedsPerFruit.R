
model {

  # PRIORS ------------------------------------------------

  for(i in 1:n_site3){

    # SEEDS PER UNDAMAGED FRUITS  -----------------------------------------

    # population-level intercept for each population
    nu_seeds[i] ~ dgamma(1, 1)
    # place prior on standard deviation on random year effects
    # page 76 https://nwfsc-timeseries.github.io/atsa-2017/Labs/Week%203%20intro%20to%20jags/intro-to-jags.pdf
    sigma0_seeds[i] ~ dnorm(0, 1) T(0,)
    # calculate precision, which is what JAGS uses
    tau0_seeds[i] <- 1/(sigma0_seeds[i]*sigma0_seeds[i])

    sigma_overdisp_seeds[i] ~ dnorm(0, 1) T(0,)
    tau_overdisp_seeds[i]<-1/( sigma_overdisp_seeds[i]*sigma_overdisp_seeds[i])

    # Proportion of seeds lost to herbivory  -----------------------------------------
    # population-level intercept for each population
    mu0[i] ~  dnorm(0, 1)
    # place prior on standard deviation on random year effects
    sigma0[i] ~ dnorm(0, 1) T(0,)
    # calculate precision, which is what JAGS uses
    tau0[i] <- 1/(sigma0[i]*sigma0[i])
  }

  # indexed by site*year combinations
  for(i in 1:n_siteYearIndex_und){

    # for each site*year combination, use the site (site_observed)
    # to get the population-level intercept and standard deviation of random year effects
    # draw the random year effects
    mu_seeds[i] ~ dlnorm(nu_seeds[site_und_observed[i]], tau0_seeds[site_und_observed[i]])
    mu_log_seeds[i] <- log(mu_seeds[i])

  }

  for(i in 1:n_siteYearIndex_dam){

    # for each site*year combination, use the site (site_observed)
    # to get the population-level intercept and standard deviation of random year effects
    # draw the random year effects
    mu[i] ~ dnorm(mu0[site_dam_observed[i]], tau0[site_dam_observed[i]])

  }

  # LIKELIHOODS -------------------------------------------------------------

  # SEEDS PER UNDAMAGED FRUITS -------------------------------------------------
  for (i in 1:n3){

    eps_tfe_seeds[i] ~ dnorm(0,tau_overdisp_seeds[site3[i]])

    log(z_sd[i]) <- mu_log_seeds[siteYearIndex_und[i]] + eps_tfe_seeds[i]
    sdno[i] ~ dpois(z_sd[i])

    # Posterior predictive and chi-squared
    y_sd.sim[i] ~ dpois(z_sd[i])
    chi2.sd.obs[i] <- pow((sdno[i] - z_sd[i]),2) / (z_sd[i])
    chi2.sd.sim[i] <- pow((y_sd.sim[i]- z_sd[i]),2) / (z_sd[i])
  }

  # # SEEDS PER DAMAGED FRUITS -------------------------------------------------
  for (i in 1:n4){

    eps_seeds[i] ~ dnorm(0,tau_overdisp_seeds[site4[i]])

    # undamaged seeds, average
    log(z_sd_und[i]) <- mu_log_seeds[siteYearIndex_dam[i]]+eps_seeds[i]

    # proportion seeds eaten
    logit(prop_dam[i]) <- mu[siteYearIndex_dam[i]]

    # calculate number of seeds per damaged fruit
    z_sd_dam[i] <- z_sd_und[i]*prop_dam[i]
    sdno_dam[i] ~ dpois(z_sd_dam[i])

    # Posterior predictive and chi-squared
    y_sd_dam.sim[i] ~ dpois(z_sd_dam[i])
    chi2.sd_dam.obs[i] <- pow((sdno_dam[i] - z_sd_dam[i]),2) / (z_sd_dam[i])
    chi2.sd_dam.sim[i] <- pow((y_sd_dam.sim[i]- z_sd_dam[i]),2) / (z_sd_dam[i])
  }


}
