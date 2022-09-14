model {

  # Priors --------------------------------------------------------------
  for(i in 1:n_siteSurvival){

    # seed survival ---------------------------------------------------------------
    for(k in 1:4){

      # WEAKLY INFORMATIVE
      mu0_s[i,k] ~ dnorm(0,1)

      # WEAKLY INFORMATIVE
      sigma0_s[i,k] ~ dnorm(0,1) T(0,)
      tau0_s[i,k] <- 1/(sigma0_s[i,k]*sigma0_s[i,k])

    }

    # age 1 seeds observed in round 1-3
    # s1
    for(j in c(1,7,11)){
      # i indexes site; j indexes year; k indexes germination variable
      mu_s[i,j] ~ dnorm(mu0_s[i,1], tau0_s[i,1])

    }

    # s2
    for(j in c(2,8,12)){
      # i indexes site; j indexes year; k indexes germination variable
      mu_s[i,j] ~ dnorm(mu0_s[i,2], tau0_s[i,2])

    }

    # age 2 seeds observed in round 1-2
    # s3
    for(j in c(3,9)){
      # i indexes site; j indexes year; k indexes germination variable
      mu_s[i,j] ~ dnorm(mu0_s[i,3], tau0_s[i,3])

    }

    # s4
    for(j in c(4,10)){
      # i indexes site; j indexes year; k indexes germination variable
      mu_s[i,j] ~ dnorm(mu0_s[i,4], tau0_s[i,4])

    }

    # age 3 seeds observed in round 1
    # s5
    for(j in 5){
      # i indexes site; j indexes year; k indexes germination variable
      mu_s[i,j] ~ dnorm(0,1)

    }

    for(j in 6){
      # i indexes site; j indexes year; k indexes germination variable
      mu_s[i,j] ~ dnorm(0,1)

    }

    # seed germination ---------------------------------------------------------------
    # population*year parameter, hierarchical logit
    for(k in 1:2){
      mu0_g[i,k] ~ dnorm(0, 1)

      # WEAKLY INFORMATIVE
      sigma0_g[i,k] ~ dnorm(0,1) T(0,)
      tau0_g[i,k] <- 1/(sigma0_g[i,k]*sigma0_g[i,k])
    }

    # age 1 seeds observed in round 1-3
    for(j in c(1,4,6)){
      # i indexes site; j indexes year; k indexes germination variable
      mu_g[i,j] ~ dnorm(mu0_g[i,1], tau0_g[i,1])

    }

    # age 2 seeds observed in round 1-2
    for(j in c(2,5)){
      # i indexes site; j indexes year; k indexes germination variable
      mu_g[i,j] ~ dnorm(mu0_g[i,2], tau0_g[i,2])

    }

    # age 3 seeds observed in round 1
    for(j in 3){
      # i indexes site; j indexes year; k indexes germination variable
      mu_g[i,j] ~ dnorm(0, 1)

    }

    # unobserved s0 ---------------------------------------------------------------
    mu0_s0[i] ~ dnorm(0,1)

    # WEAKLY INFORMATIVE
    sigma0_s0[i] ~ dnorm(0,1) T(0,)
    tau0_s0[i] <- 1/(sigma0_s0[i]*sigma0_s0[i])

    for(j in 1:2){

      mu_s0[i,j] ~ dnorm(mu0_s0[i],tau0_s0[i])


    }

  }

  # Likelihood --------------------------------------------------------------

  # seed survival --------------------------------------------------------------
  for(i in 1:n1){

    # POPULATION*YEAR MARGINAL RATES
    # Round 1: Germination
    logit(theta_g1[i,1]) <- mu_g[siteSurvival[i],1]
    logit(theta_g1[i,2]) <- mu_g[siteSurvival[i],2]
    logit(theta_g1[i,3]) <- mu_g[siteSurvival[i],3]

    # Round 1: Survival
    logit(theta_s[i,1]) <- mu_s[siteSurvival[i],1]
    logit(theta_s[i,2]) <- mu_s[siteSurvival[i],2]
    logit(theta_s[i,3]) <- mu_s[siteSurvival[i],3]
    logit(theta_s[i,4]) <- mu_s[siteSurvival[i],4]
    logit(theta_s[i,5]) <- mu_s[siteSurvival[i],5]
    logit(theta_s[i,6]) <- mu_s[siteSurvival[i],6]

    # Round 2
    logit(theta_g1[i,4]) <- mu_g[siteSurvival[i],4]
    logit(theta_g1[i,5]) <- mu_g[siteSurvival[i],5]

    # survival
    logit(theta_s[i,7]) <- mu_s[siteSurvival[i],7]
    logit(theta_s[i,8]) <- mu_s[siteSurvival[i],8]
    logit(theta_s[i,9]) <- mu_s[siteSurvival[i],9]
    logit(theta_s[i,10]) <- mu_s[siteSurvival[i],10]

    # Round 3
    logit(theta_g1[i,6]) <- mu_g[siteSurvival[i],6]

    # survival
    logit(theta_s[i,11]) <- mu_s[siteSurvival[i],11]
    logit(theta_s[i,12]) <- mu_s[siteSurvival[i],12]

    # COMPOSITE EVENT HISTORIES
    theta_c[i,1] = theta_s[i,1]                  # jan - year/round 1 - age 1
    theta_c[i,2] = theta_s[i,1]*(1-theta_g1[i,1])*theta_s[i,2] # oct - year/round 1 - age 1
    theta_c[i,3] = theta_s[i,1]*(1-theta_g1[i,1])*theta_s[i,2]*theta_s[i,3] # jan - year/round 1 - age 2
    theta_c[i,4] = theta_s[i,1]*(1-theta_g1[i,1])*theta_s[i,2]*theta_s[i,3]*(1-theta_g1[i,2])*theta_s[i,4] # oct - year/round 1 - age 2
    theta_c[i,5] = theta_s[i,1]*(1-theta_g1[i,1])*theta_s[i,2]*theta_s[i,3]*(1-theta_g1[i,2])*theta_s[i,4]*theta_s[i,5] # jan - year/round 1 - age 3
    theta_c[i,6] = theta_s[i,1]*(1-theta_g1[i,1])*theta_s[i,2]*theta_s[i,3]*(1-theta_g1[i,2])*theta_s[i,4]*theta_s[i,5]*(1-theta_g1[i,3])*theta_s[i,6] # oct - year/round 1 - age 3
    theta_c[i,7] = theta_s[i,7]                   # jan - year/round 2 - age 1
    theta_c[i,8] = theta_s[i,7]*(1-theta_g1[i,4])*theta_s[i,8] # oct - year/round 2 - age 1
    theta_c[i,9] = theta_s[i,7]*(1-theta_g1[i,4])*theta_s[i,8]*theta_s[i,9] # jan - year/round 2 - age 2
    theta_c[i,10] = theta_s[i,7]*(1-theta_g1[i,4])*theta_s[i,8]*theta_s[i,9]*(1-theta_g1[i,5])*theta_s[i,10] # oct - year/round 2 - age 2
    theta_c[i,11] = theta_s[i,11]                  # jan - year/round 3 - age 1
    theta_c[i,12] = theta_s[i,11]*(1-theta_g1[i,6])*theta_s[i,12] # oct - year/round 3 - age 1

    # NONPARAMETRIC SURVIVAL MODEL
    mu[i] <- theta_c[i,compIndex[i]]

    ## LIKELIHOOD
    y[i] ~ dbinom(mu[i], seedStart[i])

    ## POSTERIOR PREDICTIVE
    y_sim[i] ~ dbinom(mu[i], seedStart[i])

    # Chi-squared
    chi2.yobs[i] <- pow((y[i]- mu[i]*seedStart[i]),2) / (mu[i]*seedStart[i]+.001)
    chi2.ysim[i] <- pow((y_sim[i]- mu[i]*seedStart[i]),2) / (mu[i]*seedStart[i]+.001)

  }

  # seed germination --------------------------------------------------------------
  for(i in 1:n2){

    logit(g[i]) <- mu_g[siteGermination[i],germinationIndex[i]]

    # LIKELIHOOD
    seedlingJan[i] ~ dbinom(g[i], totalJan[i])

    # POSTERIOR PREDICTIVE
    seedlingJan_sim[i] ~ dbinom(g[i], totalJan[i])

    # Chi-squared
    chi2.obs[i] <- pow((seedlingJan[i]- g[i]*totalJan[i]),2) / (g[i]*totalJan[i]+.001)
    chi2.sim[i] <- pow((seedlingJan_sim[i]- g[i]*totalJan[i]),2) / (g[i]*totalJan[i]+.001)

  }

  # unobserved s0 -------------------------------------------------------------
  for(i in 1:n3){

    logit(s0_1[i]) <- mu_s0[sitePlot[i],1]
    logit(s0_2[i]) <- mu_s0[sitePlot[i],2]

    # year t-2 parameters
    logit(g_4[i]) <- mu_g[sitePlot[i],4]
    logit(g_5[i]) <- mu_g[sitePlot[i],5]
    logit(s_7[i]) <- mu_s[sitePlot[i],7]
    logit(s_8[i]) <- mu_s[sitePlot[i],8]
    logit(s_9[i]) <- mu_s[sitePlot[i],9]

    # year t-1 parameters
    logit(g_6[i]) <- mu_g[sitePlot[i],6]
    logit(s_11[i]) <- mu_s[sitePlot[i],11]

    # LIKELIHOOD
    # year t-1
    y_t1[i] ~ dbinom(s0_2[i]*s_11[i]*g_6[i],n_t1[i])
    # year t-2
    y_t2[i] ~ dbinom(s0_1[i]*s_7[i]*(1-g_4[i])*s_8[i]*s_9[i]*g_5[i],n_t2[i])

    plotSeedlings[i] ~ dsum( y_t1[i], y_t2[i])

    # POSTERIOR PREDICTIVE
    plotSeedlings1_sim[i] ~ dbinom(s0_2[i]*s_11[i]*g_6[i], n_t1[i])
    plotSeedlings2_sim[i] ~ dbinom(s0_1[i]*s_7[i]*(1-g_4[i])*s_8[i]*s_9[i]*g_5[i], n_t2[i])

    plotSeedlings_sim[i] <- plotSeedlings1_sim[i]+plotSeedlings2_sim[i]

    # Chi-squared
    # denominator
    pred_sdlgs[i] <- s0_2[i]*s_11[i]*g_6[i]*n_t1[i]+s0_1[i]*s_7[i]*(1-g_4[i])*s_8[i]*s_9[i]*g_5[i]*n_t2[i]

    chi2.plot.obs[i] <- pow((plotSeedlings[i]- pred_sdlgs[i]),2) / (pred_sdlgs[i]+.001)         # obs.
    chi2.plot.sim[i] <- pow((plotSeedlings_sim[i]- pred_sdlgs[i]),2) / (pred_sdlgs[i]+.001)         # obs.

  }


}
