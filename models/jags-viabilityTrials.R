model {

  # viability ---------------------------------------------------------------
  # i indexes site
  for(i in 1:n_siteViab){

    # hyperpriors -------------------------------------------------------------
    # k indexes seed age
    for(k in 1:2){

      # theta g
      mu0_g[i,k] ~  dnorm(0, 1)
      sigma0_g[i,k] ~ dnorm(0, 1) T(0,)
      tau0_g[i,k] <- 1/(sigma0_g[i,k]*sigma0_g[i,k])

      # theta v
      mu0_v[i,k] ~  dnorm(0, 1)
      sigma0_v[i,k] ~ dnorm(0, 1) T(0,)
      tau0_v[i,k] <- 1/(sigma0_v[i,k]*sigma0_v[i,k])

    }

    # priors ------------------------------------------------------------------
    # j indexes experimental round
    # e.g. age 1 seeds were observed in rounds 1-3
    # mu_g[siteViab[i],yearViab[i],ageViab[i]]
      for(j in 1:3){

        # theta g
        mu_g[i,j] ~ dnorm(mu0_g[i,1], tau0_g[i,1])

        # theta v
        mu_v[i,j] ~ dnorm(mu0_v[i,1], tau0_v[i,1])

      }

    # e.g. age 2 seeds were observed in rounds 1-2
    for(j in 4:5){

      # priors ------------------------------------------------------------------

      # theta g
      mu_g[i,j] ~ dnorm(mu0_g[i,2], tau0_g[i,2])

      # theta v
      mu_v[i,j] ~ dnorm(mu0_v[i,2], tau0_v[i,2])
    }

    # e.g. age 3 seeds were observed in rounds 1
    for(j in 6){

      # priors ------------------------------------------------------------------

      # theta g
      mu_g[i,j] ~  dnorm(0, 1)

      # theta v
      mu_v[i,j] ~  dnorm(0, 1)

    }
  }


  # likelihood (germination trials) --------------------------------------------------------------

  for(i in 1:n){

    # logit
    logit(theta_g[i]) <- mu_g[siteViab[i],indexViab[i]]

    # Likelihood
    germCount[i] ~ dbinom( theta_g[i] , germStart[i] )

    # POSTERIOR PREDICTIVE
    germCount_sim[i] ~ dbinom( theta_g[i] , germStart[i] )

    # Chi-squared
    chi2.germCount.obs[i] <- pow((germCount[i]- theta_g[i]*germStart[i]),2) / (theta_g[i]*germStart[i]+.001)
    chi2.germCount.sim[i] <- pow((germCount_sim[i]- theta_g[i]*germStart[i]),2) / (theta_g[i]*germStart[i]+.001)

  }

  for(i in 1:n_v){

    # logit
    logit(theta_v[i]) <- mu_v[siteViab_v[i],indexViab_v[i]]

    # Likelihood
    viabStain[i] ~ dbinom( theta_v[i] , viabStart[i] )

    # POSTERIOR PREDICTIVE
    viabStain_sim[i] ~ dbinom( theta_v[i] , viabStart[i] )

    # Chi-squared
    chi2.viabStain.obs[i] <- pow((viabStain[i]- theta_v[i]*viabStart[i]),2) / (theta_v[i]*viabStart[i]+.001)
    chi2.viabStain.sim[i] <- pow((viabStain_sim[i]- theta_v[i]*viabStart[i]),2) / (theta_v[i]*viabStart[i]+.001)

  }

}
