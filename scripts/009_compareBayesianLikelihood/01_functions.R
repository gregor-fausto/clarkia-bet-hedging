# - Functions for estimating seed bank transitions ----

# - + Likelihood function ----
# - ++ seed bag burial experiment: 2005-2008 ----

# separate estimates for all ages/rounds/seasons
# g1-6: probability for each age, round
# s1-12: probability for each age, round, season

likSeedBagModel = function(params, dat){
  g.prob = params[1:6]
  s.prob = params[7:18]

  n = dat[[1]]
  y = dat[[2]]
  n.g = dat[[3]]
  y.g = dat[[4]]

  c_ind = dat[[5]]
  g_ind = dat[[6]]

  ll.sur = c()
  for(i in 1:12){
    n.tmp = n[c_ind==i]
    y.tmp = y[c_ind==i]
    n.tmp = n.tmp[!is.na(y.tmp)]
    y.tmp = y.tmp[!is.na(y.tmp)]

    if(i == c(1)){
      prob.s = s.prob[1]
    } else if(i == 2){
      prob.s = s.prob[1]*(1-g.prob[1])*s.prob[2]
    } else if(i == 3){
      prob.s = s.prob[1]*(1-g.prob[1])*s.prob[2]*s.prob[3]
    } else if(i == 4){
      prob.s = s.prob[1]*(1-g.prob[1])*s.prob[2]*s.prob[3]*(1-g.prob[2])*s.prob[4]
    } else if(i == 5){
      prob.s = s.prob[1]*(1-g.prob[1])*s.prob[2]*s.prob[3]*(1-g.prob[2])*s.prob[4]*s.prob[5]
    } else if(i == 6){
      prob.s = s.prob[1]*(1-g.prob[1])*s.prob[2]*s.prob[3]*(1-g.prob[2])*s.prob[4]*s.prob[5]*(1-g.prob[3])*s.prob[6]
    } else if(i == 7){
      prob.s = s.prob[7]
    } else if(i == 8){
      prob.s = s.prob[7]*(1-g.prob[4])*s.prob[8]
    } else if(i == 9){
      prob.s = s.prob[7]*(1-g.prob[4])*s.prob[8]*s.prob[9]
    } else if(i == 10){
      prob.s = s.prob[7]*(1-g.prob[4])*s.prob[8]*s.prob[9]*(1-g.prob[5])*s.prob[10]
    } else if(i == c(11)){
      prob.s = s.prob[11]
    } else if(i == 12){
      prob.s = s.prob[11]*(1-g.prob[6])*s.prob[12]
    }

    ll.sur[i] = sum(dbinom(x=y.tmp, size=n.tmp, prob= prob.s, log=TRUE))
  }

  ll.ger = c()
  for(i in 1:6){
    n.tmp = n.g[g_ind==i]
    y.tmp = y.g[g_ind==i]
    n.tmp = n.tmp[!is.na(y.tmp)]
    y.tmp = y.tmp[!is.na(y.tmp)]

    prob.g = g.prob[i]

    ll.ger[i] = sum(dbinom(x=y.tmp, size=n.tmp, prob=prob.g, log=TRUE))
  }

  tmp = -sum(ll.sur,ll.ger )
  return(tmp)
}

# - Functions for estimating viabilities ----


# - + Likelihood function ----

likViabilityModel = function(params, dat){
  g.prob = params[1:6]
  v.prob = params[7:12]

  n.g = dat[[1]]
  y.g = dat[[2]]
  n.v = dat[[3]]
  y.v = dat[[4]]

  v_ind = dat[[5]]
  v_ind2 = dat[[6]]


  ll.germ = c()
  for(i in 1:6){
    n.tmp = n.g[v_ind==i]
    y.tmp = y.g[v_ind==i]
    n.tmp = n.tmp[!is.na(y.tmp)]
    y.tmp = y.tmp[!is.na(y.tmp)]

    prob.g = g.prob[i]

    ll.germ[i] = sum(dbinom(x=y.tmp, size=n.tmp, prob = prob.g, log=TRUE))
  }

  ll.viab = c()
  for(i in 1:6){
    n.tmp = n.v[v_ind2==i]
    y.tmp = y.v[v_ind2==i]
    n.tmp = n.tmp[!is.na(y.tmp)]
    y.tmp = y.tmp[!is.na(y.tmp)]

    prob.v = v.prob[i]

    ll.viab[i] = sum(dbinom(x=y.tmp, size=n.tmp, prob = prob.v, log=TRUE))
  }

  tmp = -sum(ll.germ, ll.viab )
  return(tmp)
}
