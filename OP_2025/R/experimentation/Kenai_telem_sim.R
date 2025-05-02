# The purpose of this script was originally to create a Hidden Markov model
# to estimate the respective survival and detection probabilities corresponding
# to a set of receiver arrays (stations) in the Kenai river.  The model file has
# quite a lot of unnecessary components and its current state represents several
# avenues of experimentation.

# After developing the HM model, this script was then repurposed for comparison
# between three candidate models, by means of meta-simulation.  These were:
# * The original Hidden Markov model
# * the HM model in which survival states were logically imputed when possible
# * the Multinomial model developed in R/experimentation/multinomial_model.R

# Surprisingly, comparison showed no advantage to the multinomial model, and the
# HMimputed model (#2) ran in a fraction of the time of the others.



# This script is structured as follows:
# * Defining simulation parameters (number of fish, number of stations, survival 
#   and detection probabilities)
# * Then, for (nsim) number of replicates:
#   - simulate a matrix of survival states (X) and observed data (Y)
#   - run all three candidate models
#   - store posterior medians, sd's, and within-parameter correlations for all models
# * Finally, compare all estimates, sd's, and correlations.



# A few additional statistical notes:
# * The vector of survival probabilities are defined as conditional on survival 
#   at the time of the previous event phi_j = p(survival_j | survival_(j-1))
#
# * The vector of detection probabilities are defined as conditional on 
#   survival at the time of the same event p_j = p(detection_j | survival_j)
#
# * The state matrix X, used in simulation and also within the Bayesian model,
#   is defined with a row for each fish and a column for each event (station),
#   and has values of 1 or 0 depending on whether a fish is ALIVE.
#
# * The data matrix Y, which is simulated and passed as data into the Bayesian
#   model, is defined similarly (a row for each fish and a column for each event)
#   with values of 1 or 0 depending on whether a fish is DETECTED.
#
# * Some logical imputation of the state matrix X is possible, since if a fish
#   is observed in event j, it must be alive at event j and all events preceding.
#   Defining an imputed state matrix (1: known survival, NA: unknown) and passing
#   this as a data object to JAGS will allow it to treat known survival as data
#   and NA values as unknown parameters in the HM model.  Functionally, this 
#   results in equivalent inferences (thankfully) and much shorter run-time due
#   to fewer unknown parameters.



nfish <- 500
nstations_sep <- c(3,3)

pdetection <- c(rbeta(nstations_sep[1], 18, 1), 
                rbeta(nstations_sep[2], 7, 3))

psurvival <- c(rbeta(nstations_sep[1], 18, 1),
               rbeta(nstations_sep[2], 15, 2))

# phandling <- 0.9

overall_survival <- sapply(seq_along(psurvival), \(x) prod(psurvival[1:x]))
nstations <- sum(nstations_sep)

               
               
nsim <- 150
est_detection_hm <- est_survival_hm <- est_overall_survival_hm <- matrix(nrow=nsim, ncol=nstations)
thecors_hm <- array(dim=c(2*nstations, 2*nstations, nsim))
est_detection_hmi <- est_survival_hmi <- est_overall_survival_hmi <- matrix(nrow=nsim, ncol=nstations)
thecors_hmi <- array(dim=c(2*nstations, 2*nstations, nsim))
est_detection_mn <- est_survival_mn <- est_overall_survival_mn <- matrix(nrow=nsim, ncol=nstations)
thecors_mn <- array(dim=c(2*nstations, 2*nstations, nsim))
sd_detection_hm <- sd_survival_hm <- sd_overall_survival_hm <- matrix(nrow=nsim, ncol=nstations)
sd_detection_hmi <- sd_survival_hmi <- sd_overall_survival_hmi <- matrix(nrow=nsim, ncol=nstations)
sd_detection_mn <- sd_survival_mn <- sd_overall_survival_mn <- matrix(nrow=nsim, ncol=nstations)
for(isim in 1:nsim) {
               
  # psurvival[nstations] <- 1
  
  X <- Y <- matrix(nrow=nfish, ncol=nstations)
  # X is (unobserved) survival
  # Y is survival & detection
  
  X[,1] <- rbinom(nfish, 1, psurvival[1])
  # X[,1] <- rbinom(nfish, rbinom(nfish, 1, phandling), psurvival[1])
  for(j in 2:nstations) {
    X[,j] <- rbinom(nfish, X[,j-1], psurvival[j])
  }
  
  for(j in 1:nstations) Y[,j] <- rbinom(nfish, X[,j], pdetection[j])
  
  Ximputed <- matrix(nrow=nfish, ncol=nstations)
  for(i in 1:nfish) {
    if(any(Y[i,]==1)) {
      Ximputed[i,1:max(which(Y[i,]==1))] <- 1
    }
  }
  
  colMeans(X)
  colMeans(Y)
  overall_survival
  
  # prod(psurvival)
  # mean(X[,nstations])
  
  
  library(jagshelper)
  library(jagsUI)
  
  # specify model, which is written to a temporary file
  kenai_jags <- tempfile()
  cat('model {
  for(i in 1:nfish) {
    X[i,1] ~ dbinom(psurvival[1], 1)
    Y[i,1] ~ dbinom(pdetection[1], X[i,1])
    for(j in 2:nstations) {
      X[i,j] ~ dbinom(psurvival[j], X[i,j-1])
      Y[i,j] ~ dbinom(pdetection[j], X[i,j])
    }
  }
  
  for(j in 1:nstations) {
    psurvival[j] ~ dbeta(5,1) #dnorm(mu_survival, tau_survival)T(0,1) #
    pdetection[j] ~ dbeta(5,1) #dnorm(mu_detection, tau_detection)T(0,1) #
    
    # logit(psurvival[j]) <- expit_psurvival[j] 
    # logit(pdetection[j]) <- expit_pdetection[j]
    # expit_psurvival[j] ~ dnorm(mu_survival, tau_survival)
    # expit_pdetection[j] ~ dnorm(mu_detection, tau_detection)
    
    # psurvival[j] ~ dbeta(asurvival, bsurvival)
    # pdetection[j] ~ dbeta(adetection, bdetection)
    
    overall_survival[j] <- prod(psurvival[1:j])
  }
  # psurvival[nstations] <- 1
  # pdetection[nstations] ~ dbeta(9,1)
  
  mu_survival ~ dnorm(1, 2) # dnorm(0.9, 1)
  mu_detection ~ dnorm(1, 2) # dnorm(0.9, 1)
  tau_survival <- pow(sig_survival, -2)
  tau_detection <- pow(sig_detection, -2)
  sig_survival ~ dunif(0, 5) # dunif(0, 0.1)
  sig_detection ~ dunif(0, 5) # dunif(0, 0.1)
  
  asurvival <- Msurvival*Dsurvival
  bsurvival <- Dsurvival*(1-Msurvival)
  Msurvival ~ dbeta(9,1)
  Dsurvival ~ dnorm(10,1)
  adetection <- Mdetection*Ddetection
  bdetection <- Ddetection*(1-Mdetection)
  Mdetection ~ dbeta(9,1)
  Ddetection ~ dnorm(10,1)
  
  # overall_survival <- prod(psurvival)
}', file=kenai_jags)
  
  
  
  # bundle data to pass into JAGS
  kenai_data <- list(nfish=nfish,
                     nstations=nstations,
                     # X=Ximputed,
                     Y=Y)
  
  # JAGS controls
  niter <- 4000
  # ncores <- 3
  ncores <- min(10, parallel::detectCores()-1)
  
  {
    tstart <- Sys.time()
    print(tstart)
    kenai_jags_out <- jagsUI::jags(model.file=kenai_jags, data=kenai_data,
                                   parameters.to.save=c("psurvival", "pdetection",
                                                        # "mu_survival", "mu_detection",
                                                        # "sig_survival", "sig_detection",
                                                        "overall_survival"),
                                   n.chains=ncores, parallel=T, n.iter=niter,
                                   n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
  }
  
  # nbyname(kenai_jags_out)
  # plotRhats(kenai_jags_out)
  # traceworstRhat(kenai_jags_out, parmfrow = c(3,3))
  # 
  
  # par(mfrow=c(1,3))
  # caterpillar(kenai_jags_out, p="pdetection", ylim=0:1)
  # points(pdetection)
  # caterpillar(kenai_jags_out, p="psurvival", ylim=0:1)
  # points(psurvival)
  # caterpillar(kenai_jags_out, p="overall_survival")
  # points(overall_survival)
  
  # caterpillar(kenai_jags_out, p="mu_survival")
  # caterpillar(kenai_jags_out, p="mu_detection")
  # caterpillar(kenai_jags_out, p="sig_survival")
  # caterpillar(kenai_jags_out, p="sig_detection")
  # 
  # 
  # plotcor_jags(kenai_jags_out, p=c("psurvival", "pdetection"), maxn=100)
  
  est_detection_hm[isim,] <- kenai_jags_out$q50$pdetection
  est_survival_hm[isim,] <- kenai_jags_out$q50$psurvival
  est_overall_survival_hm[isim,] <- kenai_jags_out$q50$overall_survival
  sd_detection_hm[isim,] <- kenai_jags_out$sd$pdetection
  sd_survival_hm[isim,] <- kenai_jags_out$sd$psurvival
  sd_overall_survival_hm[isim,] <- kenai_jags_out$sd$overall_survival
  thecors_hm[,,isim] <- cor(cbind(kenai_jags_out$sims.list$psurvival, kenai_jags_out$sims.list$pdetection))
  
  
  
  # bundle data to pass into JAGS
  kenai_data <- list(nfish=nfish,
                     nstations=nstations,
                     X=Ximputed,
                     Y=Y)
  
  # JAGS controls
  niter <- 4000
  # ncores <- 3
  ncores <- min(10, parallel::detectCores()-1)
  
  {
    tstart <- Sys.time()
    print(tstart)
    kenai_jags_out <- jagsUI::jags(model.file=kenai_jags, data=kenai_data,
                                   parameters.to.save=c("psurvival", "pdetection",
                                                        # "mu_survival", "mu_detection",
                                                        # "sig_survival", "sig_detection",
                                                        "overall_survival"),
                                   n.chains=ncores, parallel=T, n.iter=niter,
                                   n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
  }
  
  est_detection_hmi[isim,] <- kenai_jags_out$q50$pdetection
  est_survival_hmi[isim,] <- kenai_jags_out$q50$psurvival
  est_overall_survival_hmi[isim,] <- kenai_jags_out$q50$overall_survival
  sd_detection_hmi[isim,] <- kenai_jags_out$sd$pdetection
  sd_survival_hmi[isim,] <- kenai_jags_out$sd$psurvival
  sd_overall_survival_hmi[isim,] <- kenai_jags_out$sd$overall_survival
  thecors_hmi[,,isim] <- cor(cbind(kenai_jags_out$sims.list$psurvival, kenai_jags_out$sims.list$pdetection))
  
  
  
  # par(mfrow=c(1,3))
  # caterpillar(kenai_jags_out, p="pdetection", ylim=0:1)
  # points(pdetection)
  # caterpillar(kenai_jags_out, p="psurvival", ylim=0:1)
  # points(psurvival)
  # caterpillar(kenai_jags_out, p="overall_survival")
  # points(overall_survival)
  
  
  
  
  
  library(magrittr)
  possible_histories <- replicate(nstations, 0:1, simplify=FALSE) %>%
    expand.grid %>% as.matrix %>%
    apply(1, paste0, collapse="")
  
  Yhistories <- apply(Y, 1, paste0, collapse="") %>%
    sapply(\(x) which(x==possible_histories)) %>% unname
  
  
  empirical_props <- table(factor(Yhistories, levels=seq(length(possible_histories)))) %>%
    as.numeric/nfish
  
  
  history_mat <- replicate(nstations, 0:1, simplify=FALSE) %>%
    expand.grid %>% as.matrix
  alive <- array(dim=c(dim(history_mat), ncol(history_mat)+1))
  nstates <- rep(1, nrow(history_mat))
  alive[,,1] <- 1
  for(i in 1:nrow(history_mat)) {
    for(j in 1:ncol(history_mat)) {
      if(all(history_mat[i, (ncol(history_mat)-j+1):ncol(history_mat)] == 0)) {
        alive[i, , j+1] <- c(rep(1, ncol(history_mat)-j), rep(0, j))
        nstates[i] <- j+1
      }
    }
  }
  
  # define:
  # 1 = observed (alive)
  # 2 = unobserved (alive)
  # 3 = newly dead
  # 4 = dead
  # 5 = not considered
  states <- array(dim=c(dim(history_mat), ncol(history_mat)+1))
  # states[history_mat==1] <- 1
  # states[history_mat==0 & alive==1] <- 2
  # states[is.na(alive)] <- 4
  
  for(i in 1:dim(states)[1]) {
    for(j in 1:dim(states)[2]) {
      for(k in 1:dim(states)[3]) {
        if(is.na(alive[i,j,k])) {
          states[i,j,k] <- 5
        } else {
          if(alive[i,j,k] & history_mat[i,j]) {
            states[i,j,k] <- 1
          }
          if(alive[i,j,k] & !history_mat[i,j]) {
            states[i,j,k] <- 2
          }
          if(all(!alive[i,j:dim(states)[2],k])) {
            states[i,j,k] <- 3
          } 
          if(j>1) {
            if(!alive[i,j-1,k]) {
              states[i,j,k] <- 4
            }
          }
        }
      }
    }
  }
  
  
  library(jagshelper)
  library(jagsUI)
  
  # specify model, which is written to a temporary file
  kenai_jags <- tempfile()
  cat('model {
  for(i in 1:nfish) {
    Yhistories[i] ~ dcat(pi[])
  }
  for(i in 1:dimstates[1]) {
    for(k in 1:dimstates[3]) {
      for(j in 1:dimstates[2]) {
        pi_arr[i,j,k] <- (states[i,j,k]==1)*psurvival[j]*pdetection[j] +
                         (states[i,j,k]==2)*psurvival[j]*(1-pdetection[j]) +
                         (states[i,j,k]==3)*(1-psurvival[j]) +
                         (states[i,j,k]==4)*1 +
                         (states[i,j,k]==5)*0
      }
      pi_mat[i,k] <- prod(pi_arr[i,,k])
    }
    pi[i] <- sum(pi_mat[i,])
  }
  
  for(jstation in 1:nstations) {
    psurvival[jstation] ~ dbeta(5,1)
    pdetection[jstation] ~ dbeta(5,1)
    overall_survival[jstation] <- prod(psurvival[1:jstation])
  }
  
}', file=kenai_jags)
  
  # define:
  # 1 = observed (alive)
  # 2 = unobserved (alive)
  # 3 = newly dead
  # 4 = dead
  # 5 = not considered
  
  
  # bundle data to pass into JAGS
  kenai_data <- list(nfish=nfish,
                     nstations=nstations,
                     states=states,
                     dimstates=dim(states),
                     Yhistories=Yhistories)
  
  # JAGS controls
  niter <- 4000
  # ncores <- 3
  ncores <- min(10, parallel::detectCores()-1)
  
  {
    tstart <- Sys.time()
    print(tstart)
    kenai_jags_out <- jagsUI::jags(model.file=kenai_jags, data=kenai_data,
                                   parameters.to.save=c("psurvival", "pdetection",
                                                        "pi","overall_survival"),
                                   n.chains=ncores, parallel=T, n.iter=niter,
                                   n.burnin=niter/2, n.thin=niter/2000)
    print(Sys.time() - tstart)
  }
  
  
  # par(mfrow=c(1,3))
  # caterpillar(kenai_jags_out, p="pdetection", ylim=0:1)
  # points(pdetection)
  # caterpillar(kenai_jags_out, p="psurvival", ylim=0:1)
  # points(psurvival)
  # caterpillar(kenai_jags_out, p="overall_survival")
  # points(overall_survival)
  
  est_detection_mn[isim,] <- kenai_jags_out$q50$pdetection
  est_survival_mn[isim,] <- kenai_jags_out$q50$psurvival
  est_overall_survival_mn[isim,] <- kenai_jags_out$q50$overall_survival
  sd_detection_mn[isim,] <- kenai_jags_out$sd$pdetection
  sd_survival_mn[isim,] <- kenai_jags_out$sd$psurvival
  sd_overall_survival_mn[isim,] <- kenai_jags_out$sd$overall_survival
  thecors_mn[,,isim] <- cor(cbind(kenai_jags_out$sims.list$psurvival, kenai_jags_out$sims.list$pdetection))
  
  
  print(isim)
}

par(mfrow=c(1,3))
caterpillar(est_detection_hm, col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Detection probability")
points(pdetection)
legend("bottomleft", legend=c("true","est (inriver)", "est (nearshore)"), lwd=c(NA,3,3), col=c(1,4,2), pch=c(1,NA,NA))

caterpillar(est_detection_hmi, col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Detection probability")
points(pdetection)
legend("bottomleft", legend=c("true","est (inriver)", "est (nearshore)"), lwd=c(NA,3,3), col=c(1,4,2), pch=c(1,NA,NA))

caterpillar(est_detection_mn, col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Detection probability")
points(pdetection)
legend("bottomleft", legend=c("true","est (inriver)", "est (nearshore)"), lwd=c(NA,3,3), col=c(1,4,2), pch=c(1,NA,NA))


caterpillar(est_survival_hm, , col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Survival probability per event")
points(psurvival)
# legend("bottomleft", legend=c("true","est (inriver)", "est (nearshore)"), lwd=c(NA,3,3), col=c(1,4,2), pch=c(1,NA,NA))

caterpillar(est_survival_hmi, , col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Survival probability per event")
points(psurvival)
# legend("bottomleft", legend=c("true","est (inriver)", "est (nearshore)"), lwd=c(NA,3,3), col=c(1,4,2), pch=c(1,NA,NA))

caterpillar(est_survival_mn, , col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Survival probability per event")
points(psurvival)


caterpillar(est_overall_survival_hm, , col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Cumulative survival probability")
points(overall_survival)
# legend("bottomleft", legend=c("true","est (inriver)", "est (nearshore)"), lwd=c(NA,3,3), col=c(1,4,2), pch=c(1,NA,NA))

caterpillar(est_overall_survival_hmi, , col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Cumulative survival probability")
points(overall_survival)
# legend("bottomleft", legend=c("true","est (inriver)", "est (nearshore)"), lwd=c(NA,3,3), col=c(1,4,2), pch=c(1,NA,NA))

caterpillar(est_overall_survival_mn, , col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Cumulative survival probability")
points(overall_survival)



par(mfrow=c(1,3))
caterpillar(sd_detection_hm, col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Detection probability")
caterpillar(sd_detection_hmi, col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Detection probability")
caterpillar(sd_detection_mn, col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Detection probability")

caterpillar(sd_survival_hm, , col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Survival probability per event")
caterpillar(sd_survival_hmi, , col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Survival probability per event")
caterpillar(sd_survival_mn, , col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Survival probability per event")


caterpillar(sd_overall_survival_hm, , col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Cumulative survival probability")
caterpillar(sd_overall_survival_hmi, , col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Cumulative survival probability")
caterpillar(sd_overall_survival_mn, , col=c(rep(4, nstations_sep[1]), rep(2, nstations_sep[2])),
            main="Cumulative survival probability")

par(mfrow=c(1,3))
medcors_hm <- apply(thecors_hm, 1:2, median)
qlocors_hm <- apply(thecors_hm, 1:2, quantile, p=0.05)
qhicors_hm <- apply(thecors_hm, 1:2, quantile, p=0.95)
rownames(medcors_hm) <- colnames(medcors_hm) <- # c(paste("surv", 1:nstations), paste("det", 1:nstations))
  rownames(qlocors_hm) <- colnames(qlocors_hm) <- # c(paste("surv", 1:nstations), paste("det", 1:nstations))
  rownames(qhicors_hm) <- colnames(qhicors_hm) <- c(paste("surv", 1:nstations), paste("det", 1:nstations))
dsftools::plotcor(medcors_hm)
dsftools::plotcor(qlocors_hm)
dsftools::plotcor(qhicors_hm)

medcors_hmi <- apply(thecors_hmi, 1:2, median)
qlocors_hmi <- apply(thecors_hmi, 1:2, quantile, p=0.05)
qhicors_hmi <- apply(thecors_hmi, 1:2, quantile, p=0.95)
rownames(medcors_hmi) <- colnames(medcors_hmi) <- # c(paste("surv", 1:nstations), paste("det", 1:nstations))
  rownames(qlocors_hmi) <- colnames(qlocors_hmi) <- # c(paste("surv", 1:nstations), paste("det", 1:nstations))
  rownames(qhicors_hmi) <- colnames(qhicors_hmi) <- c(paste("surv", 1:nstations), paste("det", 1:nstations))
dsftools::plotcor(medcors_hmi)
dsftools::plotcor(qlocors_hmi)
dsftools::plotcor(qhicors_hmi)

medcors_mn <- apply(thecors_mn, 1:2, median)
qlocors_mn <- apply(thecors_mn, 1:2, quantile, p=0.05)
qhicors_mn <- apply(thecors_mn, 1:2, quantile, p=0.95)
rownames(medcors_mn) <- colnames(medcors_mn) <- # c(paste("surv", 1:nstations), paste("det", 1:nstations))
  rownames(qlocors_mn) <- colnames(qlocors_mn) <- # c(paste("surv", 1:nstations), paste("det", 1:nstations))
  rownames(qhicors_mn) <- colnames(qhicors_mn) <- c(paste("surv", 1:nstations), paste("det", 1:nstations))
dsftools::plotcor(medcors_mn)
dsftools::plotcor(qlocors_mn)
dsftools::plotcor(qhicors_mn)

rp_survival_hm <- rp_detection_hm <- rp_overall_survival_hm <- NA*psurvival
for(j in seq_along(psurvival)) {
  rp_survival_hm[j] <- dsftools::rp(est_survival_hm[,j], psurvival[j], confidence=0.95)
  rp_overall_survival_hm[j] <- dsftools::rp(est_overall_survival_hm[,j], overall_survival[j], confidence=0.95)
  rp_detection_hm[j] <- dsftools::rp(est_detection_hm[,j], pdetection[j], confidence=0.95)
}

rp_survival_hmi <- rp_detection_hmi <- rp_overall_survival_hmi <- NA*psurvival
for(j in seq_along(psurvival)) {
  rp_survival_hmi[j] <- dsftools::rp(est_survival_hmi[,j], psurvival[j], confidence=0.95)
  rp_overall_survival_hmi[j] <- dsftools::rp(est_overall_survival_hmi[,j], overall_survival[j], confidence=0.95)
  rp_detection_hmi[j] <- dsftools::rp(est_detection_hmi[,j], pdetection[j], confidence=0.95)
}

rp_survival_mn <- rp_detection_mn <- rp_overall_survival_mn <- NA*psurvival
for(j in seq_along(psurvival)) {
  rp_survival_mn[j] <- dsftools::rp(est_survival_mn[,j], psurvival[j], confidence=0.95)
  rp_overall_survival_mn[j] <- dsftools::rp(est_overall_survival_mn[,j], overall_survival[j], confidence=0.95)
  rp_detection_mn[j] <- dsftools::rp(est_detection_mn[,j], pdetection[j], confidence=0.95)
}

plotthemboth <- function(x,y,z) {
  plot(x, pch=1, col=2, ylim=range(0,x,y,z))
  points(y, pch=1, col=4)
  points(z, pch=1, col=3)
}
plotthemboth(rp_detection_hm, rp_detection_hmi, rp_detection_mn)
plotthemboth(rp_survival_hm, rp_survival_hmi, rp_survival_mn)
plotthemboth(rp_overall_survival_hm, rp_overall_survival_hmi, rp_overall_survival_mn)
