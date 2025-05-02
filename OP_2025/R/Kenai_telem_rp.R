nfish <- 500
nstations_sep <- c(3,3)
nstations <- sum(nstations_sep)


# define candidate models

library(jagsUI)
library(jagshelper)

# HIDDEN MARKOV MODEL

kenai_jags_hm <- tempfile()
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
    overall_survival[j] <- prod(psurvival[1:j])
  }
}', file=kenai_jags_hm)


# MULTINOMIAL MODEL

kenai_jags_mn <- tempfile()
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
  
}', file=kenai_jags_mn)


#### expansion of all possible capture histories ####

library(magrittr)
possible_histories <- replicate(nstations, 0:1, simplify=FALSE) %>%
  expand.grid %>% as.matrix %>%
  apply(1, paste0, collapse="")


#### expanding an array of all possible states
# dimensions are [fish, station, state]

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



nsim <- 150

# initializing matrices for parameter estimates (medians) & sd's, and correlations
est_detection_hm <- est_survival_hm <- est_overall_survival_hm <- matrix(nrow=nsim, ncol=nstations)
thecors_hm <- array(dim=c(2*nstations, 2*nstations, nsim))
est_detection_hmi <- est_survival_hmi <- est_overall_survival_hmi <- matrix(nrow=nsim, ncol=nstations)
thecors_hmi <- array(dim=c(2*nstations, 2*nstations, nsim))
est_detection_mn <- est_survival_mn <- est_overall_survival_mn <- matrix(nrow=nsim, ncol=nstations)
thecors_mn <- array(dim=c(2*nstations, 2*nstations, nsim))
sd_detection_hm <- sd_survival_hm <- sd_overall_survival_hm <- matrix(nrow=nsim, ncol=nstations)
sd_detection_hmi <- sd_survival_hmi <- sd_overall_survival_hmi <- matrix(nrow=nsim, ncol=nstations)
sd_detection_mn <- sd_survival_mn <- sd_overall_survival_mn <- matrix(nrow=nsim, ncol=nstations)

# initializing matrices for simulated "true" param values for comparison
pdetections <- psurvivals <- overall_survivals <- matrix(nrow=nsim, ncol=nstations)

for(isim in 1:nsim) {
  
  # first, simulating detection & survival probabilities
  # (separately for each simulated replicate)
  
  pdetection <- c(rbeta(nstations_sep[1], 18, 1), 
                  rbeta(nstations_sep[2], 7, 3))
  
  psurvival <- c(rbeta(nstations_sep[1], 18, 1),
                 rbeta(nstations_sep[2], 15, 2))
  
  overall_survival <- sapply(seq_along(psurvival), \(x) prod(psurvival[1:x]))
  
  # saving these
  pdetections[isim, ] <- pdetection
  psurvivals[isim, ] <- psurvival
  overall_survivals[isim, ] <- overall_survival
  
  
  # from probabilities, simulating:
  # X: unobserved state matrix of survival
  # Y: observed detection
  
  X <- Y <- matrix(nrow=nfish, ncol=nstations)
  # X is (unobserved) survival
  # Y is survival & detection
  
  X[,1] <- rbinom(nfish, 1, psurvival[1])
  # X[,1] <- rbinom(nfish, rbinom(nfish, 1, phandling), psurvival[1])
  for(j in 2:nstations) {
    X[,j] <- rbinom(nfish, X[,j-1], psurvival[j])
  }
  
  for(j in 1:nstations) Y[,j] <- rbinom(nfish, X[,j], pdetection[j])
  
  
  # imputing state matrix from observed data where logically possible
  
  Ximputed <- matrix(nrow=nfish, ncol=nstations)
  for(i in 1:nfish) {
    if(any(Y[i,]==1)) {
      Ximputed[i,1:max(which(Y[i,]==1))] <- 1
    }
  }
  
  
  # re-expressing data matrix Y in terms of possible capture history
  
  Yhistories <- apply(Y, 1, paste0, collapse="") %>%
    sapply(\(x) which(x==possible_histories)) %>% unname
  
  
  
  #### running all models
  
  # JAGS controls
  niter <- 4000
  # ncores <- 3
  ncores <- min(10, parallel::detectCores()-1)
  
  
  # RUNNING THE HIDDEN MARKOV MODEL
  
  # bundle data to pass into JAGS
  kenai_data_hm <- list(nfish=nfish,
                        nstations=nstations,
                        Y=Y)
  
  {
    tstart <- Sys.time()
    # print(tstart)
    kenai_jags_out_hm <- jagsUI::jags(model.file=kenai_jags_hm, data=kenai_data_hm,
                                      parameters.to.save=c("psurvival", "pdetection",
                                                           "overall_survival"),
                                      n.chains=ncores, parallel=T, n.iter=niter,
                                      n.burnin=niter/2, n.thin=niter/2000,
                                      verbose=FALSE)
    if(isim==1) print(Sys.time() - tstart)
  }
  
  
  # RUNNING THE SLIGHTLY-LESS-HIDDEN MARKOV MODEL
  
  # bundle data to pass into JAGS
  kenai_data_hmi <- list(nfish=nfish,
                         nstations=nstations,
                         X=Ximputed,
                         Y=Y)
  
  {
    tstart <- Sys.time()
    # print(tstart)
    kenai_jags_out_hmi <- jagsUI::jags(model.file=kenai_jags_hm, data=kenai_data_hmi,
                                       parameters.to.save=c("psurvival", "pdetection",
                                                            "overall_survival"),
                                       n.chains=ncores, parallel=T, n.iter=niter,
                                       n.burnin=niter/2, n.thin=niter/2000,
                                       verbose=FALSE)
    if(isim==1) print(Sys.time() - tstart)
  }
  
  
  # RUNNING THE MULTINOMIAL MODEL
  
  # bundle data to pass into JAGS
  kenai_data_mn <- list(nfish=nfish,
                        nstations=nstations,
                        states=states,
                        dimstates=dim(states),
                        Yhistories=Yhistories)
  
  {
    tstart <- Sys.time()
    # print(tstart)
    kenai_jags_out_mn <- jagsUI::jags(model.file=kenai_jags_mn, data=kenai_data_mn,
                                       parameters.to.save=c("psurvival", "pdetection",
                                                            "overall_survival"),
                                       n.chains=ncores, parallel=T, n.iter=niter,
                                       n.burnin=niter/2, n.thin=niter/2000,
                                       verbose=FALSE)
    if(isim==1) print(Sys.time() - tstart)
  }
  
  # storing all the things
  est_detection_hm[isim, ] <- kenai_jags_out_hm$q50$pdetection
  est_detection_hmi[isim, ] <- kenai_jags_out_hmi$q50$pdetection
  est_detection_mn[isim, ] <- kenai_jags_out_mn$q50$pdetection
  est_survival_hm[isim, ] <- kenai_jags_out_hm$q50$psurvival
  est_survival_hmi[isim, ] <- kenai_jags_out_hmi$q50$psurvival
  est_survival_mn[isim, ] <- kenai_jags_out_mn$q50$psurvival
  est_overall_survival_hm[isim, ] <- kenai_jags_out_hm$q50$overall_survival
  est_overall_survival_hmi[isim, ] <- kenai_jags_out_hmi$q50$overall_survival
  est_overall_survival_mn[isim, ] <- kenai_jags_out_mn$q50$overall_survival
  sd_detection_hm[isim, ] <- kenai_jags_out_hm$sd$pdetection
  sd_detection_hmi[isim, ] <- kenai_jags_out_hmi$sd$pdetection
  sd_detection_mn[isim, ] <- kenai_jags_out_mn$sd$pdetection
  sd_survival_hm[isim, ] <- kenai_jags_out_hm$sd$psurvival
  sd_survival_hmi[isim, ] <- kenai_jags_out_hmi$sd$psurvival
  sd_survival_mn[isim, ] <- kenai_jags_out_mn$sd$psurvival
  sd_overall_survival_hm[isim, ] <- kenai_jags_out_hm$sd$overall_survival
  sd_overall_survival_hmi[isim, ] <- kenai_jags_out_hmi$sd$overall_survival
  sd_overall_survival_mn[isim, ] <- kenai_jags_out_mn$sd$overall_survival
  thecors_hm[,,isim] <- cor(cbind(kenai_jags_out_hm$sims.list$psurvival, kenai_jags_out_hm$sims.list$pdetection))
  thecors_hmi[,,isim] <- cor(cbind(kenai_jags_out_hmi$sims.list$psurvival, kenai_jags_out_hmi$sims.list$pdetection))
  thecors_mn[,,isim] <- cor(cbind(kenai_jags_out_mn$sims.list$psurvival, kenai_jags_out_mn$sims.list$pdetection))
  
}