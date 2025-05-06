# simulation controls
nfish <- 500
nstations_sep <- c(3,3) #c(3,3)
nstations <- sum(nstations_sep)



pdetection_river_expected <- 0.95
pdetection_river_range <- seq(from=0.85, to=1, by=0.025)
pdetection_nearshore_expected <- 0.8
pdetection_nearshore_range <- seq(from=0.5, to=0.95, by=0.05)

psurvival_river_expected <- 0.95
psurvival_river_range <- seq(from=0.85, to=1, by=0.025)
psurvival_nearshore_expected <- 0.9
psurvival_nearshore_range <- seq(from=0.5, to=0.95, by=0.05)



nsim <- 3#00   # --------- 50 sim in 1.5 hrs on laptop
save_results <- TRUE



# JAGS controls
niter <- 6000   
ncores <- 6
# ncores <- min(10, parallel::detectCores()-1)



# define candidate models

library(jagsUI)
library(jagshelper)

# HIDDEN MARKOV MODEL   - this is actually the only one!

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


{
  t_overall_start <- Sys.time()
  
  est_detection_bynearshore <- 
    est_survival_bynearshore <- 
    est_overall_survival_bynearshore <- 
    sd_detection_bynearshore <- 
    sd_survival_bynearshore <- 
    sd_overall_survival_bynearshore <- 
    array(dim = c(length(pdetection_nearshore_range),
                  length(psurvival_nearshore_range),
                  nstations,
                  nsim))
  
  for(idetection in 1:length(pdetection_nearshore_range)) {
    for(isurvival in 1:length(psurvival_nearshore_range)) {
      
      pdetection <- c(rep(pdetection_river_expected, nstations_sep[1]),
                      rep(pdetection_nearshore_range[idetection], nstations_sep[2]))
      psurvival <- c(rep(psurvival_river_expected, nstations_sep[1]),
                     rep(psurvival_nearshore_range[isurvival], nstations_sep[2]))
      
      for(isim in 1:nsim) {
        
        
        # overall_survival <- sapply(seq_along(psurvival), \(x) prod(psurvival[1:x]))
        
        
        
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
          if(idetection==1 & isurvival==1) print(Sys.time() - tstart)
        }
        
        
        
        # storing all the things
        est_detection_bynearshore[idetection, isurvival, , nsim] <-
          kenai_jags_out_hmi$q50$pdetection
        est_survival_bynearshore[idetection, isurvival, , nsim] <-
          kenai_jags_out_hmi$q50$psurvival
        sd_detection_bynearshore[idetection, isurvival, , nsim] <-
          kenai_jags_out_hmi$sd$pdetection
        sd_survival_bynearshore[idetection, isurvival, , nsim] <-
          kenai_jags_out_hmi$sd$psurvival
        
        # print(isim)
      }
    }
    print(idetection)
  }
  
  
  for(idetection in 1:length(pdetection_river_range)) {
    for(isurvival in 1:length(psurvival_river_range)) {
      
      pdetection <- c(rep(pdetection_river_range[idetection], nstations_sep[1]),
                      rep(pdetection_nearshore_expected, nstations_sep[2]))
      psurvival <- c(rep(psurvival_river_range[isurvival], nstations_sep[1]),
                     rep(psurvival_nearshore_expected, nstations_sep[2]))
      
      for(isim in 1:nsim) {
        
        
        # overall_survival <- sapply(seq_along(psurvival), \(x) prod(psurvival[1:x]))
        
        
        
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
                                             n.chains=ncores, parallel=T, n.iter=niter/2,
                                             n.burnin=niter/2/2, n.thin=niter/2000/2,
                                             verbose=FALSE)
          if(idetection==1 & isurvival==1) print(Sys.time() - tstart)
        }
        
        
        
        # storing all the things
        est_detection_byriver[idetection, isurvival, , nsim] <-
          kenai_jags_out_hmi$q50$pdetection
        est_survival_byriver[idetection, isurvival, , nsim] <-
          kenai_jags_out_hmi$q50$psurvival
        sd_detection_byriver[idetection, isurvival, , nsim] <-
          kenai_jags_out_hmi$sd$pdetection
        sd_survival_byriver[idetection, isurvival, , nsim] <-
          kenai_jags_out_hmi$sd$psurvival
        
        # print(isim)
      }
    }
    print(idetection)
  }
  print(Sys.time() - t_overall_start)
  if(save_results) {
    save(pdetections, psurvivals, overall_survivals,
         est_detection_byriver, est_survival_byriver, 
         sd_detection_byriver, sd_survival_byriver,
         est_detection_bynearshore, est_survival_bynearshore, 
         sd_detection_bynearshore, sd_survival_bynearshore,
         pdetection_river_expected,
         pdetection_river_range,
         pdetection_nearshore_expected,
         pdetection_nearshore_range,
         
         psurvival_river_expected,
         psurvival_river_range,
         psurvival_nearshore_expected,
         psurvival_nearshore_range,
         file="OP_2025/data/Kenai_telem_rp_fixedprobs_simresults.Rdata")
  }
}
