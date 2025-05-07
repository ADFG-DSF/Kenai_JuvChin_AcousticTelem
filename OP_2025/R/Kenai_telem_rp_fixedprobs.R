# The purpose of this script is to further estimate the relative precision
# of the estimation of detection and survival probabilities for the Kenai River
# juvenile salmon telemetry project.

# Instead of simulating detection and survival probabilities (with uncertainty),
# these are treated as fixed, with two grids of possible values.
# * First, INRIVER detection and survival probabilities are treated as fixed and 
#   set to the _expected values in inputs, and NEARSHORE probabilities are allowed
#   to take on values according to a pre-defined grid, with outputs stored with
#   the suffix _bynearshore.
# * Next, NEARSHORE detection and survival probabilities are treated as fixed and 
#   set to the _expected values in inputs, and INRIVER probabilities are allowed
#   to take on values according to a pre-defined grid, with outputs stored with
#   the suffix _byriver.

# Following the results from Kenai_telem_rp.R, only the Hidden Markov (imputed)
# model is used for estimation, as this was substantially fastest.



# simulation controls
nfish <- 500
nstations_sep <- c(3,3) # number of NEARSHORE stations, then INRIVER stations
nstations <- sum(nstations_sep)


# defining probabilities
# _expected is held constant for river, while _range is allowed to vary for nearshore,
# and vice versa
pdetection_river_expected <- 0.95
pdetection_river_range <- seq(from=0.85, to=1, by=0.05)
pdetection_nearshore_expected <- 0.8
pdetection_nearshore_range <- seq(from=0.6, to=0.95, by=0.05)

psurvival_river_expected <- 0.95
psurvival_river_range <- seq(from=0.85, to=1, by=0.05)
psurvival_nearshore_expected <- 0.9
psurvival_nearshore_range <- seq(from=0.6, to=0.95, by=0.05)

# printing the total number of possibilities to try (to gauge runtime)
length(pdetection_nearshore_range)*length(psurvival_nearshore_range) +
  length(pdetection_river_range)*length(psurvival_river_range)
  



nsim <- 300   # 20 sim x 2k niter in 1 hr on desktop
save_results <- TRUE
run_model <- TRUE



# JAGS controls
niter <- 2000   # should make this bigger when we do it for real
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


if(run_model) {
  t_overall_start <- Sys.time()
  print(t_overall_start)
  
  # initializing everything
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
    est_detection_byriver <- 
    est_survival_byriver <- 
    est_overall_survival_byriver <- 
    sd_detection_byriver <- 
    sd_survival_byriver <- 
    sd_overall_survival_byriver <- 
    array(dim = c(length(pdetection_river_range),
                  length(psurvival_river_range),
                  nstations,
                  nsim))
  
  pdetections_bynearshore <- psurvivals_bynearshore <- overall_survivals_bynearshore <-
    array(dim=c(length(pdetection_nearshore_range),
                length(psurvival_nearshore_range),
                nstations))
  for(idetection in 1:length(pdetection_nearshore_range)) {
    for(isurvival in 1:length(psurvival_nearshore_range)) {
      
      pdetection <- c(rep(pdetection_river_expected, nstations_sep[1]),
                      rep(pdetection_nearshore_range[idetection], nstations_sep[2]))
      psurvival <- c(rep(psurvival_river_expected, nstations_sep[1]),
                     rep(psurvival_nearshore_range[isurvival], nstations_sep[2]))
      
      pdetections_bynearshore[idetection, isurvival, ] <- pdetection
      psurvivals_bynearshore[idetection, isurvival, ] <- psurvival
      overall_survivals_bynearshore[idetection, isurvival, ] <- 
        sapply(seq_along(psurvival), \(x) prod(psurvival[1:x]))
      
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
        est_detection_bynearshore[idetection, isurvival, , isim] <-
          kenai_jags_out_hmi$q50$pdetection
        est_survival_bynearshore[idetection, isurvival, , isim] <-
          kenai_jags_out_hmi$q50$psurvival
        est_overall_survival_bynearshore[idetection, isurvival, , isim] <-
          kenai_jags_out_hmi$q50$overall_survival
        sd_detection_bynearshore[idetection, isurvival, , isim] <-
          kenai_jags_out_hmi$sd$pdetection
        sd_survival_bynearshore[idetection, isurvival, , isim] <-
          kenai_jags_out_hmi$sd$psurvival
        sd_overall_survival_bynearshore[idetection, isurvival, , isim] <-
          kenai_jags_out_hmi$sd$overall_survival
        
        # print(isim)
      }
    }
    print(idetection)
  }
  
  pdetections_byriver <- psurvivals_byriver <- overall_survivals_byriver <-
    array(dim=c(length(pdetection_river_range),
                length(psurvival_river_range),
                nstations))
  for(idetection in 1:length(pdetection_river_range)) {
    for(isurvival in 1:length(psurvival_river_range)) {
      
      pdetection <- c(rep(pdetection_river_range[idetection], nstations_sep[1]),
                      rep(pdetection_nearshore_expected, nstations_sep[2]))
      psurvival <- c(rep(psurvival_river_range[isurvival], nstations_sep[1]),
                     rep(psurvival_nearshore_expected, nstations_sep[2]))
      
      pdetections_byriver[idetection, isurvival, ] <- pdetection
      psurvivals_byriver[idetection, isurvival, ] <- psurvival
      overall_survivals_byriver[idetection, isurvival, ] <- 
        sapply(seq_along(psurvival), \(x) prod(psurvival[1:x]))
      
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
        est_detection_byriver[idetection, isurvival, , isim] <-
          kenai_jags_out_hmi$q50$pdetection
        est_survival_byriver[idetection, isurvival, , isim] <-
          kenai_jags_out_hmi$q50$psurvival
        est_overall_survival_byriver[idetection, isurvival, , isim] <-
          kenai_jags_out_hmi$q50$overall_survival
        sd_detection_byriver[idetection, isurvival, , isim] <-
          kenai_jags_out_hmi$sd$pdetection
        sd_survival_byriver[idetection, isurvival, , isim] <-
          kenai_jags_out_hmi$sd$psurvival
        sd_overall_survival_byriver[idetection, isurvival, , isim] <-
          kenai_jags_out_hmi$sd$overall_survival
        
        # print(isim)
      }
    }
    print(idetection)
  }
  print(Sys.time())
  print(Sys.time() - t_overall_start)
  if(save_results) {
    save(est_detection_byriver, est_survival_byriver, est_overall_survival_byriver, 
         sd_detection_byriver, sd_survival_byriver, sd_overall_survival_byriver,
         est_detection_bynearshore, est_survival_bynearshore, est_overall_survival_bynearshore, 
         sd_detection_bynearshore, sd_survival_bynearshore, sd_overall_survival_bynearshore,
         pdetection_river_expected,
         pdetection_river_range,
         pdetection_nearshore_expected,
         pdetection_nearshore_range,
         
         psurvival_river_expected,
         psurvival_river_range,
         psurvival_nearshore_expected,
         psurvival_nearshore_range,
         
         pdetections_bynearshore,
         psurvivals_bynearshore,
         overall_survivals_bynearshore,
         pdetections_byriver,
         psurvivals_byriver,
         overall_survivals_byriver,
         
         nsim, niter, ncores, nfish,
         file="OP_2025/data/Kenai_telem_rp_fixedprobs_simresults.Rdata")
  }
} else {
  load(file="OP_2025/data/Kenai_telem_rp_fixedprobs_simresults.Rdata")
}

# investigate bias
mn_est_detection_byriver <- apply(est_detection_byriver, 1:3, mean)
mn_est_survival_byriver <- apply(est_survival_byriver, 1:3, mean)
mn_est_overall_survival_byriver <- apply(est_overall_survival_byriver, 1:3, mean)
mn_sd_detection_byriver <- apply(sd_detection_byriver, 1:3, mean)
mn_sd_survival_byriver <- apply(sd_survival_byriver, 1:3, mean)
mn_sd_overall_survival_byriver <- apply(sd_overall_survival_byriver, 1:3, mean)
mn_est_detection_bynearshore <- apply(est_detection_bynearshore, 1:3, mean)
mn_est_survival_bynearshore <- apply(est_survival_bynearshore, 1:3, mean)
mn_est_overall_survival_bynearshore <- apply(est_overall_survival_bynearshore, 1:3, mean)
mn_sd_detection_bynearshore <- apply(sd_detection_bynearshore, 1:3, mean)
mn_sd_survival_bynearshore <- apply(sd_survival_bynearshore, 1:3, mean)
mn_sd_overall_survival_bynearshore <- apply(sd_overall_survival_bynearshore, 1:3, mean)

par(mfrow=c(2,3))
for(i in 1:nstations) {
  contour(mn_est_detection_bynearshore[,,i] - pdetections_bynearshore[,,i], 
          x=pdetection_nearshore_range, y=psurvival_nearshore_range,
          xlab="nearshore detection", ylab="nearshore survival",
          main=paste("Bias detection",i))
  print(range(mn_est_detection_bynearshore[,,i] - pdetections_bynearshore[,,i]))
}
for(i in 1:nstations) {
  contour(mn_est_survival_bynearshore[,,i] - psurvivals_bynearshore[,,i], 
          x=pdetection_nearshore_range, y=psurvival_nearshore_range,
          xlab="nearshore detection", ylab="nearshore survival",
          main=paste("Bias survival",i))
  print(range(mn_est_survival_bynearshore[,,i] - psurvivals_bynearshore[,,i]))
}
for(i in 1:nstations) {
  contour(mn_est_overall_survival_bynearshore[,,i] - overall_survivals_bynearshore[,,i], 
          x=pdetection_nearshore_range, y=psurvival_nearshore_range,
          xlab="nearshore detection", ylab="nearshore survival",
          main=paste("Bias overall survival",i))
  print(range(mn_est_overall_survival_bynearshore[,,i] - overall_survivals_bynearshore[,,i]))
}

for(i in 1:nstations) {
  contour(mn_est_detection_byriver[,,i] - pdetections_byriver[,,i], 
          x=pdetection_river_range, y=psurvival_river_range,
          xlab="river detection", ylab="river survival",
          main=paste("Bias detection",i))
  print(range(mn_est_detection_byriver[,,i] - pdetections_byriver[,,i]))
}
for(i in 1:nstations) {
  contour(mn_est_survival_byriver[,,i] - psurvivals_byriver[,,i], 
          x=pdetection_river_range, y=psurvival_river_range,
          xlab="river detection", ylab="river survival",
          main=paste("Bias survival",i))
  print(range(mn_est_survival_byriver[,,i] - psurvivals_byriver[,,i]))
}
for(i in 1:nstations) {
  contour(mn_est_overall_survival_byriver[,,i] - overall_survivals_byriver[,,i], 
          x=pdetection_river_range, y=psurvival_river_range,
          xlab="river detection", ylab="river survival",
          main=paste("Bias overall survival",i))
  print(range(mn_est_overall_survival_byriver[,,i] - overall_survivals_byriver[,,i]))
}



# calculate rp
library(dsftools)
rp_detection_bynearshore <- rp_survival_bynearshore <- rp_overall_survival_bynearshore <-
  NA*pdetections_bynearshore
rp_detection_byriver <- rp_survival_byriver <- rp_overall_survival_byriver <-
  NA*pdetections_byriver

confidence <- 0.9  # should prob change this to 0.95 with more reps
relative <- FALSE  # setting this to FALSE results in absolute accuracy


for(idetection in 1:length(pdetection_nearshore_range)) {
  for(isurvival in 1:length(psurvival_nearshore_range)) { 
    for(istation in 1:nstations) {
      rp_detection_bynearshore[idetection, isurvival, istation] <- 
        rp(sim_vec=est_detection_bynearshore[idetection, isurvival, istation, ],
           true_val=pdetections_bynearshore[idetection, isurvival, istation],
           confidence=confidence, relative=relative)
      rp_survival_bynearshore[idetection, isurvival, istation] <- 
        rp(sim_vec=est_survival_bynearshore[idetection, isurvival, istation, ],
           true_val=psurvivals_bynearshore[idetection, isurvival, istation],
           confidence=confidence, relative=relative)
      rp_overall_survival_bynearshore[idetection, isurvival, istation] <- 
        rp(sim_vec=est_overall_survival_bynearshore[idetection, isurvival, istation, ],
           true_val=overall_survivals_bynearshore[idetection, isurvival, istation],
           confidence=confidence, relative=relative)
    }
  }
}
for(idetection in 1:length(pdetection_river_range)) {
  for(isurvival in 1:length(psurvival_river_range)) { 
    for(istation in 1:nstations) {
      rp_detection_byriver[idetection, isurvival, istation] <- 
        rp(sim_vec=est_detection_byriver[idetection, isurvival, istation, ],
           true_val=pdetections_byriver[idetection, isurvival, istation],
           confidence=confidence, relative=relative)
      rp_survival_byriver[idetection, isurvival, istation] <- 
        rp(sim_vec=est_survival_byriver[idetection, isurvival, istation, ],
           true_val=psurvivals_byriver[idetection, isurvival, istation],
           confidence=confidence, relative=relative)
      rp_overall_survival_byriver[idetection, isurvival, istation] <- 
        rp(sim_vec=est_overall_survival_byriver[idetection, isurvival, istation, ],
           true_val=overall_survivals_byriver[idetection, isurvival, istation],
           confidence=confidence, relative=relative)
    }
  }
}

par(mfrow=c(2,3))
for(i in 1:nstations) {
  contour(rp_detection_bynearshore[,,i], 
          x=pdetection_nearshore_range, y=psurvival_nearshore_range,
          xlab="nearshore detection", ylab="nearshore survival",
          main=paste("RP detection",i))
  print(range(rp_detection_bynearshore[,,i]))
}
for(i in 1:nstations) {
  contour(rp_survival_bynearshore[,,i], 
          x=pdetection_nearshore_range, y=psurvival_nearshore_range,
          xlab="nearshore detection", ylab="nearshore survival",
          main=paste("RP survival",i))
  print(range(rp_survival_bynearshore[,,i]))
}
for(i in 1:nstations) {
  contour(rp_overall_survival_bynearshore[,,i], 
          x=pdetection_nearshore_range, y=psurvival_nearshore_range,
          xlab="nearshore detection", ylab="nearshore survival",
          main=paste("RP overall survival",i))
  print(range(rp_overall_survival_bynearshore[,,i]))
}

par(mfrow=c(2,3))
for(i in 1:nstations) {
  contour(rp_detection_byriver[,,i], 
          x=pdetection_river_range, y=psurvival_river_range,
          xlab="river detection", ylab="river survival",
          main=paste("RP detection",i))
  print(range(rp_detection_byriver[,,i]))
}
for(i in 1:nstations) {
  contour(rp_survival_byriver[,,i], 
          x=pdetection_river_range, y=psurvival_river_range,
          xlab="river detection", ylab="river survival",
          main=paste("RP survival",i))
  print(range(rp_survival_byriver[,,i]))
}
for(i in 1:nstations) {
  contour(rp_overall_survival_byriver[,,i], 
          x=pdetection_river_range, y=psurvival_river_range,
          xlab="river detection", ylab="river survival",
          main=paste("RP overall survival",i))
  print(range(rp_overall_survival_byriver[,,i]))
}
