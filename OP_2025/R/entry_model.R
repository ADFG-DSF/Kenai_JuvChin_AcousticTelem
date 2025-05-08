# The purpose of this script is to further estimate the relative precision
# of the estimation of detection and survival probabilities for the Kenai River
# juvenile salmon telemetry project.

# This script reflects the most up-to-date version of the model, which allows
# late entry of individuals: that is, entry below the first array.  Since 
# handling mortality is no longer conflated with survival probability at the 
# first station, it is included as a parameter in simulation and in estimation.

# For the purpose of estimating relative precision, a large sequence of full 
# datasets are simulated according to the simulation controls below and the 
# Bayesian model is run for each and all parameter estimates are saved for the 
# purpose of comparing to assumed true values.  Note that uncertainty can be
# incorporated in probability parameter inputs by simulating these according
# to Beta distributions.

# Following the results from Kenai_telem_rp.R, only the Hidden Markov (imputed)
# model is used for estimation, as this was substantially fastest.



# metasimulation controls
nsim <- 10000   # 10k sim at 4k iter in 12 hr on desktop

run_model <- TRUE     # if TRUE, run the full simulation
                      # if FALSE, load the last saved results for plotting

plot_indiv_run <- FALSE  # if TRUE, produce model diagnostic plots
                         # if FALSE, suppress these

save_results <- TRUE  # if TRUE, automatically save simulation results
                      # if FALSE, suppress this

probabilities_fixed <- FALSE    # if TRUE, keep all probabilities to fixed values 
                                # if FALSE, allow all probabilities to vary
nbeta_sim <- 10 # effective sample size for beta distributions (bigger = more precise)



# simulation controls
nfish_sep <- c(200, 100, 100, 100)   # number of fish above each station, in sequence

phandling <- 0.9  # 1 minus capture & handling mortality

psurvival <- c(rep(0.95, 6), 0.9)   # 1 minus natural mortality at each station
pdetection <- c(rep(0.95, 6), 0.7)  # probability of detection at each station



# JAGS controls
niter <- 10*1000   # should make this bigger when we do it for real
ncores <- 6
# ncores <- min(10, parallel::detectCores()-1)





library(jagsUI)
library(jagshelper)



# derived quantities
nfish <- sum(nfish_sep)
nstations <- length(psurvival)
overall_survival <- sapply(seq_along(psurvival), \(x) prod(psurvival[1:x]))

if(!probabilities_fixed) {
  psurvival_input <- psurvival
  pdetection_input <- pdetection
  phandling_input <- phandling
}


# initializing storage of simulation params in case these are simulated
pdetections <- psurvivals <- overall_survivals <- matrix(nrow=nsim, ncol=nstations)
phandlings <- rep(NA, nsim)

# initializing storage of estimates
est_pdetections <- est_psurvivals <- est_overall_survivals <- matrix(nrow=nsim, ncol=nstations)
est_phandlings <- rep(NA, nsim)
sd_pdetections <- sd_psurvivals <- sd_overall_survivals <- matrix(nrow=nsim, ncol=nstations)
sd_phandlings <- rep(NA, nsim)

if(run_model) {
  t_overall_start <- Sys.time()
  for(isim in 1:nsim) {
    
    if(!probabilities_fixed) {
      psurvival <- rbeta(nstations, 
                         psurvival_input*nbeta_sim, (1-psurvival_input)*nbeta_sim)
      pdetection <- rbeta(nstations, 
                          pdetection_input*nbeta_sim, (1-pdetection_input)*nbeta_sim)
      phandling <- rbeta(1, 
                         phandling_input*nbeta_sim, (1-phandling_input)*nbeta_sim)
    }
    
    
    X <- Y <- matrix(nrow=nfish, ncol=nstations)
    # Z <- rep(NA, nfish)
    # X is (unobserved) survival
    # Y is survival & detection
    # Z is survived tagging
    
    
    # a vector of which station fish were tagged above
    entry <- unlist(lapply(seq_along(nfish_sep), \(x) rep(x, nfish_sep[x])))
    
    Z <- rbinom(nfish, size=1, prob=phandling)
    for(i in 1:nfish) {
      X[i, entry[i]] <- rbinom(1, size=Z[i], prob=psurvival[entry[i]])
      Y[i, entry[i]] <- rbinom(1, size=X[i, entry[i]], prob=pdetection[entry[i]])
      for(j in (entry[i]+1):nstations) {
        X[i,j] <- rbinom(1, size=X[i,j-1], prob=psurvival[j])
        Y[i,j] <- rbinom(1, size=X[i,j], prob=pdetection[j])
      }
    }
    Ximputed <- matrix(nrow=nfish, ncol=nstations)
    Zimputed <- rep(NA, nfish)
    for(i in 1:nfish) {
      if(any(Y[i,]==1, na.rm=TRUE)) {
        Zimputed[i] <- 1
        Ximputed[i,entry[i]:max(which(Y[i,]==1))] <- 1
      }
    }
    
    
    
    # HIDDEN MARKOV MODEL   - this is actually the only one!
    
    kenai_jags_hm_entry <- tempfile()
    cat('model {
  for(i in 1:nfish) {
    X[i, entry[i]] ~ dbinom(phandling*psurvival[entry[i]], 1)
    Y[i, entry[i]] ~ dbinom(pdetection[entry[i]], X[i,entry[i]])
    for(j in (entry[i]+1):nstations) {
      X[i, j] ~ dbinom(psurvival[j], X[i, j-1])
      Y[i, j] ~ dbinom(pdetection[j], X[i, j])
    }
  }
  
  for(j in 1:nstations) {
    psurvival[j] ~ dbeta(0.5, 0.5) 
    pdetection[j] ~ dbeta(0.5, 0.5) 
    overall_survival[j] <- prod(psurvival[1:j])
  }
  
  phandling ~ dbeta(0.5, 0.5)
  
}', file=kenai_jags_hm_entry)
    
    kenai_data_hm_entry <- list(nfish=nfish,
                                nstations=nstations,
                                X=Ximputed,
                                Y=Y,
                                entry=entry)
    
    {
      tstart <- Sys.time()
      # print(tstart)
      kenai_jags_out_hm_entry <- jagsUI::jags(model.file=kenai_jags_hm_entry, data=kenai_data_hm_entry,
                                              parameters.to.save=c("psurvival", "pdetection",
                                                                   "overall_survival",
                                                                   "phandling"),
                                              n.chains=ncores, parallel=T, n.iter=niter,
                                              n.burnin=niter/2, n.thin=niter/2000,
                                              verbose=FALSE)
      if(isim==1) print(Sys.time() - tstart)
    }
    
    print(isim)
    
    # storing sim params
    pdetections[isim, ] <- pdetection
    psurvivals[isim, ] <- psurvival
    overall_survivals[isim, ] <- overall_survival
    phandlings[isim] <- phandling
    
    # storing estimates
    est_pdetections[isim, ] <- kenai_jags_out_hm_entry$q50$pdetection 
    est_psurvivals[isim, ] <- kenai_jags_out_hm_entry$q50$psurvival
    est_overall_survivals[isim, ] <- kenai_jags_out_hm_entry$q50$overall_survival
    est_phandlings[isim] <- kenai_jags_out_hm_entry$q50$phandling
    sd_pdetections[isim, ] <- kenai_jags_out_hm_entry$sd$pdetection
    sd_psurvivals[isim, ] <- kenai_jags_out_hm_entry$sd$psurvival
    sd_overall_survivals[isim, ] <- kenai_jags_out_hm_entry$sd$overall_survival
    sd_phandlings[isim] <- kenai_jags_out_hm_entry$sd$phandling
    
    if(plot_indiv_run) {
      plotRhats(kenai_jags_out_hm_entry)
      tracedens_jags(kenai_jags_out_hm_entry, p="pdetection", parmfrow=c(3,3))
      tracedens_jags(kenai_jags_out_hm_entry, p="psurvival", parmfrow=c(3,3))
      tracedens_jags(kenai_jags_out_hm_entry, p="overall_survival", parmfrow=c(3,3))
      tracedens_jags(kenai_jags_out_hm_entry, p="phandling", parmfrow=c(3,3))
      
      parmar <- par("mar")
      par(mar=c(8,8,4,2))
      plotcor_jags(kenai_jags_out_hm_entry)
      par(mar=parmar)
      
      par(mfrow=c(2,2))
      caterpillar(kenai_jags_out_hm_entry, p="pdetection")
      points(pdetection)
      caterpillar(kenai_jags_out_hm_entry, p="psurvival")
      points(psurvival)
      caterpillar(kenai_jags_out_hm_entry, p="overall_survival")
      points(overall_survival)
      caterpillar(kenai_jags_out_hm_entry, p="phandling")
      points(phandling)
    }
  
  if(save_results & (isim %% 100 == 0)) {
    print("saving...")
    save(pdetections,
         psurvivals,
         overall_survivals,
         phandlings,
         
         est_pdetections,
         est_psurvivals,
         est_overall_survivals,
         est_phandlings,
         
         sd_pdetections,
         sd_psurvivals,
         sd_overall_survivals,
         sd_phandlings,
         
         nfish_sep, nsim, niter, ncores,
         file="OP_2025/data/Kenai_telem_hm_entry_simresults.Rdata")
  }
  }
  print(Sys.time() - t_overall_start)
} else {
  load(file="OP_2025/data/Kenai_telem_hm_entry_simresults.Rdata")
}

#### this section is only relevant if all probabilities are fixed!!
par(mfrow=c(2,2))
caterpillar(est_pdetections, main="Est detection probability")
points(pdetection)
caterpillar(est_psurvivals, main="Est survival probability")
points(psurvival)
caterpillar(est_overall_survivals, main="Est cumulative survival probability")
points(overall_survival)
caterpillar(est_phandlings, main="Est handling survival")
points(phandling)

par(mfrow=c(2,2))
caterpillar(sd_pdetections, main="SD detection probability")
caterpillar(sd_psurvivals, main="SD survival probability")
caterpillar(sd_overall_survivals, main="SD cumulative survival probability")
caterpillar(sd_phandlings, main="SD handling survival")


## Estimation bias
par(mfrow=c(2,2))
caterpillar(est_pdetections - pdetections, main="Est Bias - detection probability")
abline(h=0, lty=1)
caterpillar(est_psurvivals - psurvivals, main="Est Bias - survival probability")
abline(h=0, lty=1)
caterpillar(est_overall_survivals - overall_survivals, main="Est Bias - cumulative survival probability")
abline(h=0, lty=1)
caterpillar(est_phandlings - phandlings, main="Est Bias - handling survival")
abline(h=0, lty=1)


## RP in terms of absolute accuracy
library(dsftools)
confidence <- 0.95 # change this to 95% when we have more sims
relative <- FALSE # FALSE gives accuracy in absolute terms
rp_pdetections <- rp_psurvivals <- rp_overall_survivals <- rep(NA, nstations)
for(istation in 1:nstations) {
  rp_pdetections[istation] <- rp(sim_vec = est_pdetections[, istation],
                                 true_val = pdetections[, istation],
                                 confidence = confidence, relative = relative)
  rp_psurvivals[istation] <- rp(sim_vec = est_psurvivals[, istation],
                                 true_val = psurvivals[, istation],
                                 confidence = confidence, relative = relative)
  rp_overall_survivals[istation] <- rp(sim_vec = est_overall_survivals[, istation],
                                 true_val = overall_survivals[, istation],
                                 confidence = confidence, relative = relative)
}
rp_phandlings <- rp(sim_vec = est_phandlings,
                               true_val = phandlings,
                               confidence = confidence, relative = relative)
plot(rp_pdetections, main="Abs accuracy - detection probability")
plot(rp_psurvivals, main="Abs accuracy - survival probability")
plot(rp_overall_survivals, main="Abs accuracy - cumulative survival probability")
plot(rp_phandlings, main="Abs accuracy - handling survival")
