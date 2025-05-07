# simulation controls
nfish <- 500
nstations_sep <- c(3,3) #c(3,3)
nstations <- sum(nstations_sep)



pdetection_beta_river <- c(19, 1)
pdetection_beta_nearshore <- c(7, 3)

psurvival_beta_river <- c(18, 1)
psurvival_beta_nearshore <- c(15, 2)


nsim <- 500   # 50 sim in 1.5 hrs on laptop
save_results <- TRUE
run_model <- FALSE



# JAGS controls
niter <- 8000   # make this divisible by 4k
ncores <- 6
# ncores <- min(10, parallel::detectCores()-1)




library(jagsUI)
library(jagshelper)
library(magrittr)



if(run_model) {

xmeta <- seq(from=0, to=1, by=0.01)
ymeta <- data.frame(dbeta(xmeta, 
                          pdetection_beta_river[1], 
                          pdetection_beta_river[2]),
                    dbeta(xmeta, 
                          pdetection_beta_nearshore[1], 
                          pdetection_beta_nearshore[2]),
                    dbeta(xmeta, 
                          psurvival_beta_river[1], 
                          psurvival_beta_river[2]),
                    dbeta(xmeta, 
                          psurvival_beta_nearshore[1], 
                          psurvival_beta_nearshore[2]))
distmns <- sapply(list(pdetection_beta_river,
                       pdetection_beta_nearshore,
                       psurvival_beta_river,
                       psurvival_beta_nearshore),
                  \(x) x[1]/sum(x))
# # par(mfrow=c(2,1))
# plot(NA, xlim=0:1, ylim=range(0,ymeta[,1:2]))
# for(j in 1:2) lines(xmeta, ymeta[,j], col=c(4,2,4,2)[j], lty=c(2,2,1,1)[j])
# plot(NA, xlim=0:1, ylim=range(0,ymeta[,3:4]))
# for(j in 3:4) lines(xmeta, ymeta[,j], col=c(4,2,4,2)[j], lty=c(2,2,1,1)[j])

par(mfrow=c(1,1))
plot(NA, xlim=0:1, ylim=range(0,ymeta))
for(j in 1:4) {
  lines(xmeta, ymeta[,j], col=c(4,2,4,2)[j], lty=c(2,2,1,1)[j], lwd=2)
  segments(x0=distmns[j],   
           y0=0,
           y1=ymeta[,j][which.max(xmeta >= distmns[j])],
           col=c(4,2,4,2)[j], lty=c(2,2,1,1)[j],
           lwd=2)
}
legend("topleft", 
       legend=paste(c("river","nearshore"), rep(c("detection","survival"), each=2)),
       col=c(4,2,4,2), lty=c(2,2,1,1),lwd=2)





# define candidate models

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

{
t_overall_start <- Sys.time()
for(isim in 1:nsim) {
  
  # first, simulating detection & survival probabilities
  # (separately for each simulated replicate)
  
  pdetection <- c(rbeta(nstations_sep[1], 
                        pdetection_beta_river[1], 
                        pdetection_beta_river[2]), 
                  rbeta(nstations_sep[2], 
                        pdetection_beta_nearshore[1], 
                        pdetection_beta_nearshore[2]))
  
  psurvival <- c(rbeta(nstations_sep[1], 
                       psurvival_beta_river[1], 
                       psurvival_beta_river[2]), 
                 rbeta(nstations_sep[2], 
                       psurvival_beta_nearshore[1], 
                       psurvival_beta_nearshore[2]))
  
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
                                       n.chains=ncores, parallel=T, n.iter=niter/2,
                                       n.burnin=niter/2/2, n.thin=niter/2000/2,
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
  
  print(isim)
}
print(Sys.time() - t_overall_start)
if(save_results) {
  save(pdetections, psurvivals, overall_survivals,
       est_detection_hm, est_detection_hmi, est_detection_mn,
       est_survival_hm, est_survival_hmi, est_survival_mn,
       est_overall_survival_hm, est_overall_survival_hmi, est_overall_survival_mn,
       sd_detection_hm, sd_detection_hmi, sd_detection_mn,
       sd_survival_hm, sd_survival_hmi, sd_survival_mn,
       sd_overall_survival_hm, sd_overall_survival_hmi, sd_overall_survival_mn,
       niter, ncores, nsim,
       file="OP_2025/data/Kenai_telem_rp_simresults.Rdata")
}
}
} else {
  load(file="OP_2025/data/Kenai_telem_rp_simresults.Rdata")
}


# should look better at convergence between all (last run)
tracedens_jags(kenai_jags_out_hm, 
               p=c("pdetection","psurvival","overall_survival"),
               parmfrow=c(3,2))
tracedens_jags(kenai_jags_out_hmi, 
               p=c("pdetection","psurvival","overall_survival"),
               parmfrow=c(3,2))
tracedens_jags(kenai_jags_out_hm, 
               p=c("pdetection","psurvival","overall_survival"),
               parmfrow=c(3,2))

plotRhats(kenai_jags_out_hm)
plotRhats(kenai_jags_out_hmi)
plotRhats(kenai_jags_out_mn)

# is there bias?
par(mfrow=c(1,3))
caterpillar(est_detection_hm - pdetections, main="detection bias - hm")
abline(h=0)
caterpillar(est_detection_hmi - pdetections, main="detection bias - hmi")
abline(h=0)
caterpillar(est_detection_mn - pdetections, main="detection bias - mn")
abline(h=0)

caterpillar(est_survival_hm - psurvivals, main="survival bias - hm")
abline(h=0)
caterpillar(est_survival_hmi - psurvivals, main="survival bias - hmi")
abline(h=0)
caterpillar(est_survival_mn - psurvivals, main="survival bias - mn")
abline(h=0)

caterpillar(est_overall_survival_hm - overall_survivals, main="overall survival bias - hm")
abline(h=0)
caterpillar(est_overall_survival_hmi - overall_survivals, main="overall survival bias - hmi")
abline(h=0)
caterpillar(est_overall_survival_mn - overall_survivals, main="overall survival bias - mn")
abline(h=0)

# compare sd
caterpillar(sd_detection_hm, main="detection sd - hm")
caterpillar(sd_detection_hmi, main="detection sd - hmi")
caterpillar(sd_detection_mn, main="detection sd - mn")

caterpillar(sd_survival_hm, main="survival sd - hm")
caterpillar(sd_survival_hmi, main="survival sd - hmi")
caterpillar(sd_survival_mn, main="survival sd - mn")

caterpillar(sd_overall_survival_hm, main="overall survival sd - hm")
caterpillar(sd_overall_survival_hmi, main="overall survival sd - hmi")
caterpillar(sd_overall_survival_mn, main="overall survival sd - mn")

# can we do something like rmse? I would think so
# rmse <- function(x,y) sqrt(mean((x-y)^2, na.rm=TRUE))
rmse_detection_hm <- (est_detection_hm - pdetections)^2 %>% colMeans %>% sqrt 
rmse_detection_hmi <- (est_detection_hmi - pdetections)^2 %>% colMeans %>% sqrt 
rmse_detection_mn <- (est_detection_mn - pdetections)^2 %>% colMeans %>% sqrt
rmse_survival_hm <- (est_survival_hm - psurvivals)^2 %>% colMeans %>% sqrt 
rmse_survival_hmi <- (est_survival_hmi - psurvivals)^2 %>% colMeans %>% sqrt 
rmse_survival_mn <- (est_survival_mn - psurvivals)^2 %>% colMeans %>% sqrt 
rmse_overall_survival_hm <- (est_overall_survival_hm - overall_survivals)^2 %>% colMeans %>% sqrt 
rmse_overall_survival_hmi <- (est_overall_survival_hmi - overall_survivals)^2 %>% colMeans %>% sqrt 
rmse_overall_survival_mn <- (est_overall_survival_mn - overall_survivals)^2 %>% colMeans %>% sqrt 
plotthemboth <- function(x,y,z,...) {
  plot(x, pch="+", col=2, ylim=range(0,x,y,z), ...=...)
  points(y, pch="+", col=4)
  points(z, pch="+", col=3)
}
plotthemboth(rmse_detection_hm, rmse_detection_hmi, rmse_detection_mn, 
             main="rmse detection")
plotthemboth(rmse_survival_hm, rmse_survival_hmi, rmse_survival_mn,
             main="rmse survival")
plotthemboth(rmse_overall_survival_hm, rmse_overall_survival_hmi, rmse_overall_survival_mn,
             main="rmse overall survival")


# compare rp - this will be the kicker
relative <- FALSE  # FALSE will calculate absolute accuracy

rp_survival_hm <- rp_detection_hm <- rp_overall_survival_hm <- rep(NA, ncol(psurvivals))
for(j in seq(nstations)) {
  rp_survival_hm[j] <- dsftools::rp(est_survival_hm[,j], psurvivals[,j], 
                                    confidence=0.95, relative=relative)
  rp_overall_survival_hm[j] <- dsftools::rp(est_overall_survival_hm[,j], overall_survivals[,j], 
                                            confidence=0.95, relative=relative)
  rp_detection_hm[j] <- dsftools::rp(est_detection_hm[,j], pdetections[,j], 
                                     confidence=0.95, relative=relative)
}

rp_survival_hmi <- rp_detection_hmi <- rp_overall_survival_hmi <- rep(NA, ncol(psurvivals))
for(j in seq(nstations)) {
  rp_survival_hmi[j] <- dsftools::rp(est_survival_hmi[,j], psurvivals[,j], 
                                     confidence=0.95, relative=relative)
  rp_overall_survival_hmi[j] <- dsftools::rp(est_overall_survival_hmi[,j], overall_survivals[,j], 
                                             confidence=0.95, relative=relative)
  rp_detection_hmi[j] <- dsftools::rp(est_detection_hmi[,j], pdetections[,j], 
                                      confidence=0.95, relative=relative)
}

rp_survival_mn <- rp_detection_mn <- rp_overall_survival_mn <- rep(NA, ncol(psurvivals))
for(j in seq(nstations)) {
  rp_survival_mn[j] <- dsftools::rp(est_survival_mn[,j], psurvivals[,j], 
                                    confidence=0.95, relative=relative)
  rp_overall_survival_mn[j] <- dsftools::rp(est_overall_survival_mn[,j], overall_survivals[,j], 
                                            confidence=0.95, relative=relative)
  rp_detection_mn[j] <- dsftools::rp(est_detection_mn[,j], pdetections[,j], 
                                     confidence=0.95, relative=relative)
}
plotthemboth(rp_detection_hm, rp_detection_hmi, rp_detection_mn, 
             main="rp detection")
plotthemboth(rp_survival_hm, rp_survival_hmi, rp_survival_mn,
             main="rp survival")
plotthemboth(rp_overall_survival_hm, rp_overall_survival_hmi, rp_overall_survival_mn,
             main="rp overall survival")


# is there something like simulataneous rp??