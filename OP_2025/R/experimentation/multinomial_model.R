# The purpose of this script was to create a model that generalizes a 
# multinomial likelihood, in which the probabilities represent the probability
# of EACH POSSIBLE CAPTURE HISTORY, as a function of the vectors of survival 
# probability and detection probability.  It is worth noting that each element
# of survival probability is conditional on survival at the previous time period,
# and each element of detection probability is conditional on survival at that
# time period.

# Data is first simulated according to the numbers of fish and stations, in which
# X represents the matrix of survival as an unobserved state, structured as a 
# matrix of zeroes and ones with a row for each fish and a column for each station.
# Y represents the detection matrix (observed), with the same structure.

# Then, all possible capture histories are expanded (2^n_stations possible), and
# simulated data are matched to possible histories.

# An array is then defined, which consists of all possible states that could give
# rise to each possible capture history.  This is used within the multinomial
# model to generalize the likelihood.

# Finally, a Bayesian model is defined and run with the simulated data.
# For each possible capture history:
#   for each state that would give rise to it:
#      the product is taken of all associated conditional probabilities 
#      (survival & detection, and their complements)
#   and these are summed.



###### data simulation #####

nfish <- 500
nstations_sep <- c(3,3)

pdetection <- c(rbeta(nstations_sep[1], 18, 1), 
                rbeta(nstations_sep[2], 7, 3))

psurvival <- c(rbeta(nstations_sep[1], 18, 1),
               rbeta(nstations_sep[2], 15, 2))

# phandling <- 0.9

overall_survival <- sapply(seq_along(psurvival), \(x) prod(psurvival[1:x]))
nstations <- sum(nstations_sep)

X <- Y <- matrix(nrow=nfish, ncol=nstations)
# X is (unobserved) survival
# Y is survival & detection

X[,1] <- rbinom(nfish, 1, psurvival[1])
# X[,1] <- rbinom(nfish, rbinom(nfish, 1, phandling), psurvival[1])
for(j in 2:nstations) {
  X[,j] <- rbinom(nfish, X[,j-1], psurvival[j])
}
for(j in 1:nstations) Y[,j] <- rbinom(nfish, X[,j], pdetection[j])
colMeans(X)
colMeans(Y)
overall_survival



#### expansion of all possible capture histories ####

library(magrittr)
possible_histories <- replicate(nstations, 0:1, simplify=FALSE) %>%
  expand.grid %>% as.matrix %>%
  apply(1, paste0, collapse="")

Yhistories <- apply(Y, 1, paste0, collapse="") %>%
  sapply(\(x) which(x==possible_histories)) %>% unname


empirical_props <- table(factor(Yhistories, levels=seq(length(possible_histories)))) %>%
  as.numeric/nfish


p <- pdetection
phi <- psurvival      

## validation of the generalized method, in the case of three stations

# Var1 Var2 Var3
# [1,]    0    0    0
# [2,]    1    0    0
# [3,]    0    1    0
# [4,]    1    1    0
# [5,]    0    0    1
# [6,]    1    0    1
# [7,]    0    1    1
# [8,]    1    1    1

# pi_theo <- c(
#   phi[1]*(1-p[1]) * phi[2]*(1-p[2]) * (1-(p[3]*phi[3])) +
#     phi[1]*(1-p[1]) * (1-phi[2]) +
#     (1-phi[1]),
#   phi[1]*p[1]*(phi[2]*(1-p[2]) * (1-(p[3]*phi[3])) +
#                  (1-phi[2])),
#   phi[1]*(1-p[1]) * phi[2]*p[2] * (1-(phi[3]*p[3])),
#   phi[1]*p[1] * phi[2]*p[2] * (1-(phi[3]*p[3])),
#   phi[1]*(1-p[1]) * phi[2]*(1-p[2]) * phi[3]*p[3],
#   phi[1]*p[1] * phi[2]*(1-p[2]) * phi[3]*p[3],
#   phi[1]*(1-p[1]) * phi[2]*p[2] * phi[3]*p[3],
#   phi[1]*p[1] * phi[2]*p[2] * phi[3]*p[3]
#   )
# sum(pi_theo)
# 
# plot(pi_theo, empirical_props)
# abline(0,1)




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



##### construction of a Bayesian model #####


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
    psurvival[jstation] ~ dbeta(9,1)
    pdetection[jstation] ~ dbeta(9,1)
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

nbyname(kenai_jags_out)
plotRhats(kenai_jags_out)
traceworstRhat(kenai_jags_out, parmfrow = c(3,3))

par(mfrow=c(1,2))
caterpillar(kenai_jags_out, p="psurvival", ylim=0:1)
points(psurvival)
caterpillar(kenai_jags_out, p="pdetection", ylim=0:1)
points(pdetection)
caterpillar(kenai_jags_out, p="pi")
points(empirical_props)

caterpillar(kenai_jags_out, p="overall_survival")
points(overall_survival)

# plotcor_jags(kenai_jags_out, p=c("psurvival", "pdetection"), maxn=100)