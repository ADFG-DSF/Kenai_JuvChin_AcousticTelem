nfish_sep <- c(100, 100, 100, 100, 100) # 100 above the first 5 arrays
nfish <- sum(nfish_sep)

phandling <- 0.9  # 1 - capture & handling mortality

psurvival <- c(rep(0.95, 6), 0.9)
pdetection <- c(rep(0.95, 6), 0.7)

nstations <- length(psurvival)
overall_survival <- sapply(seq_along(psurvival), \(x) prod(psurvival[1:x]))




X <- Y <- matrix(nrow=nfish, ncol=nstations)
Z <- rep(NA, nfish)
# X is (unobserved) survival
# Y is survival & detection
# Z is survived tagging


# no idea if this will work, it should be 1 1 1 x100, 2 2 2 x100, etc
# a vector of which station fish were tagged above
tagged <- unlist(lapply(seq_along(nfish_sep), \(x) rep(x, nfish_sep[x])))

Z <- rbinom(nfish, size=1, prob=phandling)
for(i in 1:nfish) {
  X[i, tagged[i]] <- rbinom(1, size=Z[i], prob=psurvival[tagged[i]])
  Y[i, tagged[i]] <- rbinom(1, size=X[i, tagged[i]], prob=pdetection[tagged[i]])
  for(j in (tagged[i]+1):nstations) {
    X[i,j] <- rbinom(1, size=X[i,j-1], prob=psurvival[j])
    Y[i,j] <- rbinom(1, size=X[i,j], prob=pdetection[j])
  }
}
Ximputed <- matrix(nrow=nfish, ncol=nstations)
Zimputed <- rep(NA, nfish)
for(i in 1:nfish) {
  if(any(Y[i,]==1, na.rm=TRUE)) {
    Zimputed[i] <- 1
    Ximputed[i,tagged[i]:max(which(Y[i,]==1))] <- 1
  }
}


# 
# 
# X[,1] <- rbinom(nfish, 1, psurvival[1])
# # X[,1] <- rbinom(nfish, rbinom(nfish, 1, phandling), psurvival[1])
# for(j in 2:nstations) {
#   X[,j] <- rbinom(nfish, X[,j-1], psurvival[j])
# }
# 
# for(j in 1:nstations) Y[,j] <- rbinom(nfish, X[,j], pdetection[j])
# 
# Ximputed <- matrix(nrow=nfish, ncol=nstations)
# for(i in 1:nfish) {
#   if(any(Y[i,]==1)) {
#     Ximputed[i,1:max(which(Y[i,]==1))] <- 1
#   }
# }

colMeans(X)
colMeans(Y)
overall_survival