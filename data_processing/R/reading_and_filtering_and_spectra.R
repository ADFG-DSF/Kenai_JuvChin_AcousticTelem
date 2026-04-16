library(tidyverse)

### reading data

# maybe if we wrap this in a function, it will be easier to batch-process later

readJSATS <- function(filename) {
  # detecting the regions that correspond to data tables
  as_lines <- suppressWarnings(readLines(filename))
  starts <- which(substr(as_lines, 1, 8) == "Internal")
  ends <- which(substr(as_lines, 1, 8) == "File End")
  
  # reading tables and compiling as a list
  raw_data <- list()
  for(i in seq_along(starts)) {
    # raw_data[[i]] <- read.csv(filename, 
    #                           skip=starts[i]+1, 
    #                           nrows=ends[i]-starts[i]-2,
    #                           header=FALSE)
    # colnames(raw_data[[i]]) <- colnames(read.csv(filename, 
    #                                              skip=starts[i]-1, 
    #                                              nrows=0,
    #                                              header=TRUE))
    
    # let's do this tidyverse style (probably cleaner)
    raw_data[[i]] <- read_csv(filename,
                              skip = starts[i]-1,
                              n_max = ends[i] - starts[i] - 2)
    
    # converting DateTime from character 
    raw_data[[i]]$DateTime <- as.POSIXct(raw_data[[i]]$DateTime, 
                                         format="%m/%d/%Y %H:%M:%OS")
    
    ##### need to figure out how site name is handled??
    # site name seems to have changed between two data tables (only in preamble for first)
    # adding input filename (without file address)
    namesplit <- strsplit(filename, split="/")[[1]]
    raw_data[[i]]$filename <- namesplit[length(namesplit)]
  }
  
  # apply filtration function?  WRITE THIS SEPARATELY
  
  # write each data table to an external csv file
  
  # this is probably temporary for testing
  return(raw_data)
}

raw_data <- readJSATS(filename = "data_processing/test_data/BeaverCreek012_250527_145301.csv")
# write a separate function to filter, call that in the func above^^


## --------------- testing zone ----------------------

asdf <- raw_data[[2]]
summary(asdf)
head(asdf)

# asdf$DateTime <- as.POSIXct(asdf$DateTime, format="%m/%d/%Y %H:%M:%OS")


# unique_codes <- sort(unique(asdf$TagCode)) # too many  #52k
# 
# codes_table <- table(asdf$TagCode)

# unique_codes <- names(codes_table[codes_table > 100]) #18k with >2

# making a list of times for all unique tags (converted to seconds)
datetime_secs <- tapply(asdf$DateTime, asdf$TagCode, FUN=as.numeric)

# taking the subset of this for entries with >n records
datetime_secs_sub <- datetime_secs[sapply(datetime_secs, length) > 50]

# differences in times between successive observations
datetime_diffs <- lapply(datetime_secs_sub, \(x) diff(sort(x)))

# subset with diffs that are ~3 secs
datetime_diffs_3secs <- datetime_diffs[sapply(datetime_diffs, \(x) any(x > 2 & x < 4))]

par(mfrow=c(3,3))
for(i in 1:length(datetime_diffs_3secs)) plot(datetime_diffs_3secs[[i]], log="y", main=names(datetime_diffs_3secs)[i])

# # tables of possible periodic components?
# datetime_diffs_tbls


# I'm curious: how many 
npertag <- as.numeric(table(asdf$TagCode))
trialvals <- c(1:10, 
               seq(15, 25, by=5),
               30, 40,
               seq(50, 250, by=50),
               seq(300, 1000, by=100))
nperval <- sapply(trialvals, \(x) sum(npertag>=x))
plot(trialvals, nperval, log="xy")


## TRIAL FILTRATION METHOD:
## keep only entries for which another entry exists with 
## - the same tag code, and
## - is 2.5-3.5 seconds away


# if we keep this, it will be wrapped in a function and probably called
# within the data import function readJSATS()

# first sort by TagCode, then by DateTime
asdf$orig_row <- 1:nrow(asdf)
asdf$datetimenum <- as.numeric(asdf$DateTime)
asdf_sort <- asdf
asdf_sort$sortby <- paste(asdf_sort$TagCode, asdf_sort$DateTime)
asdf_sort <- asdf_sort[order(asdf_sort$sortby),]

# initialize logical vectors for conditions in which we want to keep the entry
keepthese1 <- keepthese2 <- rep(F, nrow(asdf_sort))

dt <- 0.5  # how close to 3 seconds apart do we require?

# detecting which rows to keep, then subsetting
keepthese1[1:(nrow(asdf_sort)-1)] <- 
  (asdf_sort$TagCode[1:(nrow(asdf_sort)-1)] == asdf_sort$TagCode[2:(nrow(asdf_sort))]) &
  (abs(asdf_sort$datetimenum[1:(nrow(asdf_sort)-1)] - asdf_sort$datetimenum[2:(nrow(asdf_sort))]) < (3+dt)) &
  (abs(asdf_sort$datetimenum[1:(nrow(asdf_sort)-1)] - asdf_sort$datetimenum[2:(nrow(asdf_sort))]) > (3-dt))
keepthese2[2:(nrow(asdf_sort))] <- 
  (asdf_sort$TagCode[1:(nrow(asdf_sort)-1)] == asdf_sort$TagCode[2:(nrow(asdf_sort))]) &
  (abs(asdf_sort$datetimenum[1:(nrow(asdf_sort)-1)] - asdf_sort$datetimenum[2:(nrow(asdf_sort))]) < (3+dt)) &
  (abs(asdf_sort$datetimenum[1:(nrow(asdf_sort)-1)] - asdf_sort$datetimenum[2:(nrow(asdf_sort))]) > (3-dt))
asdf_sort <- subset(asdf_sort, keepthese1 | keepthese2)



# exploring the output!
length(unique(asdf_sort$TagCode))  # how many unique tag codes
table(asdf_sort$TagCode)           # how many entries per tag code
nrow(asdf_sort)                    # how many total entries, the original was 226k


# plotting all records with respect to time.
# TagCode is on the y-axis, DateTime on the x-axis
# horizontal line segments link all records with a single TagCode
# impressions:
# - after June 9 or so, fish passed by the station relatively quickly
# - something else was going on in mid May
# - tag G7205A5CF was maybe tested early on, and deployed in a fish?
# - tag G72052B5C had over 23k records all by itself!
par(mfrow=c(1,1))
with(asdf_sort, plot(y=as.numeric(as.factor(TagCode)), x=DateTime,
                     col=as.numeric(as.factor(TagCode)),
                     ylab="", yaxt="n"))
axis(side=2, las=2, cex.axis=.5,
     at=1:length(unique(asdf_sort$TagCode)), 
     labels=levels(as.factor(asdf_sort$TagCode)))
segments(x0=tapply(asdf_sort$DateTime, asdf_sort$TagCode, min),
         x1=tapply(asdf_sort$DateTime, asdf_sort$TagCode, max),
         y0=1:length(unique(asdf_sort$TagCode)))

# plotting Tilt with respect to DateTime further suggests that late May was something else
with(asdf_sort, plot(DateTime, Tilt))

# which tag codes were identified as fish using this method?
realfish1 <- sort(unique(asdf_sort$TagCode))



## EXPERIMENTING WITH SPECTRAL METHODS
## This function investigates whether a periodic signal exists, given a sequence of trial frequencies
## inputs are:
## - x: a numeric vector, basically as.numeric(DateTime) for an indiv TagCode
## - lam_min and lam_max: min and max signal PERIODS to consider (in seconds)
## - n_lam: length of the sequence of trial periods to consider
##          large number gives better resolution, but is more time-consuming
## - log_scale: whether to construct the trial sequence on the log scale (the default) or natural
##              log_scale==TRUE will give a finer resolution at small values

spectra <- function(x, lam_min=1, lam_max=120, n_lam=1000, log_scale=TRUE) { #, lam_trial=exp(seq(from=0, to=log(10), length.out=1000))) {
  if(log_scale) {
    lam_trial <- exp(seq(from=log(lam_min), to=log(lam_max), length.out=n_lam))
  } else {
    lam_trial <- seq(from=lam_min, to=lam_max, length.out=n_lam)
  }
  strength <- pstrength <- rep(NA, length(lam_trial))
  
  # these might need to be function inputs later on
  nbin <- 100  # resolution for binning
  delta <- 2   # number of adjacent bins to include in summation
  alpha_sidak <- 1-(.99^(1/length(lam_trial))) # Sidak-corrected significance level
  
  for(i in seq_along(strength)) {
    folded <- x %% lam_trial[i]
    
    # # Using Gaussian kernel density
    # thedens <- density(folded, n=nbin)
    # whichmaxes <- (which.max(thedens$y) + (-delta):delta) %% nbin
    # whichmaxes[whichmaxes == 0] <- nbin
    # strength[i] <- sum(thedens$y[whichmaxes])/sum(thedens$y)
    # strengthvec <- thedens$y
    
    # Just binning at a regular interval
    thetab <- as.numeric(table(cut(folded, breaks=seq(from=min(folded), to=max(folded), length.out=nbin+1))))
    whichmaxes <- (which.max(thetab) + (-delta):delta) %% nbin
    whichmaxes[whichmaxes == 0] <- nbin
    strengthvec <- thetab
    
    # calculated "strength" of spectral peak
    strength[i] <- sum(strengthvec[whichmaxes])/sum(strengthvec)
    
    # approx corresponding p-value (NEEDS TO BE CORRECTED!!)
    pstrength[i] <- pbinom(q=sum(strengthvec[whichmaxes]), 
                           size=length(x),
                           prob=(1+(2*delta))/nbin,
                           lower.tail=FALSE)
  }
  
  # strength that would be "significant" at family-wise Sidak level
  strength_sidak <- qbinom(p=alpha_sidak,
                           size=length(x),
                           prob=(1+(2*delta))/nbin,
                           lower.tail=FALSE)/sum(strengthvec)
  
  return(list(lam_trial=lam_trial, 
              strength=strength, 
              pval=pstrength,
              alpha_sidak=alpha_sidak,
              strength_sidak=strength_sidak))
}

# this will be used to explore additional spectral peaks at rational multiples of 3 seconds
thegrid <- expand.grid(1:4, 1:4)
thefrac <- thegrid[,1]/thegrid[,2]



### plotting spectral peaks for tags with n>50 as identified earlier!
par(mfrow=c(3,3))
for(i in 1:length(datetime_secs_sub)) {
  thespectra <- spectra(x=datetime_secs_sub[[i]], 
                        lam_min=1, lam_max = 10, n_lam=2000)
  plot(thespectra$lam_trial, thespectra$strength, 
       log="x", type='l', main=names(datetime_secs_sub)[i], #col=0,
       xaxt="n")
  # segments(x0=thespectra$lam_trial, y0=rep(0, length(thespectra$strength)), y1=thespectra$strength)
  axis(side=1, c(3/c(3,2,1), seq(from=3, to=max(thespectra$lam_trial), by=3)))
  abline(h=thespectra$strength_sidak, lty=3)
  axis(side=1, at=3*thefrac, labels=rep("", length(thefrac)))
}


### I want to know what happens when I identify ALL tag codes with ANY periodic component
# okay, let's do it for fish with more than 20 records
datetime_secs_sub2 <- datetime_secs[sapply(datetime_secs, length) > 20]
length(datetime_secs_sub2)
spectra_list <- list()
near3 <- notnear3 <- NA
del <- 0.5
for(i in 1:length(datetime_secs_sub2)) {
  # this is actually pretty stringent
  tempresult <- spectra(datetime_secs_sub2[[i]], 
                        lam_min=1, lam_max=100, n_lam=5000)
  spectra_list[[i]] <- tempresult$lam_trial[tempresult$strength >= 2*tempresult$strength_sidak]
  near3[i] <- length(spectra_list[[i]]) > 0 & any(spectra_list[[i]] > 3-del & spectra_list[[i]] < 3+del)
  notnear3[i] <- length(spectra_list[[i]]) > 0 & !any(spectra_list[[i]] > 3-del & spectra_list[[i]] < 3+del)
  print(i)
}
names(spectra_list) <- names(datetime_secs_sub2)

# subset with a periodic peak near 3 sec
near3_sub <- near3[sapply(spectra_list, length) > 0]

# subset with a periodic peak, but somewhere else!
notnear3_sub <- notnear3[sapply(spectra_list, length) > 0]

# subset with a periodic peak somewhere
spectra_list <- spectra_list[sapply(spectra_list, length) > 0]

length(spectra_list)
sum(near3)
sum(notnear3)
sum(near3_sub)
sum(notnear3_sub)

par(mfrow=c(1,1))
plot(NA, ylim=c(0, length(spectra_list)), xlim=range(unlist(spectra_list)), 
     xlab="period (s)", ylab="",yaxt="n",
     log="x")
for(i in 1:length(spectra_list)) {
  lines(x=range(spectra_list[[i]]), y=rep(i,2), col=i)
  points(x=spectra_list[[i]], y=rep(i, length(spectra_list[[i]])), 
         col=i, pch=1+15*notnear3_sub[i])
}
abline(v=3, lty=1)
abline(v=3*thefrac, lty=3, col=adjustcolor(1, alpha.f=.2))
axis(2, seq_along(spectra_list), labels=names(spectra_list), las=2, cex.axis=.4, col=seq_along(spectra_list))


# plotting spectral peaks for all that have a spike near 3 seconds
par(mfrow=c(3,3))
for(i in which(near3)) {
  thespectra <- spectra(datetime_secs_sub2[[i]], lam_min=1, lam_max=100, n_lam=5000)
  plot(thespectra$lam_trial, thespectra$strength, 
       log="x", type='l', main=names(datetime_secs_sub2)[i], #col=0,
       xaxt="n")
  # segments(x0=thespectra$lam_trial, y0=rep(0, length(thespectra$strength)), y1=thespectra$strength)
  axis(side=1, c(3/c(3,2,1), seq(from=3, to=max(thespectra$lam_trial), by=3)))
  abline(h=thespectra$strength_sidak, lty=3)
  # axis(side=1, at=3*thefrac, labels=rep("", length(thefrac)))
}

# a bit more interesting: which have spike(s) somewhere else
par(mfrow=c(3,2))
for(i in which(notnear3)) {
  thespectra <- spectra(datetime_secs_sub2[[i]], lam_min=1, lam_max=100, n_lam=5000)
  plot(thespectra$lam_trial, thespectra$strength, 
       log="x", type='l', main=names(datetime_secs_sub2)[i], #col=0,
       xaxt="n")
  # segments(x0=thespectra$lam_trial, y0=rep(0, length(thespectra$strength)), y1=thespectra$strength)
  axis(side=1, c(3/c(3,2,1), seq(from=3, to=max(thespectra$lam_trial), by=3)))
  abline(h=thespectra$strength_sidak, lty=3)
  # axis(side=1, at=3*thefrac, labels=rep("", length(thefrac)))
}

# plotting successive differences between pairs of observations
par(mfrow=c(3,2))
for(i in which(notnear3)) {
  # plot(sort(diff(datetime_secs_sub2[[i]])), log="y", ylim=c(1,300))
  plot(diff(sort(datetime_secs_sub2[[i]])), log="y", ylim=c(1,300), yaxt="n",
       main=names(datetime_secs_sub2)[i])
  axis(2, c(3,6,9,30,60,90), las=2, cex.axis=.7)
  # abline(h=3*(1:50), col=adjustcolor(1, alpha.f=.2))
  abline(h=c(3,6,9), col=adjustcolor(4, alpha.f=.5))
  abline(h=c(30,60,90), col=adjustcolor(2, alpha.f=.5))
}

# which fish would have been identified by just the presence of a spectral peak
realfish2 <- sort(names(spectra_list))

either <- union(realfish1, realfish2)
data.frame(either, either %in% realfish1, either %in% realfish2)




### actually, maybe it's more computationally efficient to see if there are
### "significant" spikes in successive differences (e.g. lots of diffs near 3 sec)
### It's a hack, but hey, it might work!
# ------ THIS DOES NOT WORK YET ------# 
diffspikes <- function(x, lam_min=1, lam_max=120, n_lam=1000, log_scale=TRUE) { #, lam_trial=exp(seq(from=0, to=log(10), length.out=1000))) {
  if(log_scale) {
    lam_trial <- exp(seq(from=log(lam_min), to=log(lam_max), length.out=n_lam))
  } else {
    lam_trial <- seq(from=lam_min, to=lam_max, length.out=n_lam)
  }
  
  # tally all diffs, using lam_trial as bins
  xdiffs <- diff(sort(x))
  thetab <- as.numeric(table(cut(xdiffs, breaks=lam_trial)))
  
  # generate a large number of comparable vectors, fully at random
  xnull <- replicate(100, runif(n=length(x), min=min(x), max=max(x)))
  xnull_diff <- apply(xnull, 2, \(xx) diff(sort(xx)))
}