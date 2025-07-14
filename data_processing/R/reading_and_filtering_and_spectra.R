# library(tidyverse)

### reading data

# maybe if we wrap this in a function, it will be easier to batch-process later

readJSATS <- function(filename) {
  # detecting the regions that correspond to data tables
  as_lines <- readLines(filename)
  starts <- which(substr(as_lines, 1, 8) == "Internal")
  ends <- which(substr(as_lines, 1, 8) == "File End")
  # reading tables and compiling as a list
  raw_data <- list()
  for(i in seq_along(starts)) {
    raw_data[[i]] <- read.csv(filename, 
                              skip=starts[i]+1, 
                              nrows=ends[i]-starts[i]-2,
                              header=FALSE)
    colnames(raw_data[[i]]) <- colnames(read.csv(filename, 
                                                 skip=starts[i]-1, 
                                                 nrows=0,
                                                 header=TRUE))
    # converting DateTime from character 
    raw_data[[i]]$DateTime <- as.POSIXct(raw_data[[i]]$DateTime, 
                                         format="%m/%d/%Y %H:%M:%OS")
    ##### need to figure out how site name is handled??
    # site name seems to have changed between two data tables (only in preamble for first)
    # adding input filename (without file address)
    namesplit <- strsplit(filename, split="/")[[1]]
    raw_data[[i]]$filename <- namesplit[length(namesplit)]
  }
  # apply filtration function?  write this separately
  # write each data table to an external csv file
  # this is temporary for testing
  return(raw_data)
}

# raw_data <- readJSATS(filename = "flat_data/BeaverCreek012_250527_145301.csv")

raw_data <- readJSATS(filename = "experimentation/BeaverCreek012_250527_145301.csv")


# write a separate function to filter, call that in the func above^^


## testing zone

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





## ...could just sort by paste(TagCode,DateTime, and remove entries without a 3-sec neighbor)
asdf$orig_row <- 1:nrow(asdf)
asdf$datetimenum <- as.numeric(asdf$DateTime)
asdf_sort <- asdf
asdf_sort$sortby <- paste(asdf_sort$TagCode, asdf_sort$DateTime)
asdf_sort <- asdf_sort[order(asdf_sort$sortby),]
keepthese1 <- keepthese2 <- rep(F, nrow(asdf_sort))

# nr <-

dt <- 0.5
keepthese1[1:(nrow(asdf_sort)-1)] <- 
  (asdf_sort$TagCode[1:(nrow(asdf_sort)-1)] == asdf_sort$TagCode[2:(nrow(asdf_sort))]) &
  (abs(asdf_sort$datetimenum[1:(nrow(asdf_sort)-1)] - asdf_sort$datetimenum[2:(nrow(asdf_sort))]) < (3+dt)) &
  (abs(asdf_sort$datetimenum[1:(nrow(asdf_sort)-1)] - asdf_sort$datetimenum[2:(nrow(asdf_sort))]) > (3-dt))
keepthese2[2:(nrow(asdf_sort))] <- 
  (asdf_sort$TagCode[1:(nrow(asdf_sort)-1)] == asdf_sort$TagCode[2:(nrow(asdf_sort))]) &
  (abs(asdf_sort$datetimenum[1:(nrow(asdf_sort)-1)] - asdf_sort$datetimenum[2:(nrow(asdf_sort))]) < (3+dt)) &
  (abs(asdf_sort$datetimenum[1:(nrow(asdf_sort)-1)] - asdf_sort$datetimenum[2:(nrow(asdf_sort))]) > (3-dt))
asdf_sort <- subset(asdf_sort, keepthese1 | keepthese2)
length(unique(asdf_sort$TagCode))
table(asdf_sort$TagCode)
nrow(asdf_sort)
plot(tapply(asdf_sort$datetimenum, asdf_sort$TagCode, \(x) diff(range(x)))/60/60/24)#, log="y")
par(mfrow=c(1,1))
# plot(y=as.numeric(as.factor(asdf_sort$TagCode)), x=asdf_sort$DateTime,
#      col=as.numeric(as.factor(asdf_sort$TagCode)))
with(asdf_sort, plot(y=as.numeric(as.factor(TagCode)), x=DateTime,
                     col=as.numeric(as.factor(TagCode))))
segments(x0=tapply(asdf_sort$DateTime, asdf_sort$TagCode, min),
         x1=tapply(asdf_sort$DateTime, asdf_sort$TagCode, max),
         y0=1:length(unique(asdf_sort$TagCode)))

# plot(asdf_sort$DateTime, asdf_sort$Tilt)
with(asdf_sort, plot(DateTime, Tilt))





#### trying a spectral thingy

spectra <- function(x, lam_min=1, lam_max=120, n_lam=1000, log_scale=TRUE) { #, lam_trial=exp(seq(from=0, to=log(10), length.out=1000))) {
  if(log_scale) {
    lam_trial <- exp(seq(from=log(lam_min), to=log(lam_max), length.out=n_lam))
  } else {
    lam_trial <- seq(from=lam_min, to=lam_max, length.out=n_lam)
  }
  strength <- pstrength <- rep(NA, length(lam_trial))
  nbin <- 100
  delta <- 2
  alpha_sidak <- 1-(.99^(1/length(lam_trial)))
  for(i in seq_along(strength)) {
    folded <- x %% lam_trial[i]
    # thedens <- density(folded, n=nbin)
    # whichmaxes <- (which.max(thedens$y) + (-delta):delta) %% nbin
    # whichmaxes[whichmaxes == 0] <- nbin
    # strength[i] <- sum(thedens$y[whichmaxes])/sum(thedens$y)
    # strengthvec <- thedens$y
    thetab <- as.numeric(table(cut(folded, breaks=seq(from=min(folded), to=max(folded), length.out=nbin+1))))
    whichmaxes <- (which.max(thetab) + (-delta):delta) %% nbin
    whichmaxes[whichmaxes == 0] <- nbin
    strengthvec <- thetab
    strength[i] <- sum(strengthvec[whichmaxes])/sum(strengthvec)
    pstrength[i] <- pbinom(q=sum(strengthvec[whichmaxes]), 
                           size=length(x),
                           prob=(1+(2*delta))/nbin,
                           lower.tail=FALSE)
  }
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

thegrid <- expand.grid(1:4, 1:4)
thefrac <- thegrid[,1]/thegrid[,2]

par(mfrow=c(3,3))
for(i in 1:length(datetime_secs_sub)) {
  thespectra <- spectra(x=datetime_secs_sub[[i]], 
                        lam_min=1, lam_max = 10, n_lam=5000)
  plot(thespectra$lam_trial, thespectra$strength, 
       log="x", type='l', main=names(datetime_secs_sub)[i], #col=0,
       xaxt="n")
  # segments(x0=thespectra$lam_trial, y0=rep(0, length(thespectra$strength)), y1=thespectra$strength)
  axis(side=1, c(3/c(3,2,1), seq(from=3, to=max(thespectra$lam_trial), by=3)))
  abline(h=thespectra$strength_sidak, lty=3)
  axis(side=1, at=3*thefrac, labels=rep("", length(thefrac)))
}


# this didn't work like i hoped but it's something to keep trying
datetime_secs_sub2 <- datetime_secs[sapply(datetime_secs, length) > 20]
spectra_list <- list()
near3 <- notnear3 <- NA
del <- 0.5
for(i in 1:length(datetime_secs_sub2)) {
  tempresult <- spectra(datetime_secs_sub2[[i]], lam_min=1, lam_max=100, n_lam=5000)
  spectra_list[[i]] <- tempresult$lam_trial[tempresult$strength >= 2*tempresult$strength_sidak]
  near3[i] <- length(spectra_list[[i]]) > 0 & any(spectra_list[[i]] > 3-del & spectra_list[[i]] < 3+del)
  notnear3[i] <- length(spectra_list[[i]]) > 0 & !any(spectra_list[[i]] > 3-del & spectra_list[[i]] < 3+del)
  print(i)
}
names(spectra_list) <- names(datetime_secs_sub2)
near3_sub <- near3[sapply(spectra_list, length) > 0]
notnear3_sub <- notnear3[sapply(spectra_list, length) > 0]
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
par(mfrow=c(3,2))
for(i in which(notnear3)) {
  # plot(sort(diff(datetime_secs_sub2[[i]])), log="y", ylim=c(1,300))
  plot(diff(sort(datetime_secs_sub2[[i]])), log="y", ylim=c(1,300), yaxt="n")
  axis(2, c(3,6,9,30,60,90), las=2, cex.axis=.7)
  # abline(h=3*(1:50), col=adjustcolor(1, alpha.f=.2))
  abline(h=c(3,6,9), col=adjustcolor(4, alpha.f=.5))
  abline(h=c(30,60,90), col=adjustcolor(2, alpha.f=.5))
}
