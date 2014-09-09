#' Code to calculate the modern ICC version (using linear mixed-effects model).
#' A lot is borrowed from AFNI's 3dICC_REML.R. Big shout out to the AFNI folks.

#' # Setup
library(lme4)
library(plyr)
library(reshape2)
qc_anat <- read.csv("../data/corr_anat.csv")
qc_func <- read.csv("../data/corr_func.csv")


#' # Anatomical
#' 
#' Examine the reliability for the anatomical measures.
#' 
#' ## Filter Data
#' 
#' We only want subjects with multi-session data (i.e., more than one session of data).
#+ anat-filter1
qc_anat_trt <- ddply(qc_anat, .(subject), function(x) {
  if (nrow(x) == 1) {
    return(data.frame())
  } else {
    return(x)
  }
})
qc_anat_trt$site <- as.factor(as.character(qc_anat_trt$site))

#' Now we check the number of subject per session
#' and for each site, we only keep those with 10 subjects per session.
#+ anat-filter2
table(qc_anat_trt$site, qc_anat_trt$session)
qc_anat_trt <- ddply(qc_anat_trt, .(site, session), function(x) {
  if (nrow(x) < 10) {
    return(data.frame())
  } else {
    return(x)
  }
})
qc_anat_trt$site <- as.factor(as.character(qc_anat_trt$site))

#' For each site, we make the number of subjects per session the same.
#+ anat-filter3
table(qc_anat_trt$site, qc_anat_trt$session)
qc_anat_trt <- ddply(qc_anat_trt, .(site), function(x) {
  # Find the session with the minimum number of subjects
  nsubs_per_sess  <- tapply(x$subject, x$session, length)
  sess_min_subs   <- as.numeric(names(nsubs_per_sess)[which.min(nsubs_per_sess)])
  # Only use those subjects
  subs_to_use     <- x$subject[x$session == sess_min_subs]
  x_to_use        <- subset(x, subject %in% subs_to_use)
  x_to_use
})

#' Finally, a few subjects have more than one scan per site. Keep the first one.
#' There aren't enough to look at within-session reliability.
#+ anat-filter4
table(qc_anat_trt$site, qc_anat_trt$session)
qc_anat_trt <- ddply(qc_anat_trt, .(site, subject, session), function(x) {
  x[1,]
})

#' Now we are done, let's see what we have left
#+ anat-filter5
table(qc_anat_trt$site, qc_anat_trt$session)

#' ## Reshape
#' 
#' We want to do the computations for each site and each measure. So we will need to reshape the matrix.
#+ anat-reshape
qc_anat_trt2 <- melt(qc_anat_trt[,-1], id=c("subject", "session", "site"), variable.name="measure")


#' ## ICC
#' 
#' For each measure and site, we calculate the ICC.
#' We do this with our own code so we can ensure non-negative values amongst other things.
#' 
#+ calc-icc
calc_icc <- function(x) {
  # ICC Juicy Core
  nFact <- 2
  myStat<-vector(mode="numeric", length=nFact+1)
  fmAOV <- lmer(value ~ 1 + (1|subject) + (1|session), x)
  for(i in 1:nFact) myStat[i] <- VarCorr(fmAOV)[[i]][1]  # factor variances
  myStat[nFact+1] <- attr(VarCorr(fmAOV), "sc")^2  # residual variance
  myStat <- myStat/sum(myStat)
  icc <- myStat[1] 
  return(icc)
}
qc_anat_icc <- ddply(qc_anat_trt2, .(measure, site), function(x) {
  # matrix: subject x session
  x$subject <- as.factor(x$subject)
  res <- calc_icc(x)
  nsubjects <- length(unique(x$subject))
  nsessions <- length(unique(x$session))
  c(subjects=nsubjects, sessions=nsessions, icc=res)
})



#' # Functional
#' 
#' Examine the reliability for the functional measures.
#' 
#' ## Filter Data
#' 
#' We only want subjects with multi-session/scan data. Since functional data has
#' many scans within a session, we will make that clean and only look at within session
#' data for the 1st scan.
#+ func-filter1
qc_func_trt <- ddply(qc_func, .(subject), function(x) {
  if (nrow(x) == 1) {
    return(data.frame())
  } else {
    return(x)
  }
})
qc_func_trt$site <- as.factor(as.character(qc_func_trt$site))

# There seems to be some duplicates, remove them
qc_func_trt <- ddply(qc_func_trt, .(subject, session, scan), function(x) x[1,])

#' Now we check the number of subject per session
#' and for each site, we only keep those with 10 subjects per session.
#+ func-filter2
table(qc_func_trt$site, qc_func_trt$session)
qc_func_trt <- ddply(qc_func_trt, .(site, session), function(x) {
  if (nrow(subset(x, scan == 1)) < 10) {
    return(data.frame())
  } else {
    return(x)
  }
})
qc_func_trt$site <- as.factor(as.character(qc_func_trt$site))

#comp <- ddply(qc_func_trt, .(site), function(x) {
#  table(x$session, x$scan)
#})
#comp[is.na(comp)] <- 0
#comp
table(qc_func_trt$site, qc_func_trt$scan)

#' Here we get the between session data.
#' For each site, we make the number of subjects per session the same.
#' Finally, a few subjects have more than one scan per site. Keep the first one.
#+ func-filter3
table(qc_func_trt$site, qc_func_trt$session)
qc_func_trta <- ddply(qc_func_trt, .(site), function(x) {
  # Find the session with the minimum number of subjects
  nsubs_per_sess  <- tapply(x[x$scan == 1,]$subject, x[x$scan == 1,]$session, length)
  sess_min_subs   <- as.numeric(names(nsubs_per_sess)[which.min(nsubs_per_sess)])
  # Only use those subjects
  subs_to_use     <- unique(x$subject[x$session == sess_min_subs])
  x_to_use        <- subset(x, subject %in% subs_to_use)
  x_to_use
})
qc_func_trta$site <- as.factor(as.character(qc_func_trta$site))
qc_func_trta <- ddply(qc_func_trta, .(site, subject, session), function(x) {
  x[1,]
})
qc_func_trta <- ddply(qc_func_trta, .(site), function(x) {
  if (length(unique(x$session)) == 1) {
    return(data.frame())
  } else {
    return(x)
  }
})
qc_func_trta$site <- as.factor(as.character(qc_func_trta$site))

#' Here we get the within session data. We will only look at the first session
#' and keep all the possibly scans for each subject.
#+ func-filter4
table(qc_func_trt$site, qc_func_trt$session)
qc_func_trtb <- ddply(subset(qc_func_trt, session == 1), .(site), function(x) {
  # Find the scan with the minimum number of subjects
  nsubs_per_scan  <- tapply(x$subject, x$scan, length)
  scan_min_subs   <- as.numeric(names(nsubs_per_scan)[which.min(nsubs_per_scan)])
  # Only use that many scan
  subs_to_use     <- unique(x$subject[x$scan == scan_min_subs])
  x_to_use        <- subset(x, subject %in% subs_to_use)
  x_to_use
})
qc_func_trtb <- ddply(qc_func_trtb, .(site), function(x) {
  if (length(unique(x$scan)) == 1) {
    return(data.frame())
  } else {
    return(x)
  }
})
qc_func_trtb$site <- as.factor(as.character(qc_func_trtb$site))

#' Now we are done, let's see what we have left
#+ anat-filter5
table(qc_func_trta$site, qc_func_trta$session)
table(qc_func_trtb$site, qc_func_trtb$scan)


#' ## Reshape
#' 
#' We want to do the computations for each site and each measure. So we will need to reshape the matrix.
#+ anat-reshape
qc_func_trt_btw <- melt(qc_func_trta[,-1], id=c("subject", "session", "scan", "site"), variable.name="measure")
qc_func_trt_win <- melt(qc_func_trtb[,-1], id=c("subject", "session", "scan", "site"), variable.name="measure")



#' ## ICC
#' 
#' For each measure and site, we calculate the ICC.
#+ calc-icc
calc_icc <- function(x) {
  # ICC Juicy Core
  nFact <- 2
  myStat<-vector(mode="numeric", length=nFact+1)
  fmAOV <- lmer(value ~ 1 + (1|subject) + (1|session), x)
  for(i in 1:nFact) myStat[i] <- VarCorr(fmAOV)[[i]][1]  # factor variances
  myStat[nFact+1] <- attr(VarCorr(fmAOV), "sc")^2  # residual variance
  myStat <- myStat/sum(myStat)
  icc <- myStat[1] 
  return(icc)
}
qc_func_icc_btw <- ddply(qc_func_trt_btw, .(measure, site), function(x) {
  # matrix: subject x session
  x$subject <- as.factor(x$subject)
  res <- calc_icc(x)
  nsubjects <- length(unique(x$subject))
  nsessions <- length(unique(x$session))
  c(subjects=nsubjects, sessions=nsessions, icc=res)
})


#' # Save
write.csv(qc_anat_icc, file="qc_anat_icc_btw_sess.csv")
write.csv(qc_func_icc_btw, file="qc_func_icc_btw_sess.csv")
