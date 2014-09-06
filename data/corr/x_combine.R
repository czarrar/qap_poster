#!/usr/bin/env Rscript

library(plyr)

#' Combines the different QC measures in different files into one (like ABIDE)

#' Read in the data first
#+ read
qc_anat         <- read.csv("FINAL_anat_qc.csv")
qc_func         <- read.csv("qc_filt_epi.csv")
qc_func_motion  <- read.csv("FINAL_func_qc.csv")

#' Select the desired columns
#+ select
qc_anat         <- subset(qc_anat, select=c("subject", "site", "cnr", "efc", "fber", "fwhm", "qi1", "snr"))
qc_func         <- subset(qc_func, select=c("subid", "session", "scan", "site.name", "efc", "fber", "fwhm", "gsr", "dvars", "quality"))
qc_func_motion  <- subset(qc_func_motion, select=c("Subject", "Scan", "site_x", "tr", "MeanFD_Jenkinson", "NumFD_greater_than_0.20", "PercentFD_greater_than_0.20"))

#' Make the columns a bit more consistent
#' subject, session for anat
#+ consistent-anat
tmp             <- strsplit(as.character(qc_anat$subject), "_") # check: all(sapply(tmp, length)==3)
qc_anat$subject <- as.character(sapply(tmp, function(x) x[1]))
qc_anat$session <- as.numeric(sapply(tmp, function(x) x[3]))

#' Have everything have a subject, session, and scan column for func
#' Note that we rename the scans for rockland sample to be 1-3
#+ consistent-func
# Func
qc_func         <- rename(qc_func, c(subid="subject"))
qc_func$subject <- sprintf("%07i", qc_func$subject)
qc_func$scan[qc_func$scan == 645] <- 1
qc_func$scan[qc_func$scan == 1400] <- 2
qc_func$scan[qc_func$scan == 2500] <- 3
# Func Motion
qc_func_motion  <- rename(qc_func_motion, c(Subject="subject", Scan="scan", site_x="site", MeanFD_Jenkinson="mean_fd", NumFD_greater_than_0.20="num_fd", PercentFD_greater_than_0.20="perc_fd"))
tmp             <- strsplit(as.character(qc_func_motion$subject), "_") # check: all(sapply(tmp, length)==3)
qc_func_motion$subject <- as.character(sapply(tmp, function(x) x[1]))
qc_func_motion$session <- as.numeric(sapply(tmp, function(x) x[3]))
tmp             <- as.character(qc_func_motion$scan)
tmp             <- sub("645_1", "1", tmp)
tmp             <- sub("1400_1", "2", tmp)
tmp             <- sub("2500_1", "3", tmp)
tmp2            <- strsplit(tmp, "_") # check: all(sapply(tmp2, length)==3)
qc_func_motion$scan <- as.numeric(sapply(tmp2, function(x) x[2]))

#' Merge the func together
#+ merge-func
tmp             <- merge(qc_func, qc_func_motion)
tmp             <- rename(tmp, c(site="orig_site", site.name="site"))
qc_func         <- tmp

#' Filter
#+ filter
qc_anat         <- subset(qc_anat, select=c("subject", "session", "site", "cnr", "efc", "fber", "fwhm", "qi1", "snr"))
qc_func         <- subset(qc_func, select=c("subject", "session", "scan", "site", "fber", "fwhm", "gsr", "dvars", "quality", "mean_fd", "num_fd", "perc_fd"))

#+ save
write.csv(qc_anat, file="../corr_anat.csv")
write.csv(qc_func, file="../corr_func.csv")


#"anat_cnr","anat_efc","anat_fber","anat_fwhm","anat_qi1","anat_snr",
#"func_efc","func_fber","func_fwhm","func_ghost_x","func_ghost_y",
#"func_dvars","func_outlier","func_quality",
#"func_mean_fd","func_num_fd","func_perc_fd",
#"func_gsr"