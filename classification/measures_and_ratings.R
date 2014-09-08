#' Here we want to see if we can correctly predict human raters on the quality of each scan
#' given our automated measures.
#' We'll use a logistic regression model here.

# TODO: look into rater 1...who is it and why is there only 1 rating.
# then figure out how to deal with two vs three rating scales for the rater 2 and rater 3

# then see how each measure predicts the outcome measure
# and finally see about using more complicated things

# Read in the ABIDE dataset
# apply x, y, z

library(reshape2)
library(grid)
library(ggplot2)
library(plyr)
library(RColorBrewer)

qc_anat <- read.csv("../data/abide_anat.csv", row.names=1)
qc_func <- read.csv("../data/abide_func.csv", row.names=1)

# Rename columns
names(qc_anat) <- sub("qc_anat_", "qc_", names(qc_anat))
names(qc_func) <- sub("qc_func_", "qc_", names(qc_func))
names(qc_anat) <- sub("qc_rater", "rater", names(qc_anat))
names(qc_func) <- sub("qc_rater", "rater", names(qc_func))
names(qc_anat) <- sub("anat_", "", names(qc_anat))
names(qc_func) <- sub("func_", "", names(qc_func))

# Fix up raters to be consistent
## anat
qc_anat$rater_2[qc_anat$rater_2 == ""] <- NA
qc_anat$rater_3[qc_anat$rater_3 == ""] <- NA
qc_anat$rater_1 <- factor(qc_anat$rater_1, levels=c("fail", "maybe", "OK"), labels=c("fail", "maybe", "ok"))
qc_anat$rater_2 <- factor(qc_anat$rater_2, levels=c("fail", "maybe", "OK"), labels=c("fail", "maybe", "ok"))
qc_anat$rater_3 <- factor(qc_anat$rater_3, levels=c("fail", "maybe", "OK"), labels=c("fail", "maybe", "ok"))


# First let's plot everything relative to each rater

# ## Go with anat first

#' ### setup
id.vars <- c("subject", "site", "rater_1", "rater_2", "rater_3")
measure.vars <- c("cnr", "efc", "fber", "fwhm", "qi1", "snr")
qc_anat2 <- melt(qc_anat, id.vars=id.vars, measure.vars=measure.vars, variable.name="measure")
qc_anat2 <- melt(qc_anat2, measure.vars=id.vars[3:5], variable.name="rater", value.name="score")
qc_anat2$rater <- factor(qc_anat2$rater, levels=c("rater_1", "rater_2", "rater_3"), labels=1:3)
qc_anat2$score <- factor(qc_anat2$score)
head(qc_anat2)

# fix
score <- as.character(qc_anat2$score)
score[score == "maybe"] <- "fail"
qc_anat2$score <- factor(score)

#' ### plot (remove outliers)
source("../plot_measures/qa_plot_functions.R")
df <- subset(qc_anat2)
df <- ddply(df, .(measure), function(x) {
  x <- x[!is.na(x$value),]
  inds <- get_outlier_inds(x$value)
  cat(sum(inds), "\n")
  x[!inds,]
})
# create boxplot that includes outliers
p0 <- ggplot(df, aes(x=rater, y=value, color=score)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_grid(measure ~ ., scales="free_y") + 
  ylab("QC Value") +
  xlab("Rater") + 
  theme(axis.title.x      = element_text(family = "Times", face = "plain", 
                                         size=18)) +  
  theme(axis.title.y      = element_text(family = "Times", face = "plain", 
                                         size=18, angle=90, vjust=0.75)) +  
  theme(axis.text.x       = element_text(family = "Times", face = "plain", 
                                         size=14, vjust=0.95, hjust=1, angle=45)) + 
  theme(axis.text.y       = element_text(family = "Times", face = "plain", 
                                         size=16, angle=0, hjust=0.5)) + 
  theme(axis.ticks.length = unit(.15, "lines")) + 
  theme(axis.ticks.margin = unit(.15,"lines")) + 
  theme(plot.margin       = unit(c(1, 1, 0.25, 1), "lines"))
plot(p0)

# Do regressions or t-tests?
# This tells us how much the score predicts the QA measure
tmp <- ddply(qc_anat2, .(measure, rater), function(x) {
  x2 <- x[!(is.na(x$value) | is.na(x$score)),]
  summary(aov(value ~ score, data=x2))[[1]]$Pr[1]
})
# This tells us how much the QA measure predicts the rater scores
tmp2 <- ddply(qc_anat2, .(measure, rater), function(x) {
  x2 <- x[!(is.na(x$value) | is.na(x$score)),]
  summary(aov(as.numeric(score) ~ value, data=x2))[[1]]$Pr[1]
})

# Let's do a logistic regression (we will remove the maybes)
tmp3 <- ddply(qc_anat2, .(measure, rater), function(x) {
  cat(as.character(x$measure[1]), as.character(x$rater[1]), "\n")
  score <- as.character(x$score)
  score[score == "maybe"] <- "fail"
  x$score <- factor(score)
  x2 <- x[!(is.na(x$value) | is.na(x$score)),]
  mylogit <- glm(score ~ value, data=x2, family="binomial")
  summary(mylogit)$coefficients[2,4]
})

# We can also do this with things together in the same model
# ? put together in the same model
mylogit <- glm(score ~ value, data=tmp2[[18]], family="binomial")
summary(mylogit)$coefficients[2,4]

# Rater 1
## fwhm signif; qi1 and snr somewhat signif; cnr almost signif
x <- qc_anat
x$score <- x$rater_1
score <- as.character(x$score)
score[score == "maybe"] <- "fail"
x$score <- factor(score)
mylogit <- glm(score ~ cnr + efc + fwhm + fber + qi1 + snr, data=x, family="binomial", na.action=na.omit)
summary(mylogit)

# Rater 2
## cnr, efc, fwhm, qi1 signif; fber somewhat signif
x <- qc_anat
x$score <- x$rater_2
score <- as.character(x$score)
score[score == "maybe"] <- "fail"
x$score <- factor(score)
mylogit <- glm(score ~ cnr + efc + fwhm + fber + qi1 + snr, data=x, family="binomial", na.action=na.omit)
summary(mylogit)

mylm <- lm(as.numeric(score) ~ cnr + efc + fwhm + fber + qi1 + snr, data=x, na.action=na.omit)
summary(mylm)


# Rater 3
## efc and qi1 signif; cnr almost signif
x <- qc_anat
x$score <- x$rater_3
score <- as.character(x$score)
score[score == "maybe"] <- "fail"
x$score <- factor(score)
mylogit <- glm(score ~ cnr + efc + fwhm + fber + qi1 + snr, data=x, family="binomial", na.action=na.omit)
summary(mylogit)


## FUNC

# Rater 1
## efc, fber, fwhm, mean_fd signif; perc_fd and gsr somewhat
x <- qc_func
x$score <- x$rater_1
score <- as.character(x$score)
score[score == "maybe"] <- "fail"
x$score <- factor(score)
mylogit <- glm(score ~ efc + fber + fwhm + dvars + quality + mean_fd + perc_fd + gsr, data=x, family="binomial", na.action=na.omit)
summary(mylogit)

# Rater 2
## efc, fber, fwhm, quality, gsr signif; perc_fd somewhat; dvars almost
x <- qc_func
x$score <- x$rater_2
score <- as.character(x$score)
score[score == "maybe"] <- "fail"
x$score <- factor(score)
mylogit <- glm(score ~ efc + fber + fwhm + dvars + quality + mean_fd + perc_fd + gsr, data=x, family="binomial", na.action=na.omit)
summary(mylogit)

# Rater 3
## nothing significant
x <- qc_func
x$score <- x$rater_3
score <- as.character(x$score)
score[score == "maybe"] <- "fail"
x$score <- factor(score)
mylogit <- glm(score ~ efc + fber + fwhm + dvars + quality + mean_fd + perc_fd + gsr, data=x, family="binomial", na.action=na.omit)
summary(mylogit)

#' ### setup
id.vars <- c("subject", "site", "rater_1", "rater_2", "rater_3")
measure.vars <- c("efc", "fber", "fwhm", "dvars", "quality", "mean_fd", "perc_fd", "gsr")
qc_func2 <- melt(qc_func, id.vars=id.vars, measure.vars=measure.vars, variable.name="measure")
qc_func2 <- melt(qc_func2, measure.vars=id.vars[3:5], variable.name="rater", value.name="score")
qc_func2$rater <- factor(qc_func2$rater, levels=c("rater_1", "rater_2", "rater_3"), labels=1:3)
qc_func2$score <- factor(qc_func2$score)
head(qc_func2)

# fix
score <- as.character(qc_func2$score)
score[score == "maybe"] <- "fail"
qc_func2$score <- factor(score)

#' ### plot (remove outliers)
source("../plot_measures/qa_plot_functions.R")
df <- subset(qc_func2)
df <- ddply(df, .(measure), function(x) {
  x <- x[!is.na(x$value),]
  inds <- get_outlier_inds(x$value)
  cat(sum(inds), "\n")
  x[!inds,]
})
# create boxplot that includes outliers
p1 <- ggplot(df, aes(x=rater, y=value, color=score)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_grid(measure ~ ., scales="free_y") + 
  ylab("QC Value") +
  xlab("Rater") + 
  theme(axis.title.x      = element_text(family = "Times", face = "plain", 
                                         size=18)) +  
  theme(axis.title.y      = element_text(family = "Times", face = "plain", 
                                         size=18, angle=90, vjust=0.75)) +  
  theme(axis.text.x       = element_text(family = "Times", face = "plain", 
                                         size=14, vjust=0.95, hjust=1, angle=45)) + 
  theme(axis.text.y       = element_text(family = "Times", face = "plain", 
                                         size=16, angle=0, hjust=0.5)) + 
  theme(axis.ticks.length = unit(.15, "lines")) + 
  theme(axis.ticks.margin = unit(.15,"lines")) + 
  theme(plot.margin       = unit(c(1, 1, 0.25, 1), "lines"))
plot(p1)


library(aod)
wald.test(b = coef(mylogit), Sigma = vcov(mylogit), Terms = 4:6)
