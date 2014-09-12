#+ setup
source("../plot_measures/qa_plot_functions.R")
library(ggplot2)
library(grid)
library(RColorBrewer)

#' Let's plot the previously calculated reliability measures
#' We want to show each measure is reliable
#' 
#+ read
anat_icc <- read.csv("qc_anat_icc_btw_sess.csv")
func_icc <- read.csv("qc_func_icc_btw_sess.csv")

#+ format
anat_icc$measure <- sub("^anat_", "", anat_icc$measure)
func_icc$measure <- sub("^func_", "", func_icc$measure)
func_icc <- func_icc[func_icc$measure != "num_fd",]

#+ relabel
anat_icc$measure <- factor(anat_icc$measure, 
                           levels=c("cnr", "efc", "fber", "fwhm", "qi1", "snr"), 
                           labels=c("CNR", "EFC", "FBER", "FWHM", "Qi1", "SNR"))
func_icc$measure <- factor(func_icc$measure, 
                           levels=c( "efc","fber","fwhm","gsr","dvars","quality","mean_fd","perc_fd"), 
                           labels=c("EFC", "FBER", "FWHM", "GSR", "DVARS", "Quality", "Mean FD", "Percent FD"))


#+ plot-anat-icc
p1 <- ggplot(anat_icc, aes(x=measure, y=icc, color=measure, fill=measure)) +
  geom_boxplot(lwd=1.5) + 
  stat_summary(geom="crossbar", width=0.65, fatten=2, color="white", fun.data = function(x) c(y=median(x), ymin=median(x), ymax=median(x))) + 
  ylab("ICC") +
  xlab("") + 
  theme_bw() + 
  theme(panel.border=element_blank()) +
  theme(panel.grid.major.x= element_blank()) + 
  theme(panel.grid.minor  = element_blank()) + 
  theme(panel.grid.major.y= element_line(color="grey70")) + 
  theme(axis.title.x      = element_blank()) +  
  theme(axis.title.y      = element_text(family = "Century Gothic", face = "plain", 
                                         size=18, angle=90, vjust=1)) +  
  theme(axis.text.x       = element_text(family = "Century Gothic", face = "plain", 
                                         size=16, vjust=1, hjust=1, angle=45)) + 
  theme(axis.text.y       = element_text(family = "Times", face = "plain", 
                                         size=16, angle=0, hjust=0.5)) + 
  theme(axis.ticks.y = element_line(color="grey")) +
  theme(axis.ticks.length = unit(.15, "lines")) + 
  theme(axis.ticks.margin = unit(.45,"lines")) + 
  theme(axis.ticks.x      = element_blank()) +
  theme(plot.margin       = unit(c(1, 1, 0.25, 1), "lines")) + 
  theme(legend.position="none")
plot(p1)
ggsave("anat_icc_btw.png", width=4, height=3, dpi=100)

#+ plot-func-icc
p3 <- ggplot(func_icc, aes(x=measure, y=icc, color=measure, fill=measure)) +
  geom_boxplot(lwd=1.5) + 
  stat_summary(geom="crossbar", width=0.65, fatten=2, color="white", fun.data = function(x) c(y=median(x), ymin=median(x), ymax=median(x))) + 
  ylab("ICC") +
  xlab("") + 
  theme_bw() + 
  theme(panel.border=element_blank()) +
  theme(panel.grid.major.x= element_blank()) + 
  theme(panel.grid.minor  = element_blank()) + 
  theme(panel.grid.major.y= element_line(color="grey70")) + 
  theme(axis.title.x      = element_blank()) +  
  theme(axis.title.y      = element_text(family = "Century Gothic", face = "plain", 
                                         size=18, angle=90, vjust=1)) +  
  theme(axis.text.x       = element_text(family = "Century Gothic", face = "plain", 
                                         size=16, vjust=1, hjust=1, angle=45)) + 
  theme(axis.text.y       = element_text(family = "Times", face = "plain", 
                                         size=16, angle=0, hjust=0.5)) + 
  theme(axis.ticks.y = element_line(color="grey")) +
  theme(axis.ticks.length = unit(.15, "lines")) + 
  theme(axis.ticks.margin = unit(.45,"lines")) + 
  theme(axis.ticks.x      = element_blank()) +
  theme(plot.margin       = unit(c(1, 1, 0.25, 1), "lines")) + 
  theme(legend.position="none")
plot(p3)
ggsave("func_icc_btw.png", width=5, height=3, dpi=100)
