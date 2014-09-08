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
anat_ken <- read.csv("qc_anat_kendall_btw_sess.csv")
func_icc <- read.csv("qc_func_icc_btw_sess.csv")
func_ken <- read.csv("qc_func_kendall_btw_sess.csv")

#+ format
anat_icc$measure <- sub("^anat_", "", anat_icc$measure)
anat_ken$measure <- sub("^anat_", "", anat_ken$measure)
func_icc$measure <- sub("^func_", "", func_icc$measure)
func_ken$measure <- sub("^func_", "", func_ken$measure)
func_icc <- func_icc[func_icc$measure != "num_fd",]
func_ken <- func_ken[func_ken$measure != "num_fd",]

#+ plot-anat-icc
p1 <- ggplot(anat_icc, aes(x=measure, y=icc)) +
  geom_boxplot() + 
  geom_hline(yintercept=0.5, lty=2) + 
  geom_hline(yintercept=0, lty=3) + 
  ylab("ICC") +
  xlab("") + 
  ggtitle("Reliability of Anatomical Measures") +
  theme(axis.title.x      = element_text(family = "Times", face = "plain", 
                                         size=18)) +  
  theme(axis.title.y      = element_text(family = "Times", face = "plain", 
                                         size=18, angle=90, vjust=0.75)) +  
  theme(axis.text.x       = element_text(family = "Times", face = "plain", 
                                         size=16, vjust=0.95, hjust=0.5, angle=0)) + 
  theme(axis.text.y       = element_text(family = "Times", face = "plain", 
                                         size=16, angle=0, hjust=0.5)) + 
  theme(axis.ticks.length = unit(.15, "lines")) + 
  theme(axis.ticks.margin = unit(.15,"lines")) + 
  theme(plot.margin       = unit(c(1, 1, 0.25, 1), "lines"))
plot(p1)
ggsave("anat_icc_btw.png", width=8, height=5, dpi=100)

#+ plot-anat-kendall
p2 <- ggplot(anat_ken, aes(x=measure, y=wt)) +
  geom_boxplot() + 
  geom_hline(yintercept=0.75, lty=2) + 
  geom_hline(yintercept=0.5, lty=3) + 
  ylab("Kendall's W") +
  xlab("") + 
  ylim(c(0,1)) + 
  ggtitle("Consistency of Anatomical Measures") +
  theme(axis.title.x      = element_text(family = "Times", face = "plain", 
                                         size=18)) +  
  theme(axis.title.y      = element_text(family = "Times", face = "plain", 
                                         size=18, angle=90, vjust=0.75)) +  
  theme(axis.text.x       = element_text(family = "Times", face = "plain", 
                                         size=16, vjust=0.95, hjust=0.5, angle=0)) + 
  theme(axis.text.y       = element_text(family = "Times", face = "plain", 
                                         size=16, angle=0, hjust=0.5)) + 
  theme(axis.ticks.length = unit(.15, "lines")) + 
  theme(axis.ticks.margin = unit(.15,"lines")) + 
  theme(plot.margin       = unit(c(1, 1, 0.25, 1), "lines"))
plot(p2)
ggsave("anat_ken_btw.png", width=8, height=5, dpi=100)

#+ plot-func-icc
p3 <- ggplot(func_icc, aes(x=measure, y=icc)) +
  geom_boxplot() + 
  geom_hline(yintercept=0.5, lty=2) + 
  geom_hline(yintercept=0, lty=3) + 
  ylab("ICC") +
  xlab("") + 
  ggtitle("Reliability of Functional Measures") +
  theme(axis.title.x      = element_text(family = "Times", face = "plain", 
                                         size=18)) +  
  theme(axis.title.y      = element_text(family = "Times", face = "plain", 
                                         size=18, angle=90, vjust=0.75)) +  
  theme(axis.text.x       = element_text(family = "Times", face = "plain", 
                                         size=16, vjust=0.95, hjust=0.5, angle=0)) + 
  theme(axis.text.y       = element_text(family = "Times", face = "plain", 
                                         size=16, angle=0, hjust=0.5)) + 
  theme(axis.ticks.length = unit(.15, "lines")) + 
  theme(axis.ticks.margin = unit(.15,"lines")) + 
  theme(plot.margin       = unit(c(1, 1, 0.25, 1), "lines"))
plot(p3)
ggsave("func_icc_btw.png", width=8, height=5, dpi=100)

#+ plot-func-kendall
p4 <- ggplot(func_ken, aes(x=measure, y=wt)) +
  geom_boxplot() + 
  geom_hline(yintercept=0.75, lty=2) + 
  geom_hline(yintercept=0.5, lty=3) + 
  ylab("Kendall's W") +
  xlab("") + 
  ylim(c(0,1)) + 
  ggtitle("Consistency of Functional Measures") +
  theme(axis.title.x      = element_text(family = "Times", face = "plain", 
                                         size=18)) +  
  theme(axis.title.y      = element_text(family = "Times", face = "plain", 
                                         size=18, angle=90, vjust=0.75)) +  
  theme(axis.text.x       = element_text(family = "Times", face = "plain", 
                                         size=16, vjust=0.95, hjust=0.5, angle=0)) + 
  theme(axis.text.y       = element_text(family = "Times", face = "plain", 
                                         size=16, angle=0, hjust=0.5)) + 
  theme(axis.ticks.length = unit(.15, "lines")) + 
  theme(axis.ticks.margin = unit(.15,"lines")) + 
  theme(plot.margin       = unit(c(1, 1, 0.25, 1), "lines"))
plot(p4)
ggsave("func_ken_btw.png", width=8, height=5, dpi=100)
