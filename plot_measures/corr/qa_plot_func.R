#' # CoRR: Functional Data Quality Metrics
#' 
#' Here we plot all the functional QA measures.
#' 
#' ## Load Dependencies
#' 
#' Functions and libraries that are needed for plotting needs
#' Please see https://github.com/preprocessed-connectomes-project/abide/blob/master/plot/qa_plot_functions.R for actual code.
#+ func-spat-source
source("../qa_plot_functions.R")

#' ## Read in Data
#' Along with reading the data, we setup descriptions that will be associated
#' with each column and used as the label for the y-axis.
#+ func-spat-read
df 					<- read.csv("../../data/corr_func.csv")
df$site     <- factor(sub("_", " ", as.character(df$site)))
qa.measures <- c(
  "func_efc", "func_fber", "func_fwhm", "func_gsr", 
  "func_dvars", "func_quality", "func_mean_fd", "func_perc_fd" #"func_outlier"
)
qa.descs    <- list(
  func_efc  = "Entropy Focus Criteria", 
  func_fber = "Foreground to Background\nEnergy Ratio", 
  func_fwhm = "Smoothness of Voxels", 
  func_gsr  = "Ghost to Signal Ratio", 
  func_dvars    = "Standardized DVARS", 
#  func_outlier  = "Fraction of Outlier Voxels", 
  func_quality  = "Mean Distance to Median Volume", 
  func_mean_fd  = "Mean Framewise Displacement (FD)", 
  func_perc_fd  = "Percent FD greater than 0.2mm"
)

#' ## Plot each measure
#' Now we plot the data. Note that here we are removing outliers with values
#' greater than 3 times the IQR relative to the 25% or 75% mark.
#+ func-spat-plot, fig.width=10, fig.height=5, dpi=100
for (measure in qa.measures) {
  desc <- qa.descs[[measure]]
  plot_measure(df, measure, desc, site.col="site", plot=TRUE, 
               outfile=NULL, rm.outlier=TRUE)
}
