#' This calculates the collinearity between the different QC measures
#' 
#' We first load the abide and corr dataset
#' then we compute the correlation between the measures
#' and finally we plot that difference
#' 
#+ setup
library(corrplot)
library(RColorBrewer)

#' # ABIDE
#' 
#' ## Setup
#+ abide-setup
# Read
qc_anat <- read.csv("../data/abide_anat.csv", row.names=1)
qc_func <- read.csv("../data/abide_func.csv", row.names=1)
# Fetch relevant columns
anat_cols <- grep("^anat_", names(qc_anat))
func_cols <- grep("^func_", names(qc_func))
# Rename relevant columns for ease
names(qc_anat) <- sub("anat_", "", names(qc_anat))
names(qc_func) <- sub("func_", "", names(qc_func))

#' ## Correlate
#+ abide-correlate
abide_anat <- cor(qc_anat[,anat_cols], use="pairwise.complete.obs")
abide_func <- cor(qc_func[,func_cols], use="pairwise.complete.obs")

#' # CoRR
#' 
#' ## Setup
#+ corr-setup
# Read
qc_anat <- read.csv("../data/corr_anat.csv", row.names=1)
qc_func <- read.csv("../data/corr_func.csv", row.names=1)
# Fetch relevant columns
anat_cols <- grep("^anat_", names(qc_anat))
func_cols <- grep("^func_", names(qc_func))
# Rename relevant columns for ease
names(qc_anat) <- sub("anat_", "", names(qc_anat))
names(qc_func) <- sub("func_", "", names(qc_func))
# Only include session 1 and scan 1
qc_anat <- subset(qc_anat, session == 1)
qc_func <- subset(qc_func, session == 1 & scan == 1)

#' ## Correlate
#+ corr-correlate
corr_anat <- cor(qc_anat[,anat_cols], use="pairwise.complete.obs")
corr_func <- cor(qc_func[,func_cols], use="pairwise.complete.obs")


#' # Visualize
#' 
#' We look at the abide and corr results together.
#' 
#' ## Fixes
#' Since CoRR is missing outliers, we remove that from the abide results
#+ viz-fix-rmcol
rmcol <- which(colnames(abide_func) == "outlier")
abide_func <- abide_func[-rmcol,-rmcol]

#' We also want to ensure the column names are in the same order.
#' The anatomical is actually fine so only do func cols.
#' We match in reference to the abide functional colnames
#+ viz-fix-ordercol
# anat
all(colnames(abide_anat) == colnames(corr_anat))
# func
ref <- colnames(abide_func)
trg <- colnames(corr_func)
inds <- sapply(ref, function(x) which(trg == x))
corr_func0 <- corr_func # just in case
corr_func <- corr_func0[inds,inds]

#' Now we can combine the two together for one big happy family
#'+ viz-fix-combine
# anat
anat <- abide_anat
anat[upper.tri(anat)] <- corr_anat[upper.tri(anat)]
# func (remove the num_fd column since redundant)
func <- abide_func
func[upper.tri(func)] <- corr_func[upper.tri(func)]
rmind <- which(colnames(func) == "num_fd")
func <- func[-rmind,-rmind]

#' ## Plot
#' Generate the corrplot
#+ viz-plot
cols <- rev(brewer.pal(10, "RdBu"))
# the anatomical
corrplot(anat, method="shade", diag=F, outline=T, order="FPC", col=cols, cl.length=length(cols)+1)
quartz.save("all_anat_corrplot.png", width=5, height=5)
# the functional
corrplot(func, method="shade", diag=F, outline=T, order="FPC", col=cols, cl.length=length(cols)+1)
quartz.save("all_func_corrplot.png", width=5, height=5)

#+ viz-plot-thresh
cols <- rev(brewer.pal(10, "RdBu"))
cols[5:6] <- "white"
# the anatomical
corrplot(anat, method="shade", diag=F, outline=T, order="FPC", col=cols, cl.length=length(cols)+1, addgrid.col="grey")
quartz.save("all_thr_anat_corrplot.png", width=5, height=5)
# the functional
corrplot(func, method="shade", diag=F, outline=T, order="FPC", col=cols, cl.length=length(cols)+1, addgrid.col="grey")
quartz.save("all_thr_func_corrplot.png", width=5, height=5)
