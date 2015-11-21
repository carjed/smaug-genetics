#!/usr/bin/Rscript

##############################################################################
# Script for running genome-wide models
# Built under R version 3.2.2
#
# Usage:
# ./mod_shell.r --summfile="/path/to/summaryfile"
#				--binfile="/path/to/binfile"
#				--binw=binwidth (integer)
# 				--adj=number of bases in either direction to look
# 				--run_agg=TRUE: aggregate data
# 				--pcs=FALSE: use principal components for modeling
# 				--negbin_model==TRUE: run negbin regression
# 				--log_model=FALSE: run logit regression
# 				--run_predict=FALSE: get predicted values
# 				--categ="NN_NN" category to use for logit model
##############################################################################

##############################################################################
# Setup: process args, load function helper script, load packages
##############################################################################
ptm <- proc.time()

parentdir<-dirname(getwd())
cat("Loading functions and packages...\n")

source("get_functions.R")

# Get args from command line; defaults defined below
args <- getArgs(
	defaults=list(adj=3,
		binw=1000000,
		# summfile=paste0(parentdir, "/output/7bp_1000k/full_j.summary"),
		# binfile=paste0(parentdir, "/output/7bp_1000k/full_bin.txt"),
		summfile=paste0(parentdir, "/output/7bp_1000k/chrX.expanded.summary"),
		binfile=paste0(parentdir, "/output/7bp_1000k/chrX.bin_out.txt"),
		run_agg=TRUE,
		pcs=FALSE,
		categ="AT_CG",
		negbin_model=TRUE,
		log_model=FALSE,
		run_predict=FALSE))

# The usePackage function loads packages if they already exist,
# otherwise installs from default CRAN repository
suppressMessages(usePackage(ggplot2))
suppressMessages(usePackage(dplyr))
suppressMessages(usePackage(tidyr))
suppressMessages(usePackage(reshape2))
suppressMessages(usePackage(RColorBrewer))
suppressMessages(usePackage(MASS))
suppressMessages(usePackage(speedglm))
suppressMessages(usePackage(boot))
suppressMessages(usePackage(devtools))
suppressMessages(usePackage(ggbio))

# Manual toggle for installing ggbio package
# Uses the install_github() function from devtools to pull latest version,
# due to issue on some clusters where using
# "source("https://bioconductor.org/biocLite.R")"
# does not properly update the installer
#
# If this is not an issue, can simply run:
# source("https://bioconductor.org/biocLite.R")
# biocLite("ggbio", suppressUpdates=TRUE)
install_ggbio <- 0
if(install_ggbio){
	install_github("Bioconductor/BiocInstaller")
	biocLite("ggbio", suppressUpdates=TRUE)
}

# Get ideogram data for plotting track under chromosome
data(hg19IdeogramCyto, package = "biovizBase")

# Parse args--command line options become objects, so instead of using
# "args$summfile", we can just use "summfile"
# Works as a lightweight version of attach(args), but explicitly defines objects
# in future version, automatically parse numeric variables properly
# Modified from: http://www.r-bloggers.com/extract-objects-from-a-list/

cat("Script will run with the following parameters:\n")
for(i in 1:length(args)){
	##first extract the object value
	tempobj=unlist(args[i])
	varname=names(args[i])

	# optional: print args
	cat(varname, ":", tempobj, "\n")

	##now create a new variable with the original name of the list item
	eval(parse(text=paste(names(args)[i],"= tempobj")))
}
cat("\n")

# Need to manually coerce binwidth and adj args to numeric
binw <- as.numeric(binw)
adj <- as.numeric(adj)

# Define additional variables for cleaner strings, etc.
bink <- binw/1000
nbp <- adj*2+1

datadir <- dirname(summfile)

# Define color palettes
myPaletteCat <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
rb <- c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g <- myPaletteG(6)[1:3]

tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

##############################################################################
# Read in data
##############################################################################
ptm <- proc.time()

if(!file.exists(summfile)){
	cat("Merged summary/bin files do not exist---Merging now...\n")
	combinecmd <- paste0(
		"awk 'FNR==1 && NR!=1{while(/^CHR/) getline; } 1 {print} ' ",
		datadir, "/chr*.expanded.summary > ", datadir, "/full.summary")
	combinecmd2 <- paste0(
		"awk 'FNR==1 && NR!=1{while(/^CHR/) getline; } 1 {print} ' ",
		datadir, "/chr*.bin_out.txt > ", datadir, "/full_bin.txt")
	system(combinecmd)
	system(combinecmd2)
}

cat("Reading summary file:", summfile, "...\n")

summ_5bp_100k <- read.table(summfile, header=F, stringsAsFactors=F, skip=1)
names(summ_5bp_100k)<-c(
	"CHR", "POS", "REF", "ALT", "DP", "AN", "SEQ", "ALTSEQ", "GC")

# summ_5bp_100k <- read.table(summfile, header=F, stringsAsFactors=F)
# names(summ_5bp_100k)<-c(
# 	"CHR", "POS", "REF", "ALT", "DP", "AN", "SEQ", "ALTSEQ", "GC", "SAMPLE")
summ_5bp_100k$BIN <- ceiling(summ_5bp_100k$POS/binw)

cat("Reading bin file:", binfile, "...\n")
bins_5bp_100k <- read.table(binfile, header=T, stringsAsFactors=F, check.names=F)

tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

##############################################################################
# Update data
##############################################################################
ptm <- proc.time()
cat("Updating data...\n")
source("update_dat.r")
dat_5bp_100k <- updateData(summ_5bp_100k, bins_5bp_100k, adj)
rm(summ_5bp_100k)
rm(bins_5bp_100k)
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

##############################################################################
# Run aggregation script to get:
# -genome-wide heatmaps of relative rates
# -negative binomial model to summarize spatial distribution
# -table of relative rates
##############################################################################
if(run_agg){
	ptm <- proc.time()
	cat("Aggregating data...\n")
	source("agg_dat.r")
	aggV <- aggData(dat_5bp_100k, adj) #<-modify the adj value for 3bp data

	agg_5bp_100k <- aggV$oe
	rates5 <- aggV$agg
	summagg2 <- aggV$summagg2

	ratefile <- paste0(parentdir, "/output/", nbp, "bp_", bink, "k_rates.txt")
	write.table(rates5, ratefile, col.names=T, row.names=F, quote=F, sep="\t")

	tottime <- (proc.time()-ptm)[3]
	cat("Done (", tottime, "s)\n")
}

##############################################################################
# Build file of covariates and run PCA
##############################################################################
ptm <- proc.time()
mutcov2file<-paste0(parentdir, "/output/logmod_data/", bink, "kb_mut_cov2.txt")
if(!file.exists(mutcov2file)){
	cat("Building covariate data...\n")
	source("get_covs.r")
} else {
	cat("Reading existing covariate datafile:", mutcov2file, "...\n")
	mut_cov<-read.table(mutcov2file, header=F, stringsAsFactors=F)
}

if(pcs==1){
	names(mut_cov) <- c("CHR", "BIN", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
		"PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13")
} else {
	names(mut_cov) <- c("CHR", "BIN", "H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3",
		"H3K27ac", "H3K27me3", "H3K36me3", "CPGI", "EXON", "TIME", "RATE",
		"prop_GC", "LAMIN")
}

danames <- names(mut_cov)
covnames <- danames[-c(1:2, 14)]

tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

##############################################################################
# Run negative binomial regression models and plot
##############################################################################
if(negbin_model){
	cat("Initializing negbin regression model...\n")
	ptm <- proc.time()
	source("negbin_mod.r")
	tottime <- (proc.time()-ptm)[3]
	cat("Done. Finished in", tottime, "seconds \n")
}

##############################################################################
# Run logistic regression model and obtain predicted probabilities
##############################################################################
if(log_model){
	ptm <- proc.time()
	cat("Initializing logistic regression model...\n")
	source("log_mod.r")
	tottime <- (proc.time()-ptm)[3]
	cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
