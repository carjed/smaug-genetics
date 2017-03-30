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

parentdir <- dirname(getwd())

cat("Loading functions and packages...\n")
scriptdir <- dirname(sys.frame(1)$ofile)
# source(paste0(scriptdir, "/get_functions.r"))
source("./R/get_functions.r")

# Get args from command line; defaults defined below
args <- getArgs(
	defaults=list(adj=4,
		binw=1000000,
		summfile=paste0(parentdir, "/output/9bp_1000k_singletons_full/full.summary"),
		binfile=paste0(parentdir, "/output/9bp_1000k_singletons_full/full_bin.txt"),
		# summfile=paste0(parentdir, "/output/5bp_100k/full.summary"),
		# binfile=paste0(parentdir, "/output/5bp_100k/full_bin.txt"),
		# summfile=paste0(parentdir, "/output/7bp_1000k/chrX.expanded.summary"),
		# binfile=paste0(parentdir, "/output/7bp_1000k/chrX.bin_out.txt"),
		run_agg=TRUE,
		pcs=FALSE,
		categ="AT_CG",
		negbin_model=FALSE,
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
suppressMessages(usePackage(psych))
# Install the bedr package from github, if needed

# install_github('carjed/bedr')
# require(bedr)
# suppressMessages(usePackage(ggbio))

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
# data(hg19IdeogramCyto, package = "biovizBase")

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

	# Change ^SEQ to ^CHR--needed to fix bug in header of common variant data
	combinecmd <- paste0(
		"awk 'FNR==1 && NR!=1{while(/^CHR/) getline; } 1 {print} ' ",
		datadir, "/chr*.expanded.summary > ", datadir, "/full.summary")
	combinecmd2 <- paste0(
		"awk 'FNR==1 && NR!=1{while(/^CHR/) getline; } 1 {print} ' ",
		datadir, "/chr*.bin_out.txt > ", datadir, "/full_bin.txt")
	system(combinecmd)
	system(combinecmd2)
}

cat("Prepping summary file:", summfile, "...\n")

# Data and function for collapsing categories
categs <- c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA")
catdf <- data.frame(Category=categs, stringsAsFactors=F) %>%
	mutate(c1=paste0(substr(Category, 1, 1), substr(Category, 4, 4)),
		c2=paste0(substr(Category, 2, 2), substr(Category, 5, 5))) %>%
		gather(class, sub, c1:c2)

cd2 <- catdf$Category
names(cd2) <- catdf$sub

collapseCat <- function(CAT){
	Category <- cd2[CAT]
	return(Category)
}

# Add category and sequence information
prepSites <- function(summfile){
	sites <- read.table(summfile, header=F, stringsAsFactors=F, skip=1)
	names(sites) <- c("CHR", "POS", "REF", "ALT", "DP", "AN", "SEQ", "ALTSEQ")

	sites <- sites %>%
		mutate(BIN=ceiling(POS/binw),
			Sequence=ifelse(
				substr(SEQ,adj+1,adj+1)<substr(ALTSEQ,adj+1,adj+1),
				paste0(SEQ,"(",ALTSEQ,")"),
				paste0(ALTSEQ,"(",SEQ,")"),
			# CAT=paste0(REF, ALT),
			Category=collapseCat(CAT),

			# SEQMIN=pmin(SEQ, ALTSEQ),
			Category2=ifelse(
				substr(Sequence,adj+1,adj+2)=="CG",
				paste0("cpg_",Category),
				Category)))

	return(sites)
}

sites <- prepSites(summfile)

# Read in motif counts per chromosome
cat("Reading bin file:", binfile, "...\n")
bins <- read.table(binfile, header=T, stringsAsFactors=F, check.names=F)

# summarize motif counts genome-wide
mct <- bins %>%
	dplyr::select(CHR, Sequence=MOTIF, nMotifs=COUNT) %>%
	group_by(Sequence) %>%
	summarise(nMotifs=sum(nMotifs))

# Get relative mutation rates per subtype
cat("Generating motif relative rates...\n")
aggseq <- sites %>%
	group_by(Sequence, Category2, BIN) %>%
	summarise(n=n()) %>%
	summarise(nERVs=sum(n))
aggseq <- merge(aggseq, mct, by="Sequence")
aggseq$rel_prop <- aggseq$nERVs/aggseq$nMotifs
rates5 <- aggseq

cbp <- adj+1
aggseq <- aggseq %>%
	mutate(SEQ1=substr(Sequence, cbp, cbp),
		SEQ3=substr(Sequence, cbp-1, cbp+1),
		SEQ5=substr(Sequence, cbp-2, cbp+2),
		SEQ7=substr(Sequence, cbp-3, cbp+3))

rates7 <- aggseq %>%
	mutate(Category=gsub("cpg_", "", Category2)) %>%
	group_by(Category, SEQ7) %>%
	summarise(nERVs=sum(nERVs),
		nMotifs=sum(nMotifs),
		ERV_rel_rate=nERVs/nMotifs)

ratefile <- paste0(parentdir, "/output/", nbp, "bp_", bink, "k_rates.txt")
write.table(rates5, ratefile, col.names=T, row.names=F, quote=F, sep="\t")

tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

##############################################################################
# Build file of covariates and run PCA
##############################################################################
ptm <- proc.time()
mutcov2file <- paste0(parentdir, "/output/logmod_data/", bink, "kb_mut_cov2.txt")
if(!file.exists(mutcov2file)){
	cat("Building covariate data...\n")
	source("./R/get_covs.r")
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
	source("./R/negbin_mod.r")
	tottime <- (proc.time()-ptm)[3]
	cat("Done. Finished in", tottime, "seconds \n")
}

##############################################################################
# Run logistic regression model and obtain predicted probabilities
##############################################################################
if(log_model){
	ptm <- proc.time()
	cat("Initializing logistic regression model...\n")
	source("./R/log_mod.r")
	tottime <- (proc.time()-ptm)[3]
	cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
