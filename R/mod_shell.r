#!/usr/bin/Rscript

##############################################################################
# Script for running genome-wide models
# Built under R version 3.3.2
##############################################################################
# Setup: load function helper scripts and packages
# The usePackage function loads packages if they already exist,
# otherwise installs from default CRAN repository
##############################################################################
cat("Loading functions and packages...\n")
source("./R/get_functions.r")

packages <- c("tidyverse", "broom", "RColorBrewer", "MASS", "boot", "speedglm",
	"psych", "lmtest", "fmsb", "hexbin", "cowplot", "grid", "gtable", "gridExtra",
	"yaml", "devtools", "openxlsx", "Biostrings", "svglite", "NMF")

sapply(packages, function(x) suppressMessages(usePackage(x)))

# Load predefined color palettes, once RColorBrewer package is loaded
source("./R/palettes.r")

##############################################################################
# Setup: get arguments from _config.yaml file
##############################################################################
args <- yaml.load_file("./_config.yaml")
attach(args)

cat("Script will run with the following parameters:\n")
print(data.frame(n=paste0(names(args), ": ", unlist(args))), right=F)

# Install the bedr package from github, if needed
# install_github('carjed/bedr')

# Load forked bedr package from custom library path specified in config file
suppressMessages(require(bedr, lib.loc=libpath, quietly=T))

# suppressMessages(usePackage(ggbio))

# Install ggrepel from github for advanced options
# install_github('slowkow/ggrepel')
# suppressMessages(require(ggrepel))

# Define additional variables for cleaner strings, etc.
bink <- binw/1000
nbp <- adj*2+1

datadir <- paste0(parentdir,
	"/output/", nbp, "bp_", bink, "k_singletons_", data)

summfile <- paste0(parentdir, "/summaries/full.summary")
singfile <- paste0(parentdir, "/singletons/full.singletons")
bindir <- paste0(parentdir, "/motif_counts/", nbp, "-mers/full")

##############################################################################
# Read and preprocess data
# returns list of 4 elements:
# - full_data$sites (summarized info for each singleton)
# - full_data$bins (raw bin file)
# - full_data$mct (motif counts genome-wide)
# -full_data$aggseq (initial relative mutation rates per K-mer subtype)
##############################################################################
full_data <- getData(summfile, singfile, bindir)

##############################################################################
# Drop singletons in individuals with abnormal mutation signatures
##############################################################################
ptm <- proc.time()
cat("Analyzing sample mutation signatures...\n")
source("./R/ind_sigs.r")
tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

full_data$sites <- full_data$sites %>%
	filter(ID %in% keep_ids$ID)

##############################################################################
# Prepare singleton input for logit model
##############################################################################
if(build_logit){
	cat("Preparing data for logistic regression model...\n")

	mut_cats <- c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA")

	i<-3
	for(chr in 1:22){
		posfile <- paste0(parentdir,
			"/output/logmod_data/chr", chr, "_sites.txt")
		dat <- full_data$sites %>%
		# summfile1 <- sites %>%
			# filter(Category==categ) %>%
			filter(CHR==chr) %>%
			mutate(Type=gsub("cpg_", "", Category2),
				SEQA=substr(Motif, cbp-i, cbp+i),
				SEQB=substr(Motif, cbp*3-i, cbp*3+i),
				Sequence=paste0(SEQA, "(", SEQB, ")")) %>%
			dplyr::select(CHR, POS, Sequence, Type) %>%
			mutate(mut=1) %>%
			spread(Type, mut, fill=0)

		write.table(dat, posfile, col.names=F, row.names=F, quote=F, sep="\t")
	}
}

##############################################################################
# Get relative mutation rates per subtype; plot as heatmap
##############################################################################
ptm <- proc.time()
cat("Analyzing K-mer mutation rates...\n")
source("/R/kmer_analysis.r")
tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

##############################################################################
# Compare 7-mer rates between ERVs and Aggarwala & Voight rates
##############################################################################
ptm <- proc.time()
cat("Comparing with AV model...\n")
source("/R/av_comp.r")
tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

##############################################################################
# Scripts below assume 7-mers+features model is already complete
##############################################################################
sitefile <- paste0(parentdir, "/output/predicted/validation_sites.txt")
if(file.exists(sitefile)){

	# Validation model
	# -must have parentdir/output/rocdat.7bp.2.txt containing list of sites output
	# from this model and pre-annotated with 7-mer rates from
	# 	-MAC10+ (MU_C)
	# 	-ERVs (MU_S)
	# 	-AV (MU_A)
	ptm <- proc.time()
	cat("Validating models on de novo mutations...\n")
	source("./R/validation.r")
	tottime <- (proc.time()-ptm)[3]
	cat("Done (", tottime, "s)\n")

	# Analyze effects of genomic features
	ptm <- proc.time()
	cat("Analyzing genomic features...\n")
	source("./R/coef_summary.r")
	tottime <- (proc.time()-ptm)[3]
	cat("Done (", tottime, "s)\n")

	# Run region-based models
	if(negbin_model){
		cat("Initializing negbin regression model...\n")
		ptm <- proc.time()
		source("./R/negbin_mod.r")
		tottime <- (proc.time()-ptm)[3]
		cat("Done. Finished in", tottime, "seconds \n")
	}
} else {
	cat("Looks like you haven't run the logit models yet.
	Refer to the following README files for instructions:

	-data_mgmt/per-site_dp/README.md
	-data_mgmt/logit_scripts/README.md
	-data_mgmt/process_predicted/README.md \n")
}
