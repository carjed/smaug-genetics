#!/usr/bin/Rscript

##############################################################################
# Script for running genome-wide models
# Built under R version 3.3.2
##############################################################################

##############################################################################
# Setup: load functions, packages, and arguments
##############################################################################
cat("Loading packages...\n")
# source("./R/get_functions.r")
# source("./R/palettes.r")

# load yaml package to parse args
install.packages("yaml", quiet=TRUE)
require(yaml)

# parse args
args <- yaml.load_file("./_config.yaml")
attach(args)

# add user libpath, if it exists (useful if running on a different node)
.libPaths(c(.libPaths(), libpath))

# install/load devtools packages
require(devtools)
install_github('carjed/smaug', quiet=TRUE)
install_github('slowkow/ggrepel', quiet=TRUE)
gh_packages <- c("smaug", "ggrepel")
invisible(sapply(gh_packages, function(x)
	suppressMessages(require(x, character.only=TRUE))))

# install/load CRAN packages
packages <- c("tidyverse", "broom", "RColorBrewer", "MASS", "boot", "speedglm",
	"psych", "lmtest", "fmsb", "hexbin", "cowplot", "grid", "gtable", "gridExtra",
	"yaml", "openxlsx", "Biostrings", "svglite", "NMF")
invisible(sapply(packages, function(x)
	suppressMessages(usePackage(x))))

# load Bioconductor packages
# suppressMessages(usePackage(ggbio))

cat("Script will run with the following parameters:\n")
print(data.frame(ARGLIST = paste0(names(args), ": ", unlist(args))), right=F)

# Define additional variables for cleaner strings, etc.
bink <- binw/1000
nbp <- adj*2+1

datadir <- paste0(parentdir,
	"/output/", nbp, "bp_", bink, "k_singletons_", data)

summfile <- paste0(parentdir, "/summaries/", mac, ".", data, ".summary")
singfile <- paste0(parentdir, "/singletons/full.singletons")
bindir <- paste0(parentdir, "/motif_counts/", nbp, "-mers/full")

# Read and preprocess data
full_data <- getData(summfile, singfile, bindir)
gc()

##############################################################################
# Drop singletons in individuals with abnormal mutation signatures
##############################################################################
cat("Analyzing sample mutation signatures...\n")
timefun("./R/ind_sigs.r")

full_data$sites <- full_data$sites %>%
	filter(ID %in% ped$ID)
dim(full_data$sites)

full_data$sites <- full_data$sites %>%
	filter(ID %in% keep_ids$ID)
dim(full_data$sites)

##############################################################################
# Get relative mutation rates per subtype; plot as heatmap
##############################################################################
cat("Analyzing K-mer mutation rates...\n")
cat("Calculating K-mer rates...\n")
full_data$aggseq <- get_aggseq(full_data$sites, full_data$mct)
timefun("./R/kmer_analysis.r")

##############################################################################
# Compare full data to strict ancestral data
##############################################################################
cat("Comparing full data with high-confidence ancestral data...\n")
timefun("./R/ancestral.r")

# Choose to drop ERVs with mispolarized/unknown ancestral state?
# Recalculates rates
if(drop_mispolarized){
	full_data$sites <- sites_c_hc
	full_data$aggseq <- get_aggseq(full_data$sites, full_data$mct)
	timefun("./R/kmer_analysis.r")
}

##############################################################################
# Compare full data to masked data
##############################################################################
cat("Comparing full data with strict mask...\n")
timefun("./R/full_vs_mask.r")

# Additional filtering?

##############################################################################
# Compare 7-mer rates between ERVs and Aggarwala & Voight rates
##############################################################################
cat("Comparing with AV model...\n")
timefun("./R/av_comp.r")

##############################################################################
# Prepare singleton input for logit model
##############################################################################
if(build_logit){
	cat("Preparing data for logistic regression model...\n")

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

	cat("Validating models on de novo mutations...\n")
	timefun("./R/validation.r")

	# Analyze effects of genomic features
	cat("Analyzing genomic features...\n")
	timefun("./R/coef_summary.r")

	# Run region-based models
	if(negbin_model){
		cat("Initializing negbin regression model...\n")
		timefun("./R/negbin_mod.r")
	}
} else {
	cat("Looks like you haven't run the logit models yet.
	Refer to the following README files for instructions:

	-data_mgmt/per-site_dp/README.md
	-data_mgmt/logit_scripts/README.md
	-data_mgmt/process_predicted/README.md \n")
}
