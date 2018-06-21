#!/usr/bin/Rscript

##############################################################################
# Logistic regression model
##############################################################################
args <- commandArgs(trailingOnly=TRUE)
# args <- c("AT", "/net/bipolar/jedidiah/mutation", "/net/snowwhite/home/jedidiah/R/x86_64-pc-linux-gnu-library/3.3", 16)
catopt <- args[1]
parentdir <- args[2]
libpath <- args[3]
jobid <- args[4]
jobid <- as.numeric(jobid)

.libPaths(c(.libPaths(), libpath))

options(useHTTPS=FALSE)
options(scipen = 8)

suppressPackageStartupMessages(library("speedglm", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("smaug", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("dplyr", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("boot", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("yaml", quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library("bedr", quietly=TRUE, warn.conflicts=FALSE))

# source(paste0(parentdir, "/smaug-genetics/R/get_functions.r"))

yaml_args <- yaml.load_file(paste0(parentdir, "/smaug-genetics/_config.yaml"))
attach(yaml_args)

nbp_run <- 7

# Fast list of 6 basic categories from agg_5bp_100k data
mut_cats <- c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA")

# subset reference data to only AT or GC bases
# catopt <- substr(categ,0,2)

run_cats <- mut_cats[grepl(paste0("^", catopt), mut_cats)]

motiffile <- paste0(parentdir, "/output/7bp_1000k_rates.txt")
incmd <- paste0("cut -f1-3 ", motiffile)
motifdat <- read.table(pipe(incmd), header=T, stringsAsFactors=F)
motifs <- motifdat %>%
	mutate(Category=gsub("cpg_", "", Category2)) %>%
	filter(grepl(paste0("^", catopt), Category)) %>%
	dplyr::select(Sequence) %>%
	unlist %>% unique

runmotif <- motifs[jobid]
escmotif <- substr(runmotif, 0, nbp_run)

# stopif(file.exists(paste0(parentdir, "output/logmod_data/coefs")))

sitefile <- paste0(parentdir, "/output/logmod_data/motifs3/", escmotif, ".txt")
sites <- read.table(sitefile, header=F, stringsAsFactors=F)
names(sites) <- c("CHR", "POS", "Sequence", mut_cats, "DP")
sites <- sites %>%
	arrange(CHR, POS)

# Add histone marks to site data
hists <- c("H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3",
	"H3K27ac", "H3K27me3", "H3K36me3")
dflist <- list()
for(i in 1:length(hists)){
  mark <- hists[i]
  file <- paste0(parentdir, "/reference_data/sort.E062-",
		mark, ".bed")
  hist <- binaryCol(sites, file)
  dflist[[i]] <- hist
}

df <- as.data.frame(do.call(cbind, dflist))
names(df) <- hists
sites <- cbind(sites, df)

# Add other features
sites$EXON <- binaryCol(sites,
	paste0(parentdir, "/reference_data/GRCh37_RefSeq_sorted.bed"))
sites$CpGI <- binaryCol(sites,
	paste0(parentdir, "/reference_data/cpg_islands_sorted.bed"))
sites$RR <- rcrCol(sites,
	paste0(parentdir, "/reference_data/recomb_rate.bed"))
sites$LAMIN <- binaryCol(sites,
	paste0(parentdir, "/reference_data/lamin_B1_LADS2.bed"))
sites$DHS <- binaryCol(sites,
	paste0(parentdir, "/reference_data/DHS.bed"))
sites$TIME <- repCol(sites,
	paste0(parentdir, "/reference_data/lymph_rep_time.txt"))
sites$GC <- gcCol(sites,
	paste0(parentdir, "/reference_data/gc10kb.bed"))

write.table(sites[,!(names(mtcars) %in% c("CHR", "POS", "Sequence"))],
	paste0(parentdir, "/output/logmod_data/annotated/", escmotif, "_annotated.txt"))