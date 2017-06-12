#!/usr/bin/Rscript

##############################################################################
# Logistic regression model
##############################################################################
args <- commandArgs(trailingOnly=TRUE)

catopt <- args[1]
parentdir <- args[2]
libpath <- args[3]
jobid <- args[4]
jobid <- as.numeric(jobid)

.libPaths(c(.libPaths(), libpath))

options(useHTTPS=FALSE)
options(scipen = 8)

invisible(suppressMessages(require(speedglm, quietly=T, warn.conflicts=FALSE)))
invisible(suppressMessages(require(smaug, quietly=T, warn.conflicts=FALSE)))
invisible(suppressMessages(require(dplyr, quietly=T, warn.conflicts=FALSE)))
invisible(suppressMessages(require(boot, quietly=T, warn.conflicts=FALSE)))
invisible(suppressMessages(require(yaml, quietly=T, warn.conflicts=FALSE)))
invisible(suppressMessages(require(bedr, quietly=T, warn.conflicts=FALSE)))

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

sitefile <- paste0(parentdir, "/output/logmod_data/motifs/", escmotif, ".txt")
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


##############################################################################
# Function for running logit model--given input motif,
# writes predicted mutation rates and returns list of coefficient estimates
##############################################################################
logitMod <- function(sites, categ, split){

	# sites$GC <- gcContentCalc(sites_for_GC, 10000, organism=Hsapiens)

	# Run logit model for categories with >10 singletons, return coefficients
	# Otherwise, returns single marginal rate
	predicted <- sites[,1:2]

	coefs <- data.frame()
	mut_cats <- c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA")
	ind <- match(categ, mut_cats)+3

	if(sum(sites[,ind])>10){
		log_mod_formula <- as.formula(paste(categ, "~",
			paste(names(sites)[-(1:9)], collapse="+")))
		log_mod <- speedglm(log_mod_formula, data=sites, family=binomial(), maxit=50)

		predicted$mu <- round(inv.logit(predict(log_mod, newdata=sites)), 6)

		predicted$mu[is.na(predicted$mu)] <- round(sum(sites$mut)/nrow(sites), 6)

		# Get coefficients from model summary, clean up data formats
		coefs <- data.frame(summary(log_mod)$coefficients, stringsAsFactors=F)
		coefs <- cbind(Cov = rownames(coefs), coefs)
		coefs$Cov <- as.character(coefs$Cov)
		rownames(coefs) <- NULL
		names(coefs) <- c("Cov", "Estimate", "SE", "Z", "pval")
		coefs$pval <- as.numeric(as.character(coefs$pval))

		coefs$Sequence <- escmotif

	} else {
		cat("Not enough data--predicted rates will be marginal rate only\n")
		predicted$mu <- round(sum(sites$mut)/nrow(sites), 6)
	}

	# New dir for each chromosome (faster to write, slower to sort?)
	if(split){
		cat(paste0("Writing predicted rates to: ", parentdir, "/output/predicted/", categ, "/chr*/", escmotif, ".txt\n"))
		chr.split <- split(predicted, predicted$CHR)
		for(i in 1:length(chr.split)){
			chr <- unique(chr.split[[i]]$CHR)
			preddir <- paste0(parentdir, "/output/predicted/", categ, "/chr", chr, "/")
			dir.create(preddir, recursive=T, showWarnings = FALSE)

			predfile <- paste0(preddir, categ, "_", escmotif, ".txt")
			write.table(chr.split[[i]], predfile,
				col.names=F, row.names=F, quote=F, sep="\t")
		}
	}

	return(coefs)
}

##############################################################################
# Run models
##############################################################################
cat("Running model on", runmotif, "sites...\n")

for(categ in run_cats){
	coefs <- logitMod(sites=sites, categ=categ, split=TRUE)

	coefdir <- paste0(parentdir, "/output/logmod_data/coefs/", categ, "/")
	dir.create(coefdir, recursive=T, showWarnings = FALSE)
	coeffile <- paste0(coefdir, categ, "_", escmotif, "_coefs.txt")
	cat(paste0("Writing coefficient estimates to: ", coeffile, "\n"))
	write.table(coefs, coeffile, col.names=F, row.names=F, quote=F, sep="\t")
}
