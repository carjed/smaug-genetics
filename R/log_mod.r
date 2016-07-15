#!/usr/bin/Rscript

##############################################################################
# Logistic regression model
##############################################################################

.libPaths(c(.libPaths(), "/exports/home/jedidiah/R/x86_64-pc-linux-gnu-library/2.13"))

options(useHTTPS=FALSE)
suppressMessages(require(speedglm, quietly=T))
suppressMessages(require(bedr, quietly=T))
suppressMessages(require(dplyr, quietly=T))
suppressMessages(require(BSgenome.Hsapiens.UCSC.hg19, quietly=T))
suppressMessages(require(Repitools, quietly=T))
suppressMessages(require(boot, quietly=T))

source("./get_functions.r")

args <- getArgs(
	defaults=list(
		jobid=25,
		categ="AT_CG",
		nmotifs=4096,
		nodes=10))

nbp <- 7
parentdir <- "/net/bipolar/jedidiah/mutation"

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

jobid <- as.numeric(jobid)

# cluster <- makeCluster(nodes, type = "SOCK", outfile="/net/bipolar/jedidiah/mutation/snow.log")
# registerDoSNOW(cluster)

# Target mutation rate
mu <- 1e-8
nbases <- 5.8e9
nind <- 3612
indmuts <- mu*nbases
nsing <- 10000
numgens <- ceiling(nsing/indmuts)
meangens <- mean(1:numgens)
denom <- nind*meangens

# Fast list of 6 basic categories from agg_5bp_100k data
mut_cats <- c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA")

# subset reference data to only AT or GC bases
catopt <- substr(categ,0,2)

motiffile <- "/net/bipolar/jedidiah/mutation/output/7bp_1000k_rates.txt"
incmd <- paste0("cut -f1-3 ", motiffile)
motifdat <- read.table(pipe(incmd), header=T, stringsAsFactors=F)
motifs <- motifdat %>%
	mutate(Category=gsub("cpg_", "", Category2)) %>%
	filter(Category==categ) %>%
	dplyr::select(Sequence) %>%
	unlist

runmotif<-motifs[jobid]

# comb <- function(x, ...) {
#       mapply(rbind,x,...,SIMPLIFY=FALSE)
# }

##############################################################################
# Run models
##############################################################################
cat("Running model on", motif, "sites...\n")
coefs <- logitMod(motif=runmotif, nbp=nbp, parentdir=parentdir, categ=categ)

# covlist <- clusterApply(cluster, motifs[1:nmotifs], logitMod, nbp=nbp, parentdir=parentdir, categ=categ)
# fullcoef <- rbind_all(covlist)
escmotif <- substr(runmotif, 0, nbp)

coefdir <- paste0(parentdir, "/output/logmod_data/coefs/", categ, "/")
coeffile <- paste0(coefdir, categ, "_", escmotif, "_coefs.txt")
write.table(coefs, coeffile, col.names=F, row.names=F, quote=F, sep="\t")

##############################################################################
# Function for running logit model--given input motif,
# writes predicted mutation rates and returns list of coefficient estimates
##############################################################################
logitMod <- function(motif, nbp, parentdir, categ){

	escmotif <- substr(motif, 0, nbp)

	# source("./get_functions.r")

	# Merge per-chromosome motif files to single file
	# Define name of temporary file for motif i
	# sitefile <- paste0(parentdir, "/output/logmod_data/motifs/", categ, "/",
	# 	categ, "_", escmotif, ".txt")
	sitefile <- paste0(parentdir, "/output/logmod_data/motifs/", categ, "/dp/",
			categ, "_", escmotif, "_dp.txt")

	if(!(file.exists(sitefile))){
		cat("Merging ", motif, " files...\n")

		perchrtmp <- paste0(parentdir,
			"/output/logmod_data/chr*/chr*_", categ, "_", motif, ".txt")

		catcmd1 <- paste0("find ", parentdir, "/output/logmod_data/chr* -name '*",
			escmotif, "*.txt' | sort -V | xargs cat >> ", sitefile)
		system(catcmd1)
	}

	incmd <- paste0("cut -f1-5,18,19 ", sitefile)
	sites <- read.table(pipe(incmd), header=F, stringsAsFactors=F)
	names(sites) <- c("POS", "CHR", "BIN", "Sequence", "mut", "EXON", "DP")

	# Initialize data for calculating GC content
	# sites_for_GC <- data.frame(position=sites$POS, chr=paste0("chr", sites$CHR))

	# Add histone marks to site data
	hists <- c("H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3", "H3K27ac", "H3K27me3", "H3K36me3")
	dflist <- list()
	for(i in 1:length(hists)){
	  mark <- hists[i]
	  file <- paste0(parentdir, "/reference_data/histone_marks/broad/sort.E062-", mark, ".bed")
	  hist <- binaryCol(sites, file)
	  dflist[[i]] <- hist
	}

	df <- as.data.frame(do.call(cbind, dflist))
	names(df) <- hists
	sites <- cbind(sites, df)

	# Add other features
	sites$CpGI <- binaryCol(sites,
		"/net/bipolar/jedidiah/mutation/reference_data/cpg_islands_sorted.bed")
	sites$RR <- rcrCol(sites,
		"/net/bipolar/jedidiah/mutation/reference_data/recomb_rate.bed")
	sites$LAMIN <- binaryCol(sites,
		"/net/bipolar/jedidiah/mutation/reference_data/lamin_B1_LADS2.bed")
	sites$DHS <- binaryCol(sites,
		"/net/bipolar/jedidiah/mutation/reference_data/DHS.bed")
	sites$TIME <- repCol(sites,
		"/net/bipolar/jedidiah/mutation/reference_data/lymph_rep_time.txt")
	sites$GC <- gcCol(sites, "/net/bipolar/jedidiah/mutation/output/3bp_10k/full_bin.txt")
	# sites$GC <- gcContentCalc(sites_for_GC, 10000, organism=Hsapiens)

	# Run logit model for categories with >10 singletons, return coefficients
	# Otherwise, returns single marginal rate
	predicted <- sites[,c(2,1,3)]
	coefs <- data.frame()
	if(sum(sites$mut)>10){
		log_mod_formula <- as.formula(paste("mut~",
			paste(names(sites)[-(1:5)], collapse="+")))
		log_mod <- speedglm(log_mod_formula, data=sites, family=binomial(), maxit=50)

		predicted$mu <- round(inv.logit(predict(log_mod, newdata=sites)), 6)

		# Get coefficients from model summary, clean up data formats
		coefs <- data.frame(summary(log_mod)$coefficients, stringsAsFactors=F)
		coefs <- cbind(Cov = rownames(coefs), coefs)
		coefs$Cov <- as.character(coefs$Cov)
		rownames(coefs) <- NULL

		coefs[,-1] <- data.frame(apply(coefs[,-1], 2,
			function(x) as.numeric(as.character(x))))
		names(coefs) <- c("Cov", "Estimate", "SE", "Z", "pval")
		coefs$Sequence <- escmotif

	} else {
		cat("Not enough data--using marginal rate only\n")
		predicted$mu <- round(sum(sites$mut)/nrow(sites), 6)
	}

# New dir for each bin (slower to write, faster to sort)
	# chr.split<-split(sites, sites$CHR)
	# for(i in 1:length(chr.split)){
	# 	bin.split<-split(chr.split[[i]], chr.split[[i]]$BIN)
	# 	for(j in 1:length(bin.split)){
	# 		chr<-unique(bin.split[[j]]$CHR)
	# 		bin<-unique(bin.split[[j]]$BIN)
	# 		preddir <- paste0(parentdir,
					# "/output/predicted/", categ, "/chr", chr, "/bin", bin, "/")
	# 		dir.create(preddir, recursive=T)
	#
	# 		predfile<-paste0(preddir, categ, "_", escmotif, ".txt")
			# write.table(bin.split[[j]], predfile,
			# 	col.names=F, row.names=F, quote=F, sep="\t")
	# 	}
	# }

# New dir for each chromosome (faster to write, slower to sort?)
	chr.split<-split(predicted, predicted$CHR)
	for(i in 1:length(chr.split)){
		chr<-unique(chr.split[[i]]$CHR)
		preddir <- paste0(parentdir, "/output/predicted/", categ, "/chr", chr, "/")
		dir.create(preddir, recursive=T)

		predfile<-paste0(preddir, categ, "_", escmotif, ".txt")
		write.table(chr.split[[i]], predfile,
			col.names=F, row.names=F, quote=F, sep="\t")
	}

	# list(coefs, predicted) #<- use to output predicted data as list element (high mem. usage!)
	# list(coefs)
	return(coefs)
}
