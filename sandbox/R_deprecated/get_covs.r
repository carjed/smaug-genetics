##############################################################################
# Build file of covariates and run PCA
# (previously run from mod_shell.r)
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
# Read in covariate data
##############################################################################

getHistBins <- function(file, mark){

	data <- read.table(file, header=F, stringsAsFactors=F, sep="\t")

	names(data) <- c("CHR", "Start", "End", "mark")
	names(data)[4] <- mark

	data$CHR <- as.integer(substring(data$CHR, 4))
	data$BIN <- ceiling(data$End/binw)

	# Old version--parse default output from bedops with "|" as delimiter
	# data$count <- as.integer(gsub(".*\\|", "", data$End))
	# data$End <- as.integer(gsub("\\|.*", "", data$End))

	return(data[,c(1,4,5)])
}

H3K4me1file <- paste0(parentdir, "/reference_data/histone_marks/binned/H3K4me1.", bink, "kb.bed")
H3K4me3file <- paste0(parentdir, "/reference_data/histone_marks/binned/H3K4me3.", bink, "kb.bed")
H3K9acfile <- paste0(parentdir, "/reference_data/histone_marks/binned/H3K9ac.", bink, "kb.bed")
H3K9me3file <- paste0(parentdir, "/reference_data/histone_marks/binned/H3K9me3.", bink, "kb.bed")
H3K27acfile <- paste0(parentdir, "/reference_data/histone_marks/binned/H3K27ac.", bink, "kb.bed")
H3K27me3file <- paste0(parentdir, "/reference_data/histone_marks/binned/H3K27me3.", bink, "kb.bed")
H3K36me3file <- paste0(parentdir, "/reference_data/histone_marks/binned/H3K36me3.", bink, "kb.bed")

# Number of reads per histone mark
H3K4me1 <- getHistBins(H3K4me1file, "H3K4me1")
H3K4me3 <- getHistBins(H3K4me3file, "H3K4me3")
H3K9ac <- getHistBins(H3K9acfile, "H3K9ac")
H3K9me3 <- getHistBins(H3K9me3file, "H3K9me3")
H3K27ac <- getHistBins(H3K27acfile, "H3K27ac")
H3K27me3 <- getHistBins(H3K27me3file, "H3K27me3")
H3K36me3 <- getHistBins(H3K36me3file, "H3K36me3")

# Pct in CpG island
CPGIfile <- paste0(parentdir, "/reference_data/cpg_islands_", bink, "kb.bed")
CPGI <- getHistBins(CPGIfile, "CPGI")
CPGI$CPGI <- CPGI$CPGI/binw

# Pct exonic
EXONfile <- paste0(parentdir, "/reference_data/coding_exon_", bink, "kb.bed")
EXON <- getHistBins(EXONfile, "EXON")
EXON$EXON <- EXON$EXON/binw

# Replication timing
repfile <- paste0(parentdir, "/reference_data/lymph_rep_time.txt")
reptime <- read.table(repfile, header=F, stringsAsFactors=F, sep="\t")
names(reptime) <- c("CHR", "POS", "TIME")
reptime$BIN <- ceiling(reptime$POS/binw)
rtagg <- aggregate(TIME~CHR+BIN, reptime, mean)

# Recombination rate
rcfile <- paste0(parentdir, "/reference_data/recomb_rate.bed")
rcrate <- read.table(rcfile, header=T, stringsAsFactors=F, sep="\t")
rcrate$CHR <- as.integer(substring(rcrate$CHR, 4))
rcrate$POS <- (rcrate$START+rcrate$END)/2
rcrate$BIN <- ceiling(rcrate$POS/binw)
rcagg <- aggregate(RATE~CHR+BIN, rcrate, mean)

# GC content
pctgc <- dat_5bp_100k$bin[,c(1,5,4)]
pctgc$CHR <- as.integer(substring(pctgc$CHR, 4))

# Lamin B1
lamfile <- paste0(parentdir, "/reference_data/laminB1m.bed")
lam <- read.table(lamfile, header=F, stringsAsFactors=F, sep="\t")
names(lam)<-c("CHR", "START", "END", "VAL")
lam$CHR <- as.integer(substring(lam$CHR, 4))
lam$POS <- (lam$START+lam$END)/2
lam$BIN <- ceiling(lam$POS/binw)
lamagg <- aggregate(VAL~CHR+BIN, lam, mean)

# Merge all data and sort

# Old Version--reduce-merge
# mut_cov <- Reduce(function(x, y) merge(x, y, all=TRUE),
					  # list(H3K4me1, H3K4me3, H3K9ac, H3K9me3, H3K27ac, H3K27me3, H3K36me3, CPGI, EXON))

# New version--cbind histone marks, cpgi, exon, then merge with rep timing and rc rate
mut_cov <- cbind(H3K4me1, H3K4me3[,2], H3K9ac[,2], H3K9me3[,2],
			   H3K27ac[,2], H3K27me3[,2], H3K36me3[,2], CPGI[,2], EXON[,2])
names(mut_cov) <- c("CHR", "H3K4me1", "BIN", "H3K4me3", "H3K9ac", "H3K9me3",
				  "H3K27ac", "H3K27me3", "H3K36me3", "CPGI", "EXON")
mut_cov <- merge(mut_cov, rtagg, by=c("CHR", "BIN"))
mut_cov <- merge(mut_cov, rcagg, by=c("CHR", "BIN"))
mut_cov <- merge(mut_cov, pctgc, by=c("CHR", "BIN"))
mut_cov <- merge(mut_cov, lamagg, by=c("CHR", "BIN"))

if(pcs==1){
	pc.dat <- prcomp(mut_cov[,3:ncol(mut_cov)], center=T, scale=T)
	scores <- as.data.frame(round(pc.dat$x, 3))

	mut_cov<-cbind(mut_cov[,1:2], scores)
}

write.table(mut_cov, mutcov2file, sep="\t", col.names=F, row.names=F, quote=F)
