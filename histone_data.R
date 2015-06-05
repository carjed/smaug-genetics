getHistBins <- function(file, mark){

	data<-read.table(file, header=F, stringsAsFactors=F, sep="\t")

	names(data) <- c("CHR", "Start", "End", "mark")
	names(data)[4] <- mark
	
	data$CHR <- as.integer(substring(data$CHR, 4))
	data$BIN <- ceiling(data$End/binw)

	# Old version--parse default output from bedops with "|" as delimiter
	# data$count <- as.integer(gsub(".*\\|", "", data$End))
	# data$End <- as.integer(gsub("\\|.*", "", data$End))

	return(data[,c(1,4,5)])
}

# Number of reads per histone mark
H3K4me1<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K4me1.100kb.bed", "H3K4me1")
H3K4me3<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K4me3.100kb.bed", "H3K4me3")
H3K9ac<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K9ac.100kb.bed", "H3K9ac")
H3K9me3<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K9me3.100kb.bed", "H3K9me3")
H3K27ac<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K27ac.100kb.bed", "H3K27ac")
H3K27me3<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K27me3.100kb.bed", "H3K27me3")
H3K36me3<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K36me3.100kb.bed", "H3K36me3")

# Pct in CpG island
CPGI<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/cpg_islands_100kb.bed", "CPGI")
CPGI$CPGI<-CPGI$CPGI/binw

# Pct exonic 
EXON<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/coding_exon_100kb.bed", "EXON")
EXON$EXON<-EXON$EXON/binw

# Replication timing
reptime <- read.table("/net/bipolar/jedidiah/mutation/reference_data/lymph_rep_time.txt", header=F, stringsAsFactors=F, sep="\t")
names(reptime) <- c("CHR", "POS", "TIME")
reptime$BIN <- ceiling(reptime$POS/binw)
rtagg <- aggregate(TIME~CHR+BIN, reptime, mean)

# Recombination rate
rcrate <- read.table("/net/bipolar/jedidiah/mutation/reference_data/recomb_rate.bed", header=T, stringsAsFactors=F, sep="\t")
rcrate$CHR <- as.integer(substring(rcrate$CHR, 4))
rcrate$POS <- (rcrate$START+rcrate$END)/2
rcrate$BIN <- ceiling(rcrate$POS/binw)
rcagg <- aggregate(RATE~CHR+BIN, rcrate, mean)

# GC content
pctgc <- bins_5bp_100k[,c(1,5,4)]
pctgc$CHR <- as.integer(substring(pctgc$CHR, 4))

##############################################################################
# Merge all data and sort
##############################################################################
# Old Version--reduce-merge
# mut_cov <- Reduce(function(x, y) merge(x, y, all=TRUE),
					  # list(H3K4me1, H3K4me3, H3K9ac, H3K9me3, H3K27ac, H3K27me3, H3K36me3, CPGI, EXON))

# New version--cbind histone marks, cpgi, exon, then merge with replication timing and recombination rate
mut_cov<-cbind(H3K4me1, H3K4me3[,2], H3K9ac[,2], H3K9me3[,2], H3K27ac[,2], H3K27me3[,2], H3K36me3[,2], CPGI[,2], EXON[,2])
names(mut_cov)<-c("CHR", "H3K4me1", "BIN", "H3K4me3", "H3K9ac", "H3K9me3", "H3K27ac", "H3K27me3", "H3K36me3", "CPGI", "EXON")
mut_cov<-merge(mut_cov, rtagg, by=c("CHR", "BIN"))
mut_cov<-merge(mut_cov, rcagg, by=c("CHR", "BIN"))
mut_cov<-merge(mut_cov, pctgc, by=c("CHR", "BIN"))

mut_cov <- merge(ct.ord, mut_cov, by=c("CHR", "BIN"))

# mut_cov_sort<-mut_cov[order(mut_cov$CHR, mut_cov$Start),]
# cor(mut_cov_sort$H3K4me1, mut_cov_sort$H3K27ac)
