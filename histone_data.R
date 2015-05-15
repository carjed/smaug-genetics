##############################################################################
# Read histone data and update columns
##############################################################################
getHistBins <- function(file, mark){

	data<-read.table(file, header=F, stringsAsFactors=F, sep="\t")

	names(data) <- c("CHR", "Start", "End")

	data$count <- as.integer(gsub(".*\\|", "", data$End))
	names(data)[4] <- mark
	data$End <- as.integer(gsub("\\|.*", "", data$End))

	return(data)
}

H3K4me1<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K4me1.100kb.bed", "H3K4me1")
H3K4me3<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K4me3.100kb.bed", "H3K4me3")
H3K9ac<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K9ac.100kb.bed", "H3K9ac")
H3K9me3<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K9me3.100kb.bed", "H3K9me3")
H3K27ac<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K27ac.100kb.bed", "H3K27ac")
H3K27me3<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K27me3.100kb.bed", "H3K27me3")
H3K36me3<-getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K36me3.100kb.bed", "H3K36me3")

##############################################################################
# Merge all data and sort
##############################################################################
histone_dat <- Reduce(function(x, y) merge(x, y, all=TRUE),
					  list(H3K4me1, H3K4me3, H3K9ac, H3K9me3, H3K27ac, H3K27me3, H3K36me3))
histone_dat_sort<-histone_dat[order(histone_dat$CHR, histone_dat$Start),]