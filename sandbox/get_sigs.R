##############################################################################
# Read in .mpf data, subset by 6 main categories, and plot 5bp mutation sigs
# using pmsignature packaged developed by Shiraishi et al. (2015)
# http://biorxiv.org/content/early/2015/06/01/019901
##############################################################################

source("http://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges", "BSgenome.Hsapiens.UCSC.hg19"))
require(pmsignature)

positions <- read.table("/net/bipolar/jedidiah/mutation/output/full_summ.mpf", header=F, stringsAsFactors=F)
positions$V2 <- paste0("chr", positions$V2)

positions$CAT <- paste(positions$V4, positions$V5, sep="")

# Manually remove bins near chr20 centromere
# positions <- positions[ which(positions$BIN<260 | positions$BIN>300),]

positions$Category[positions$CAT=="AC" | positions$CAT=="TG"] <- "AT_CG"
positions$Category[positions$CAT=="AG" | positions$CAT=="TC"] <- "AT_GC"
positions$Category[positions$CAT=="AT" | positions$CAT=="TA"] <- "AT_TA"
positions$Category[positions$CAT=="GA" | positions$CAT=="CT"] <- "GC_AT"
positions$Category[positions$CAT=="GC" | positions$CAT=="CG"] <- "GC_CG"
positions$Category[positions$CAT=="GT" | positions$CAT=="CA"] <- "GC_TA"

cats<-unique(positions$Category)
for(i in 1:length(cats)){
	cat<-cats[i]
	pos_dat <- positions[positions$Category==cat, 1:5]
	write.table(pos_dat, "pos_dat.mpf", col.names=F, row.names=F, sep="\t", quote=F)

	G<-readMPFile("/net/bipolar/jedidiah/mutation/smaug-genetics/pos_dat.mpf", numBases=5)

	Param <- getPMSignature(G, K = 1)

	imgfile<-paste0("/net/bipolar/jedidiah/mutation/images/",cat,"_sig.png")
	png(imgfile)
	visPMSignature(Param, 1)
	dev.off()
}