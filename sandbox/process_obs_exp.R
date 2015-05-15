##############################################################################
# This script processes the singleton summary and bin file for a given
# resolution (3bp or 5bp) and binwidth, and writes a file with
# the expected and observed counts per bin, ranked by O/E odds
#
#   Category BIN       exp obs      odds rk
# 1    AT_CG 296  2.039292   1 0.4903662  1
# 2    AT_CG 261  1.742047   1 0.5740374  2
# 3    AT_CG 445 55.093663  32 0.5808291  3
# 4    AT_CG 597 77.246581  49 0.6343323  4
##############################################################################

##############################################################################
# Process Options/Args
# Define color palettes
# Define functions
##############################################################################
suppressMessages(require(ggplot2))
suppressMessages(require(plyr))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))
suppressMessages(require(grid))

source("/net/bipolar/jedidiah/mutation/smaug-genetics/get_functions.R")

adj<-2
nbp<-adj*2+1
binw<-100000
write<-1

summfile<-paste0("/net/bipolar/jedidiah/mutation/output/chr20.",nbp,"bp.expanded.summary")
binfile<-paste0("/net/bipolar/jedidiah/mutation/output/chr20.bin_out_",nbp,"bp.txt")
chr22 <- read.table(summfile, header=T, stringsAsFactors=F)
bins <- read.table(binfile, header=T, stringsAsFactors=F)

# chr22 <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.expanded.summary", header=T, stringsAsFactors=F)
# bins <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.bin_out.txt", header=T, stringsAsFactors=F)

chr22 <- chr22[-grep(",", chr22$ALT),]

##############################################################################
# Add columns to data
##############################################################################
chr22$BIN <- ceiling(chr22$POS/binw)
chr22$CAT <- paste(chr22$REF, chr22$ALT, sep="")

# Manually remove bins near chr20 centromere
# chr22 <- chr22[ which(chr22$BIN<260 | chr22$BIN>300),]
chr22$Category[chr22$CAT=="AC" | chr22$CAT=="TG"] <- "AT_CG"
chr22$Category[chr22$CAT=="AG" | chr22$CAT=="TC"] <- "AT_GC"
chr22$Category[chr22$CAT=="AT" | chr22$CAT=="TA"] <- "AT_TA"
chr22$Category[chr22$CAT=="GA" | chr22$CAT=="CT"] <- "GC_AT"
chr22$Category[chr22$CAT=="GC" | chr22$CAT=="CG"] <- "GC_CG"
chr22$Category[chr22$CAT=="GT" | chr22$CAT=="CA"] <- "GC_TA"

chr22$Sequence <- ifelse(
	substr(chr22$SEQ,adj+1,adj+1)<substr(chr22$ALTSEQ,adj+1,adj+1),
	paste0(chr22$SEQ,"(",chr22$ALTSEQ,")"),
	paste0(chr22$ALTSEQ,"(",chr22$SEQ,")")
)

# get complement of sequence columns in bin file and remove duplicates
for(i in 5:((4^(adj*2+1))+4)){
	names(bins)[i] <- paste0(names(bins)[i], "(", revcomp(names(bins)[i]), ")" )
}

bins2 <- bins[,names(bins)%in%unique(chr22$Sequence)]
bins <- cbind(bins[,1:4],bins2)
xmax <- floor(max(chr22$BIN)/100)*100

bins2 <- melt(bins[,4:((4^(adj*2+1))/2+4)], id="BIN")
bins2 <- aggregate(data=bins2, value ~ variable, sum)
names(bins2) <- c("Sequence", "COUNT")
bins2$Sequence <- sub("[.]", "(", bins2$Sequence)
bins2$Sequence <- sub("[.]", ")", bins2$Sequence)
bins2$SEQ1 <- substr(bins2$Sequence, 0, adj*2+1)
bins2$SEQ2 <- substr(bins2$Sequence, (adj*2+1)+2, (adj*2+2)+(adj*2+1))
bins2$SEQMIN <- pmin(bins2$SEQ1, bins2$SEQ2)
bins2 <- data.frame(bins2$COUNT, bins2$SEQMIN)
names(bins2) <- c("COUNT", "SEQMIN")

chr22$SEQMIN <- pmin(chr22$SEQ, chr22$ALTSEQ)
chr22 <- merge(chr22, bins2, by="SEQMIN")

aggseq <- count(chr22, c("Sequence", "Category", "COUNT"))
aggseq$rel_prop <- aggseq$freq/aggseq$COUNT



##############################################################################
# Get observed counts for each category/bin
##############################################################################
ct<-count(chr22,c("BIN","Category"))
ct.ord<-ct[order(ct$Category),]

##############################################################################
# Get expected counts for each 
##############################################################################

# Transpose bin data to merge with 
binsT <- setNames(data.frame(t(bins[,-c(1:4)])), paste0("Bin",seq(1,nrow(bins))))
binsT$Sequence<-rownames(binsT)

# merge counts per sequence/category with counts of mutable motifs
aggseq.m<-merge(aggseq, binsT, by="Sequence")

# get expected counts for each row
aggseq.m[,6:ncol(aggseq.m)]<-aggseq.m$rel_prop*aggseq.m[,6:ncol(aggseq.m)]

# sum all sequence combinations for each category/bin
atcg<-colSums(aggseq.m[aggseq.m$Category=="AT_CG",-c(1:5)])
atgc<-colSums(aggseq.m[aggseq.m$Category=="AT_GC",-c(1:5)])
atta<-colSums(aggseq.m[aggseq.m$Category=="AT_TA",-c(1:5)])
gcat<-colSums(aggseq.m[aggseq.m$Category=="GC_AT",-c(1:5)])
gccg<-colSums(aggseq.m[aggseq.m$Category=="GC_CG",-c(1:5)])
gcta<-colSums(aggseq.m[aggseq.m$Category=="GC_TA",-c(1:5)])

# vector of expected counts
exp<-c(atcg, atgc, atta, gcat, gccg, gcta)

# vector of categories
catrep<-c(rep("AT_CG",ncol(aggseq.m)-5), rep("AT_GC",ncol(aggseq.m)-5), rep("AT_TA",ncol(aggseq.m)-5), 
          rep("GC_AT",ncol(aggseq.m)-5), rep("GC_CG",ncol(aggseq.m)-5), rep("GC_TA",ncol(aggseq.m)-5))
		  
# create data frame with observed, expected, and bin
oe2<-data.frame(exp, Category=catrep, BIN=rep(1:(ncol(aggseq.m)-5),6))
oe2<-merge(oe2,ct.ord, by=c("Category","BIN"))
names(oe2)[4]<-"obs"

# get odds ratio for each bin/category
oe2$odds<-oe2$obs/oe2$exp

# order by OR and add column of ranks
oe2.ord<-oe2[order(oe2$Category, oe2$odds),]
oe2.ord <- ddply(oe2.ord, .(Category), transform, rk=seq_along(Category))
oe2.ord$res<-paste0(nbp,"bp")

if(write==1){
	data_out<-paste0("oe_",nbp,"bp.txt")
	write.table(oe2.ord, data_out, sep="\t", quote=F, row.names=F, col.names=T)
}



# if(write==1){
	# data_out<-paste0("oe_",nbp,"bp.txt")
	# write.table(oe1bp.ord, data_out, sep="\t", quote=F, row.names=F, col.names=T)
# }