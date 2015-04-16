suppressMessages(require(ggplot2))
suppressMessages(require(plyr))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))
suppressMessages(require(grid))
adj<-1
myPaletteCat <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
rb<-c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g<-myPaletteG(6)[1:3]
chr22 <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.3bp.expanded.summary", header=T, stringsAsFactors=F)
bins <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.bin_out_3bp.txt", header=T, stringsAsFactors=F)
chr22 <- chr22[-grep(",", chr22$ALT),]
source("get_functions.R")
binw<-100000
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
# z<-as.matrix(bins[,5:((4^(adj*2+1))/2+4)]) %*% aggseq2$rel_prop
# zdf<-data.frame(exp=c(z), BIN=seq(1,631))
# oe<-merge(ct, zdf, "BIN")
# oe$oe<-oe$freq/oe$exp
# names(oe)<-c("BIN", "obs", "exp", "ratio")
# oe.ord<-oe[order(oe$ratio),]
# oe.ord$rk<-seq(1,599)
# ggplot(oe.ord, aes(x=rk, y=ratio))+geom_point()+facet_wrap(~Category)
# ggsave("3bpres.png")
# oe$res<-"3bp"
# write.table(oe, "oe_3bp.txt", sep="\t", quote=F, row.names=F, col.names=T)
##############################################################################
# Get expected counts for each 
##############################################################################
# Transpose bin data to merge with 
binsT = setNames(data.frame(t(bins[,-c(1:4)])), paste0("Bin",seq(1,nrow(bins))))
binsT$Sequence<-rownames(binsT)
# merge counts per sequence/category with counts of mutable motifs
aggseq.m<-merge(aggseq, binsT, by="Sequence")
# get expected counts for each row
aggseq.m[,6:ncol(aggseq.m)]<-aggseq.m$rel_prop*aggseq.m[,6:ncol(aggseq.m)]
# sum all sequence combinations f
atgc<-colSums(aggseq.m[aggseq.m$Cat
atcg<-colSums(aggseq.m[aggseq.m$Categor
atta<-colSums(aggseq.m[aggseq.m
gcat<-
gccg<-colSums(aggseq.m[aggseq.m$Category=="GC_CG",-c(1:5)])
gcta<-colSums(
# vector of expected counts
exp<-c(atcg, atgc, atta, gcat, gccg, gcta)
# 
catrep<-c(rep("AT_CG",ncol(aggseq.m)-5), rep("AT_GC",ncol(aggseq.m)-5), 
          rep("GC_AT",ncol(aggseq.m)-5), rep("GC_CG",ncol(aggseq.m)-5),
  
# create data frame with observed, expected, and bin
oe2<-data.frame
oe2<-merge(oe2,ct.ord, by=c("Category","BIN"
names(oe2)[4]<-"obs"
# get odds ratio for each bin/category
oe2$odds<-oe2$obs/oe2$exp
# order by OR and add column of ranks
oe2.ord<-oe2[order(oe2$Category, oe2$odds),]
oe2.ord$rk<-rep(1:length(unique(oe2.ord$BIN)),6)
# Plot OR profiles of each cat
ggplot(oe2.ord, aes(x=rk, y=or))+geom_
nbp<-adj*2+1
oe2.ord$res
data_out<-pa
suppressMessages(require(ggplot2))
suppressMessages(require(plyr))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))
suppressMessages(require(grid))
adj<-1
myPaletteCat <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
rb<-c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g<-myPaletteG(6)[1:3]
chr22 <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.3bp.expanded.summary", header=T, stringsAsFactors=F)
bins <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.bin_out_3bp.txt", header=T, stringsAsFactors=F)
chr22 <- chr22[-grep(",", chr22$ALT),]
source("get_functions.R")
binw<-100000
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
# Transpose bin data to merge with 
binsT = setNames(data.frame(t(bins[,-c(1:4)])), paste0("Bin",seq(1,nrow(bins))))
binsT$Sequence<-rownames(binsT)
# merge counts per sequence/category with counts of mutable motifs
aggseq.m<-merge(aggseq, binsT, by="Sequence")
# get expected counts for each row
aggseq.m[,6:ncol(aggseq.m)]<-aggseq.m$rel_prop*aggseq.m[,6:ncol(aggseq.m)]
# sum all sequence combinations for each category/bin
atgc<-colSums(aggseq.m[aggseq.m$Category=="AT_GC",-c(1:5)])
atcg<-colSums(aggseq.m[aggseq.m$Category=="AT_CG",-c(1:5)])
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
oe2.ord$rk<-rep(1:length(unique(oe2.ord$BIN)),6)
# Plot OR profiles of each category
ggplot(oe2.ord, aes(x=rk, y=or))+geom_point()+facet_wrap(~Category)
nbp<-adj*2+1
oe2.ord$res<-paste0(nbp,"bp")
data_out<-paste0("oe_",nbp,"bp.txt")
write.table(oe2.ord, data_out, sep="\t", quote=F, row.names=F, col.names=T)
dim(oe2.ord)
head(oe2.ord)
3592/6
oe2.ord$rk <- with(oe2.ord, ave(Category, Category, odds, FUN=seq))
oe2.ord$rk <- with(oe2.ord, ave(Category, odds, FUN=seq))
require(plyr)
oe2.ord <- ddply(oe2.ord, .(Category), transform, rk=seq_along(Category))
head(oe2.ord)
tail(oe2.ord)
max(oe2.ord$rk)
ggplot(oe2.ord, aes(x=rk, y=or))+geom_point()+facet_wrap(~Category)
nbp<-adj*2+1
oe2.ord$res<-paste0(nbp,"bp")
data_out<-paste0("oe_",nbp,"bp.txt")
write.table(oe2.ord, data_out, sep="\t", quote=F, row.names=F, col.names=T)
oe5<-read.table("oe_5bp.txt", sep="\t", header=T, stringsAsFactors=F)
oe3<-read.table("oe_3bp.txt", sep="\t", header=T, stringsAsFactors=F)
oe.full<-rbind(oe3, oe5)
oe.full2<-oe.full[(oe.full$or>1.1),]
ggplot(oe.full, aes(x=rk, y=odds, group=res, colour=res))+geom_line()+facet_wrap(~Category)
oe5<-read.table("oe_5bp.txt", sep="\t", header=T, stringsAsFactors=F)
oe3<-read.table("oe_3bp.txt", sep="\t", header=T, stringsAsFactors=F)
oe.full<-rbind(oe3, oe5)
oe.full2<-oe.full[(oe.full$or>1.1),]
ggplot(oe.full, aes(x=rk, y=odds, group=res, colour=res))+geom_line()+facet_wrap(~Category)
oe5<-read.table("oe_5bp.txt", sep="\t", header=T, stringsAsFactors=F)
head(oe5)
head(oe3)
suppressMessages(require(ggplot2))
suppressMessages(require(plyr))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))
suppressMessages(require(grid))
adj<-2
myPaletteCat <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
rb<-c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g<-myPaletteG(6)[1:3]
chr22 <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.5bp.expanded.summary", header=T, stringsAsFactors=F)
bins <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.bin_out_5bp.txt", header=T, stringsAsFactors=F)
chr22 <- chr22[-grep(",", chr22$ALT),]
source("get_functions.R")
binw<-100000
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
# z<-as.matrix(bins[,5:((4^(adj*2+1))/2+4)]) %*% aggseq2$rel_prop
# zdf<-data.frame(exp=c(z), BIN=seq(1,631))
# oe<-merge(ct, zdf, "BIN")
# oe$oe<-oe$freq/oe$exp
# names(oe)<-c("BIN", "obs", "exp", "ratio")
# oe.ord<-oe[order(oe$ratio),]
# oe.ord$rk<-seq(1,599)
# ggplot(oe.ord, aes(x=rk, y=ratio))+geom_point()+facet_wrap(~Category)
# ggsave("3bpres.png")
# oe$res<-"3bp"
# write.table(oe, "oe_3bp.txt", sep="\t", quote=F, row.names=F, col.names=T)
##############################################################################
# Get expected counts for each 
##############################################################################
# Transpose bin data to merge with 
binsT = setNames(data.frame(t(bins[,-c(1:4)])), paste0("Bin",seq(1,nrow(bins))))
binsT$Sequence<-rownames(binsT)
# merge counts per sequence/category with counts of mutable motifs
aggseq.m<-merge(aggseq, binsT, by="Sequence")
# get expected counts for each row
aggseq.m[,6:ncol(aggseq.m)]<-aggseq.m$rel_prop*aggseq.m[,6:ncol(aggseq.m)]
# sum all sequence combinations for each category/bin
atgc<-colSums(aggseq.m[aggseq.m$Category=="AT_GC",-c(1:5)])
atcg<-colSums(aggseq.m[aggseq.m$Category=="AT_CG",-c(1:5)])
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
oe2<-data.frame(exp, Category=catrep, BI
oe2<-merge(oe2,ct.ord, by=c("Category"
names(oe2)[4]<-"obs"
# g
oe2$odds<-oe
# order by OR and add column of ranks
oe2.ord<-oe2[order(oe2$Category, oe2$odd
# oe2.ord$rk<-rep(1:l
oe2.ord <- ddply(oe2.ord, .(Category), transform, rk=seq_along(Category))
# P
ggplot(oe2.ord, aes(x=rk, y=odds))+g
nbp<-adj*2+1
oe2.ord$res<-paste0(nbp,"bp")
data
head(oe2.ord)
nbp<-adj*2+1
oe2.ord$res<-paste0(nbp,"bp")
data_out<-paste0("oe_",nbp,"bp.txt")
write.table(oe2.ord, data_out, sep="\t", quote=F, row.names=F, col.names=T)
oe5<-read.table("oe_5bp.txt", sep="\t", header=T, stringsAsFactors=F)
oe3<-read.table("oe_3bp.txt", sep="\t", header=T, stringsAsFactors=F)
oe.full<-rbind(oe3, oe5)
oe.full2<-oe.full[(oe.full$or>1.1),]
ggplot(oe.full, aes(x=rk, y=odds, group=res, colour=res))+geom_line()+facet_wrap(~Category)
ggsave("3bp_vs_5bp.png")
tail(oe5[oe5$Category=="GC_CG",])
oe.full2<-oe.full[(oe.full$odds>0.5 & oe.full$odds<2),]
ggplot(oe.full2, aes(x=rk, y=odds, group=res, colour=res))+geom_line()+facet_wrap(~Category)
ggsave("3bp_vs_5bp.png")
choose(16,2)
96*120
0.05/11520
ls()
dim(oe.full2)
dim(oe.full)
dim(aggseq.m)
dim(aggseq)
head(aggseq)
choose(4,2)
6*4
24*2
48*96
120*96


# Test for uniformity of 5bp motifs that share a 3bp motif
b<-c("A", "C", "G", "T")
cats<-unique(aggseq$Category)

for(i in 1:6){

	aggcat<-aggseq[aggseq$Category==cats[i],]
	# print(head(aggcat))

	for(j in 1:4){
		for(k in 1:4){
			ref<-unique(substr(aggcat$Sequence,3,3))
			motif<-paste0(b[j],ref,b[k])
			dat<-aggcat[substr(aggcat$Sequence,2,4)==motif,]
			
			dat$exp<-dat$COUNT*(sum(dat$freq)/sum(dat$COUNT))
			# print(head(dat))
			cat<-cats[i]

			test<-chisq.test(dat$freq, p=dat$exp/sum(dat$exp))
			if(test$p.value>0.05/96){
				print(cat)
				print(motif)
				print(test$p.value)
				# print(dat)
			}
		}
	}
}

aggseq.atgc<-aggseq[aggseq$Category=="AT_GC",]
aggseq.atgc.AT<-aggseq.atgc[substr(aggseq.atgc$Sequence,2,4)=="AAT",]
aggseq.atgc.CA<-aggseq.atgc[substr(aggseq.atgc$Sequence,2,4)=="CAA",]

# Calculate expected counts
aggseq.atgc.CA$exp<-aggseq.atgc.CA$COUNT*(sum(aggseq.atgc.CA$freq)/sum(aggseq.atgc.CA$COUNT))

chisq.test(freq, p=exp/sum(exp))
savehistory("test_5bp.r")
