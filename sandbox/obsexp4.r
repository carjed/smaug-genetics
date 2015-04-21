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
source("/net/bipolar/jedidiah/mutation/smaug-genetics/get_functions.R")
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

# for(i in 1:120){
	# if()
# }

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
oe2<-data.frame(exp, Category=catrep, BIN=rep(1:(ncol(aggseq.m)-5),6))
oe2<-merge(oe2,ct.ord, by=c("Category","BIN"))
names(oe2)[4]<-"obs"

# get odds ratio for each bin/category
oe2$odds<-oe2$obs/oe2$exp

# order by OR and add column of ranks
oe2.ord<-oe2[order(oe2$Category, oe2$odds),]
# oe2.ord$rk<-rep(1:length(unique(oe2.ord$BIN)),6)
oe2.ord <- ddply(oe2.ord, .(Category), transform, rk=seq_along(Category))

# Plot odds profiles of each category
ggplot(oe2.ord, aes(x=rk, y=odds))+
	geom_point()+
	facet_wrap(~Category)
ggsave("/net/bipolar/jedidiah/mutation/images/odds_profiles.png")

# nbp<-adj*2+1
# oe2.ord$res<-paste0(nbp,"bp")
# data_out<-paste0("oe_",nbp,"bp.txt")
# write.table(oe2.ord, data_out, sep="\t", quote=F, row.names=F, col.names=T)

# oe5<-read.table("oe_5bp.txt", sep="\t", header=T, stringsAsFactors=F)
# oe3<-read.table("oe_3bp.txt", sep="\t", header=T, stringsAsFactors=F)
# oe.full<-rbind(oe3, oe5)

# oe.full2<-oe.full[(oe.full$odds>0.5 & oe.full$odds<2),]
# ggplot(oe.full2, aes(x=rk, y=odds, group=res, colour=res))+geom_line()+facet_wrap(~Category)
# ggsave("3bp_vs_5bp.png")

cats<-unique(aggseq$Category)

# Test for uniformity 
# order by bins
oe2.ord2<-oe2[order(oe2$BIN),]

for(i in 1:6){

	# aggcat<-oe5[oe5$Category==cats[i],]
	aggcat<-oe2.ord2[oe2.ord2$Category==cats[i],]
	aggcat<-aggcat[aggcat$obs>25,]
	test<-chisq.test(aggcat$obs, p=aggcat$exp/sum(aggcat$exp))
	cat<-cats[i]
	print(cat)
	print(test$p.value)
}

# Test for uniformity of 5bp motifs that share a 3bp motif


b<-c("A", "C", "G", "T")

for(i in 1:6){
	
	aggcat<-aggseq[aggseq$Category==cats[i],]
	aggcat<-aggcat[aggcat$freq>25,]
	# test<-chisq.test(aggcat$obs, p=aggcat$exp/sum(aggcat$exp))
	cat<-cats[i]
	msg<-paste0("checking category ",cat)
	print(msg)
	# print(cat)
	# print(test$p.value)

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

oe2.ord2$diff<-oe2.ord2$obs-oe2.ord2$exp

# Plot residual spatial variation after 5bp estimate
ggplot(oe2.ord2, aes(x=BIN, y=diff))+geom_bar(stat="identity", position="identity")+facet_wrap(~Category)
ggsave("/net/bipolar/jedidiah/mutation/images/diffs.png")

# [1] "AT_CG"
# [1] 1.060231e-05
# [1] "AT_GC"
# [1] 9.193585e-60
# [1] "AT_TA"
# [1] 8.038046e-54
# [1] "GC_AT"
# [1] 0
# [1] "GC_CG"
# [1] 5.570329e-33
# [1] "GC_TA"
# [1] 4.950446e-271
# [1] "checking category AT_CG"
# [1] "checking category AT_GC"
# [1] "checking category AT_TA"
# [1] "AT_TA"
# [1] "TAC"
# [1] 0.1405975
# [1] "checking category GC_AT"
# [1] "checking category GC_CG"
# [1] "checking category GC_TA"
# [1] "GC_TA"
# [1] "CCG"
# [1] 0.03218077




