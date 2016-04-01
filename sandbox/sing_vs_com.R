#! /usr/bin/env Rscript

##############################################################################
# Compare correlation between 100kb windowed distribution of singletons
# and common variants.
# -Should be merged with mutation_method_corr.R
##############################################################################

##############################################################################
# Load packages and define color palettes
##############################################################################
suppressMessages(require(ggplot2))
suppressMessages(require(plyr))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
rb<-c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g<-myPaletteG(6)[1:3]
rbg<-c(myPaletteB(6)[1:3],myPaletteR(6)[1:3], myPaletteG(6)[1:3])

##############################################################################
# Read data
##############################################################################
sing<-read.table("/net/bipolar/jedidiah/mutation/output/supp/chr10.cpgi.expanded.summary_S", stringsAsFactors=F)
names(sing)<-c("CHR", "POS", "REF", "ALT", "DP", "AN", "ANNO", "SEQ", "ALTSEQ", "GC", "cpgi")
spec<-read.table("/net/bipolar/jedidiah/mutation/output/supp/chr10.cpgi.expanded.summary_C", stringsAsFactors=F)
names(spec)<-c("CHR", "POS", "REF", "ALT", "DP", "AN", "ANNO", "SEQ", "ALTSEQ", "GC", "cpgi")

sing_bin<-read.table("/net/bipolar/jedidiah/mutation/output/supp/chr10.bin_out.txt_S", header=T, stringsAsFactors=F, check.names=F)
spec_bin<-read.table("/net/bipolar/jedidiah/mutation/output/supp/chr10.bin_out.txt_C", header=T, stringsAsFactors=F, check.names=F)

# chr10<-read.table("chr10.expanded.summary", header=T)
# chr10$gp<-sample(0:1, nrow(chr10), replace=T)

# sing<-chr10[chr10$gp==0,]
# spec<-chr10[chr10$gp==1,]

# sing_bin<-read.table("chr10.bin_out.txt", header=T, check.names=F)
# spec_bin<-sing_bin

##############################################################################
# Add category and bin columns
##############################################################################
binw<-100000
adj<-1
sing$BIN<-ceiling(sing$POS/binw)

chr22<-sing

chr22$Sequence<-ifelse(
	substr(chr22$SEQ,adj+1,adj+1)<substr(chr22$ALTSEQ,adj+1,adj+1),
	paste0(chr22$SEQ,"(",chr22$ALTSEQ,")"),
	paste0(chr22$ALTSEQ,"(",chr22$SEQ,")")
)

# Function to get reverse complement
revcomp = function(DNAstr) {
	step1 = chartr("ACGT","TGCA",DNAstr)
	step2 = unlist(strsplit(step1, split=""))
	step3 = rev(step2)
	step4 = paste(step3, collapse="")
	return(step4)
}

# get complement of sequence columns in bin file and remove duplicates
for(i in 5:((4^(adj*2+1))+4)){
	names(sing_bin)[i]<-paste0(names(sing_bin)[i], "(", revcomp(names(sing_bin)[i]), ")" )
}

sing_bin2<-sing_bin[,names(sing_bin)%in%unique(chr22$Sequence)]
sing_bin<-cbind(sing_bin[,1:4],sing_bin2)

chr22$CAT<-paste(chr22$REF, chr22$ALT, sep="")
chr22$Category[chr22$CAT=="AC" | chr22$CAT=="TG"]<-"AT_CG"
chr22$Category[chr22$CAT=="AG" | chr22$CAT=="TC"]<-"AT_GC"
chr22$Category[chr22$CAT=="AT" | chr22$CAT=="TA"]<-"AT_TA"
chr22$Category[chr22$CAT=="GA" | chr22$CAT=="CT"]<-"GC_AT"
chr22$Category[chr22$CAT=="GC" | chr22$CAT=="CG"]<-"GC_CG"
chr22$Category[chr22$CAT=="GT" | chr22$CAT=="CA"]<-"GC_TA"

chr22$Category<-ifelse((chr22$cpgi==0 & substr(chr22$Sequence,adj+1,adj+2)=="CG"), paste0("cpg_",chr22$Category), chr22$Category)

xmax<-floor(max(chr22$BIN)/100)*100

spec$BIN<-ceiling(spec$POS/binw)
spec$CAT<-paste(spec$REF, spec$ALT, sep="")

spec$Sequence<-ifelse(
	substr(spec$SEQ,adj+1,adj+1)<substr(spec$ALTSEQ,adj+1,adj+1),
	paste0(spec$SEQ,"(",spec$ALTSEQ,")"),
	paste0(spec$ALTSEQ,"(",spec$SEQ,")")
)

# get complement of sequence columns in bin file and remove duplicates
for(i in 5:((4^(adj*2+1))+4)){
	names(spec_bin)[i]<-paste0(names(spec_bin)[i], "(", revcomp(names(spec_bin)[i]), ")" )
}

spec_bin2<-spec_bin[,names(spec_bin)%in%unique(spec$Sequence)]
spec_bin<-cbind(spec_bin[,1:4],spec_bin2)


spec$Category[spec$CAT=="AC" | spec$CAT=="TG"]<-"AT_CG"
spec$Category[spec$CAT=="AG" | spec$CAT=="TC"]<-"AT_GC"
spec$Category[spec$CAT=="AT" | spec$CAT=="TA"]<-"AT_TA"
spec$Category[spec$CAT=="GA" | spec$CAT=="CT"]<-"GC_AT"
spec$Category[spec$CAT=="GC" | spec$CAT=="CG"]<-"GC_CG"
spec$Category[spec$CAT=="GT" | spec$CAT=="CA"]<-"GC_TA"

spec$Category<-ifelse((spec$cpgi==0 & substr(spec$Sequence,adj+1,adj+2)=="CG"), paste0("cpg_",spec$Category), spec$Category)

##############################################################################
# Calculate subtype-specific relative rates per bin
##############################################################################
chr22b<-merge(chr22, sing_bin, by="BIN", all=TRUE)
countsing<-count(chr22b, c("BIN", "AT", "CG", "prop_GC", "Category"))
countsing<-merge(countsing, aggregate(freq~BIN, data=countsing, sum), by="BIN", all=TRUE)
countsing$rel_prop<-countsing$freq.x/countsing$freq.y
countsing<-countsing[order(countsing$BIN),]
countsing$prop<-countsing$freq.x/100000

chr22c<-merge(spec, spec_bin, by="BIN", all=TRUE)
countspec<-count(chr22c, c("BIN", "AT", "CG", "prop_GC", "Category"))
countspec<-merge(countspec, aggregate(freq~BIN, data=countspec, sum), by="BIN", all=TRUE)
countspec$rel_prop<-countspec$freq.x/countspec$freq.y
countspec<-countspec[order(countspec$BIN),]
countspec$prop<-countspec$freq.x/100000

names(countsing)[names(countsing)=="prop"]<-"sing_prop"
names(countspec)[names(countspec)=="prop"]<-"spec_prop"

names(countsing)[names(countsing)=="freq.x"]<-"sing_count"
names(countspec)[names(countspec)=="freq.x"]<-"spec_count"

##############################################################################
# Merge and sort combined data
##############################################################################
fullcount<-merge(countsing, countspec, by=c("BIN", "Category"))
fullcount<-fullcount[order(fullcount$BIN),]
fullcount<-fullcount[complete.cases(fullcount),]

##############################################################################
# Calculate correlation in each subtype and create df
##############################################################################
cors<-round(c(
		# cor(fullcount[fullcount$Category=="AT_CG",]$sing_prop, fullcount[fullcount$Category=="AT_CG",]$spec_prop),
		# cor(fullcount[fullcount$Category=="AT_GC",]$sing_prop, fullcount[fullcount$Category=="AT_GC",]$spec_prop),
		# cor(fullcount[fullcount$Category=="AT_TA",]$sing_prop, fullcount[fullcount$Category=="AT_TA",]$spec_prop),
		# cor(fullcount[fullcount$Category=="GC_AT",]$sing_prop, fullcount[fullcount$Category=="GC_AT",]$spec_prop),
		# cor(fullcount[fullcount$Category=="GC_CG",]$sing_prop, fullcount[fullcount$Category=="GC_CG",]$spec_prop),
		# cor(fullcount[fullcount$Category=="GC_TA",]$sing_prop, fullcount[fullcount$Category=="GC_TA",]$spec_prop)
		cor(fullcount[fullcount$Category=="AT_CG",]$sing_count, fullcount[fullcount$Category=="AT_CG",]$spec_count),
		cor(fullcount[fullcount$Category=="AT_GC",]$sing_count, fullcount[fullcount$Category=="AT_GC",]$spec_count),
		cor(fullcount[fullcount$Category=="AT_TA",]$sing_count, fullcount[fullcount$Category=="AT_TA",]$spec_count),
		cor(fullcount[fullcount$Category=="GC_AT",]$sing_count, fullcount[fullcount$Category=="GC_AT",]$spec_count),
		cor(fullcount[fullcount$Category=="GC_CG",]$sing_count, fullcount[fullcount$Category=="GC_CG",]$spec_count),
		cor(fullcount[fullcount$Category=="GC_TA",]$sing_count, fullcount[fullcount$Category=="GC_TA",]$spec_count),
		cor(fullcount[fullcount$Category=="cpg_GC_AT",]$sing_count, fullcount[fullcount$Category=="cpg_GC_AT",]$spec_count),
		cor(fullcount[fullcount$Category=="cpg_GC_CG",]$sing_count, fullcount[fullcount$Category=="cpg_GC_CG",]$spec_count),
		cor(fullcount[fullcount$Category=="cpg_GC_TA",]$sing_count, fullcount[fullcount$Category=="cpg_GC_TA",]$spec_count)
	), 2)

cors<-paste0("r=",cors)

xpos<-rep(median(fullcount$sing_count), 9)
ypos<-0.9*rep(max(fullcount$sing_count), 9)





fullcount$Category <- factor(fullcount$Category, levels = c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA", "cpg_GC_AT", "cpg_GC_CG", "cpg_GC_TA"))
# cats<-sort(unique(fullcount$Category))

dat<-data.frame(xpos, ypos, levels(fullcount$Category), cors)
names(dat)<-c("x", "y", "Category", "val")

# dat$Category <- factor(dat$Category, levels = c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA", "cpg_GC_AT", "cpg_GC_CG", "cpg_GC_TA"))

##############################################################################
# Scatter plot of relative rates in singleton vs. divergent sites
##############################################################################
# my_grob = grobTree(textGrob("This text stays in place!", x=0.1,  y=0.95, hjust=0,
  # gp=gpar(col="blue", fontsize=12, fontface="italic")))

ggplot(data=fullcount,
		aes(x=sing_count, y=spec_count, colour=Category, group=Category))+
	geom_point(alpha=0.4)+
	geom_smooth(method=lm, se=FALSE, colour="black")+
	scale_colour_manual(values=rbg)+
	xlab("Singletons")+
	ylab("Common Variants (MAC>10)")+
	# xlab("Group 1 Relative Rate")+
	# ylab("Group 2 Relative Rate")+
	theme_bw()+
	theme(legend.position="none")+
	geom_text(data=dat, aes(x=-Inf, y=Inf, label=val, hjust=0, vjust=1), colour="black")+
	facet_wrap(~Category, scales="free")
ggsave("/net/bipolar/jedidiah/mutation/images/sing_com_corr.png", width=6.2, height=4.2)

ggplot()+
	geom_point(data=fullcount[fullcount$Category=="AT_GC",], aes(x=sing_count, y=spec_count, colour=prop_GC.y, size=5), alpha=0.6)+
	scale_colour_gradientn(colours=myPalette(9))+
	xlab("Singletons")+
	ylab("Common Variants (MAC>10)")+
	# xlab("Group 1 Relative Rate")+
	# ylab("Group 2 Relative Rate")+
	theme_bw()+
	geom_text(data=dat, aes(x=-Inf,y=Inf, label=val, hjust=0, vjust=1))+
ggsave("/net/bipolar/jedidiah/mutation/images/sing_com_corr_bgc.png", width=12.4, height=8.4)
