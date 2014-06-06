#! /usr/bin/env Rscript

#####################################################
#Plots mutation types as proportions
#1. Relative to all corresponding nucleotides per bin
#2. Relative to observed total per bin
#####################################################

##############################################################################
# Process Options/Args
##############################################################################

options(warn=-1)
sink("R.log")

suppressMessages(require(ggplot2))
suppressMessages(require(plyr))
suppressMessages(require(RColorBrewer))

args<-commandArgs(TRUE)

chr<-as.character(args[1])
macl<-as.character(args[2])
binw<-as.numeric(args[3])
cpg_flag<-as.character(args[4])
summ<-as.character(args[5])

if (macl=="singletons") mac<-"Singleton"
if (macl=="doubletons") mac<-"Doubleton"

##############################################################################
# Initialize titles/output strings
##############################################################################

main_rel_title<-print(paste0("Chr",chr," ",mac, " Relative Mutation Rate"))
main_rel_out<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_mutation_prop.png"))

main_scale_title<-print(paste0("Chr",chr," ",mac, " Mutations--scaled proportions"))
main_scale_out<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_mutation_prop2.png"))

main_dist_title<-print(paste0("Chr",chr," ",mac, " Distribution by Mutation Type"))
main_ident_out<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_dist_ident.png"))
main_dodge_out<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_dist_dodge.png"))

##############################################################################
#read in summary file and bins file
#
#summ<-print(paste0("/net/bipolar/jedidiah/bcftools/summaries/",macl,"/all/chr",chr,".",macl,".summary.txt"))
#chr22<-read.table("/net/bipolar/jedidiah/testpipe/summaries/chr22.summary", header=F, stringsAsFactors=F)
#names(chr22)<-c("CHR", "POS", "REF", "ALT", "DP", "AN", "ANNO")
##############################################################################

chr22<-read.table(summ, header=T, stringsAsFactors=F)
bins<-read.table("bin_out.txt", header=T, stringsAsFactors=F)

chr22$DP<-as.numeric(chr22$DP)
chr22$DP[is.na(chr22$DP)]=mean(chr22$DP, na.rm=T)
chr22$AVGDP<-chr22$DP/(chr22$AN/2)
chr22$BIN<-ceiling(chr22$POS/binw)
chr22$CAT<-paste(chr22$REF, chr22$ALT, sep="")

chr22$Category[chr22$CAT=="AC" | chr22$CAT=="TG"]<-"AT to CG"
chr22$Category[chr22$CAT=="AG" | chr22$CAT=="TC"]<-"AT to GC"
chr22$Category[chr22$CAT=="AT" | chr22$CAT=="TA"]<-"AT to TA"
chr22$Category[chr22$CAT=="GA" | chr22$CAT=="CT"]<-"GC to AT"
chr22$Category[chr22$CAT=="GC" | chr22$CAT=="CG"]<-"GC to CG"
chr22$Category[chr22$CAT=="GT" | chr22$CAT=="CA"]<-"GC to TA"

##############################################################################
# Heatmaps
##############################################################################

depth_heat_out<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_mutation_vs_depth_heatmap.png"))
hotspot_heat_out<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_mutation_vs_hotspot_heatmap.png"))

#long form aggregate--average depth per bin per category
aggdata<-aggregate(AVGDP ~ BIN+Category, data=chr22, mean)

#subset to exclude bins with very high avg. depth
agg2<-aggdata[aggdata$BIN>250,]

#Graphics--heatmap of avg. depth per bin for the 6 mutation categories
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

xmax<-floor(max(aggdata$BIN)/100)*100
ggplot(agg2, aes(x=BIN, y=Category, fill=AVGDP))+geom_tile()+scale_fill_gradientn(colours=myPalette(4))+scale_x_continuous(breaks=seq(0,xmax,50))
ggsave(depth_heat_out)

hotspot_agg<-aggregate(DIST ~ BIN+Category, data=chr22, mean)

ggplot(hotspot_agg, aes(x=BIN, y=Category, fill=DIST))+geom_tile()+scale_fill_gradientn(colours=myPalette(4))+scale_x_continuous(breaks=seq(0,xmax,50))
ggsave(hotspot_heat_out)

##############################################################################
# Local Sequence Analysis
##############################################################################

chr22$Sequence<-paste0(pmin(chr22$SEQ, chr22$ALTSEQ),"(",pmax(chr22$SEQ, chr22$ALTSEQ),")")
cats<-factor(chr22$Category)

for (i in 1:6) {
	cat<-levels(cats)[i]
	chr22s<-chr22[chr22$Category==cat,]
	title<-paste0("Chr22: ",cat," Mutations by Local Sequence")
	ggplot(chr22s, aes(x=POS))+geom_histogram(binwidth=binw, position="identity")+ggtitle(title)+facet_wrap(~Sequence)
	out<-paste0("/net/bipolar/jedidiah/images/chr",chr," ",cat," by local seq.png")
	ggsave(out)
}

aggseq<-count(chr22, c("Sequence", "Category"))
aggseq_a<-aggseq[grep("^A", aggseq$Category),]
aggseq_g<-aggseq[grep("^G", aggseq$Category),]

at_seq_title<-paste0("Chr",chr,": AT Mutations by Local Sequence")
gc_seq_title<-paste0("Chr",chr,": GC Mutations by Local Sequence")
at_seq_out<-paste0("/net/bipolar/jedidiah/images/chr",chr,"_AT_seq.png")
gc_seq_out<-paste0("/net/bipolar/jedidiah/images/chr",chr,"_GC_seq.png")

ggplot(aggseq_a, aes(x=Category, y=freq, fill=Sequence))+geom_bar(position="dodge", stat="identity")+ggtitle(at_seq_title)
ggsave(at_seq_out)
ggplot(aggseq_g, aes(x=Category, y=freq, fill=Sequence))+geom_bar(position="dodge", stat="identity")+ggtitle(gc_seq_title)
ggsave(gc_seq_out)

##############################################################################
#CpG Analysis
##############################################################################
if (cpg_flag=="on") {
	#names(chr22a)<-c("POS","PAIR","CPGI")
	cpg_scale_title<-print(paste0("Chr",chr," ",mac, " CpG Mutations--scaled proportions"))
	cpg_scale_out<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_cpg_mutation_prop2.png"))

	cpg_dist_title<-print(paste0("Chr",chr," ",mac, " CpG Distribution by Mutation Type"))
	cpg_ident_out<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_cpg_dist_ident.png"))
	cpg_dodge_out<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_cpg_dist_dodge.png"))
	
	chr22m$CpG[chr22m$PAIR=="CG" | chr22m$PAIR=="GC"]<-1
	chr22m$CpG[chr22m$PAIR!="CG" & chr22m$PAIR!="GC"]<-0

	chr22cpg<-chr22m[(chr22m$CpG==1 & chr22m$CPGI==0),]
	chr22cpg$CAT<-paste(chr22cpg$REF, chr22cpg$ALT, sep="")

	chr22cpg$Category[chr22cpg$CAT=="GA" | chr22cpg$CAT=="CT"]<-"GC to AT"
	chr22cpg$Category[chr22cpg$CAT=="GC" | chr22cpg$CAT=="CG"]<-"GC to CG"
	chr22cpg$Category[chr22cpg$CAT=="GT" | chr22cpg$CAT=="CA"]<-"GC to TA"
	
	#Plot distribution of 3 CpG categories
	ggplot(chr22cpg, aes(x=POS, colour=Category, fill=Category, group=Category, alpha=0.5))+geom_histogram(binwidth=binw, position="identity")+ggtitle(cpg_dist_title)
	ggsave(cpg_ident_out)
	ggplot(chr22cpg, aes(x=POS, colour=Category, fill=Category, group=Category, alpha=0.5))+geom_histogram(binwidth=binw, position="dodge")+ggtitle(cpg_dist_title)
	ggsave(cpg_dodge_out)

	#CpG--Merge bins + summary files and process
	countcpg<-count(chr22cpg, c("Category", "BIN", "CHR"))
	count2cpg<-merge(countcpg, bins, by="BIN")
	count3cpg<-count2cpg[grep("^GC", count2cpg$Category),]
	count3cpg$prop<-count3cpg$freq/count3cpg$CG
	
	count4cpg<-aggregate(freq~BIN, data=count3cpg, sum)
	names(count4cpg)<-c("BIN", "total")
	count5cpg<-merge(count4cpg, count3cpg, by="BIN")
	count5cpg$rel_prop<-count5cpg$freq/count5cpg$total
	count5cpg<-count5cpg[order(count5cpg$BIN, count5cpg$Category),]
	
	suppressMessages(ggplot(count5cpg, aes(x=factor(BIN), y=rel_prop, colour=Category, fill=Category, alpha=0.5))+geom_bar(position="stack", stat="identity")+ scale_x_discrete(breaks=seq(0,xmax,50))+xlab("Bin")+ylab("Proportion")+ggtitle(cpg_scale_title))
	suppressMessages(ggsave(cpg_scale_out))
} 

#All--Merge bins + summary files and process
count<-count(chr22, c("Category", "BIN", "CHR"))
count2<-merge(count, bins, by="BIN")
countAT<-count2[grep("^AT", count2$Category),]
countGC<-count2[grep("^GC", count2$Category),]
countAT$prop<-countAT$freq/countAT$AT
countGC$prop<-countGC$freq/countGC$CG
count3<-rbind(countAT, countGC)

#Plot proportion relative to all possible sites in each bin
suppressMessages(ggplot(count3, aes(x=factor(BIN), y=prop, colour=Category, fill=Category, group=Category, alpha=0.5))+geom_bar(position="dodge", stat="identity")+ scale_x_discrete(breaks=seq(0,xmax,50))+xlab("Bin")+ylab("Proportion")+ggtitle(main_rel_title))
suppressMessages(ggsave(main_rel_out))

#Obtain proportion of each mutation category relative to total hits in each bin 
count4<-aggregate(freq~BIN, data=count3, sum)
names(count4)<-c("BIN", "total")
count5<-merge(count4, count3, by="BIN")
count5$rel_prop<-count5$freq/count5$total
count5<-count5[order(count5$BIN, count5$Category),]

out4<-print(paste0("chr",chr,"_",mac,"_mutation_prop2.txt"))
#write.table(count5, out4, col.names=T, row.names=F, quote=F, sep="\t")

#Plot proportions as stacked bars adding up to 1
suppressMessages(ggplot(count5, aes(x=factor(BIN), y=rel_prop, colour=Category, fill=Category, alpha=0.5))+geom_bar(position="stack", stat="identity")+ scale_x_discrete(breaks=seq(0,xmax,50))+xlab("Bin")+ylab("Proportion")+ggtitle(main_scale_title))
suppressMessages(ggsave(main_scale_out))

#Plot dodged and stacked distributions of 6 main categories
ggplot(chr22, aes(x=POS, colour=Category, fill=Category, group=Category, alpha=0.5))+geom_histogram(binwidth=binw, position="identity")+ggtitle(main_dist_title)
ggsave(main_ident_out)
ggplot(chr22, aes(x=POS, colour=Category, fill=Category, group=Category, alpha=0.5))+geom_histogram(binwidth=binw, position="dodge")+ggtitle(main_dist_title)
ggsave(main_dodge_out)
