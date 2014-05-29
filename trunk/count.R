#! /usr/bin/env Rscript

#########################################################################################
#Plots simple histograms for mutation types, with separate plots for CpG categories
#########################################################################################
options(warn=-1)
sink("R.log")

require(ggplot2)

args<-commandArgs(TRUE)

chr<-as.character(args[1])
macl<-as.character(args[2])
binw<-as.numeric(args[3])
cpg_flag<-as.character(args[4])

if (macl=="singletons") mac<-"Singleton"
if (macl=="doubletons") mac<-"Doubleton"

loc<-print(paste0("/net/bipolar/jedidiah/bcftools/summaries/",macl,"/all/chr",chr,".",macl,".summary.txt"))
title1<-print(paste0("Chr",chr," ",mac, " CpG Distribution by Mutation Type"))
title2<-print(paste0("Chr",chr," ",mac, " Distribution by Mutation Type"))
out1<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_cpg_dist_ident.png"))
out2<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_cpg_dist_dodge.png"))
out3<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_dist_ident.png"))
out4<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_dist_dodge.png"))

#read in summary file and update columns
chr22<-read.table(loc, header=F, stringsAsFactors=F)
names(chr22)<-c("CHR", "POS", "REF", "ALT", "ANNO")

chr22$CAT<-paste(chr22$REF, chr22$ALT, sep="")
chr22$Category[chr22$CAT=="AC" | chr22$CAT=="TG"]<-"AT to CG"
chr22$Category[chr22$CAT=="AG" | chr22$CAT=="TC"]<-"AT to GC"
chr22$Category[chr22$CAT=="AT" | chr22$CAT=="TA"]<-"AT to TA"
chr22$Category[chr22$CAT=="GA" | chr22$CAT=="CT"]<-"GC to AT"
chr22$Category[chr22$CAT=="GC" | chr22$CAT=="CG"]<-"GC to CG"
chr22$Category[chr22$CAT=="GT" | chr22$CAT=="CA"]<-"GC to TA"

#Read in perl script output containing local sequence and bin number
chr22a<-read.table("cpg_out.txt", header=F, stringsAsFactors=F)

if (cpg_flag=="on") {
	names(chr22a)<-c("POS","PAIR","BIN","CPGI")
	chr22m<-merge(chr22, chr22a, by="POS")
	chr22m$CpG[chr22m$PAIR=="CG" | chr22m$PAIR=="GC"]<-1
	chr22m$CpG[chr22m$PAIR!="CG" & chr22m$PAIR!="GC"]<-0

	chr22cpg<-chr22m[(chr22m$CpG==1 & chr22m$CPGI==0),]
	chr22cpg$CAT<-paste(chr22cpg$REF, chr22cpg$ALT, sep="")

	chr22cpg$Category[chr22cpg$CAT=="GA" | chr22cpg$CAT=="CT"]<-"GC to AT"
	chr22cpg$Category[chr22cpg$CAT=="GC" | chr22cpg$CAT=="CG"]<-"GC to CG"
	chr22cpg$Category[chr22cpg$CAT=="GT" | chr22cpg$CAT=="CA"]<-"GC to TA"

	#Plot distribution of 3 CpG categories
	ggplot(chr22cpg, aes(x=POS, colour=Category, fill=Category, group=Category, alpha=0.5))+geom_histogram(binwidth=binw, position="identity")+ggtitle(title1)
	ggsave(out1)
	ggplot(chr22cpg, aes(x=POS, colour=Category, fill=Category, group=Category, alpha=0.5))+geom_histogram(binwidth=binw, position="dodge")+ggtitle(title1)
	ggsave(out2)
} else {
	names(chr22a)<-c("POS","PAIR","BIN")
	chr22m<-merge(chr22, chr22a, by="POS")
}

#Plot distribution of 6 main categories
ggplot(chr22, aes(x=POS, colour=Category, fill=Category, group=Category, alpha=0.5))+geom_histogram(binwidth=binw, position="identity")+ggtitle(title2)
ggsave(out3)
ggplot(chr22, aes(x=POS, colour=Category, fill=Category, group=Category, alpha=0.5))+geom_histogram(binwidth=binw, position="dodge")+ggtitle(title2)
ggsave(out4)
