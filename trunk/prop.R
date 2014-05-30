#! /usr/bin/env Rscript

#####################################################
#Plots mutation types as proportions
#1. Relative to all corresponding nucleotides per bin
#2. Relative to observed total per bin
#####################################################
options(warn=-1)
sink("R.log")

suppressMessages(require(ggplot2))
suppressMessages(require(plyr))

args<-commandArgs(TRUE)

chr<-as.character(args[1])
macl<-as.character(args[2])
binwidth<-as.numeric(args[3])
cpg_flag<-as.character(args[4])
summ<-as.character(args[5])

if (macl=="singletons") mac<-"Singleton"
if (macl=="doubletons") mac<-"Doubleton"

#summ<-print(paste0("/net/bipolar/jedidiah/bcftools/summaries/",macl,"/all/chr",chr,".",macl,".summary.txt"))
title1<-print(paste0("Chr",chr," ",mac, " Relative Mutation Rate"))
title2<-print(paste0("Chr",chr," ",mac, " Mutations--scaled proportions"))
title3<-print(paste0("Chr",chr," ",mac, " CpG Mutations--scaled proportions"))
out1<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_mutation_prop.png"))
out2<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_mutation_prop2.png"))
out3<-print(paste0("/net/bipolar/jedidiah/images/chr",chr,"_",mac,"_cpg_mutation_prop2.png"))

#read in summary file and update columns
chr22<-read.table(summ, header=F, stringsAsFactors=F)
names(chr22)<-c("CHR", "POS", "REF", "ALT", "ANNO")

chr22$CAT<-paste(chr22$REF, chr22$ALT, sep="")
chr22$Category[chr22$CAT=="AC" | chr22$CAT=="TG"]<-"AT to CG"
chr22$Category[chr22$CAT=="AG" | chr22$CAT=="TC"]<-"AT to GC"
chr22$Category[chr22$CAT=="AT" | chr22$CAT=="TA"]<-"AT to TA"
chr22$Category[chr22$CAT=="GA" | chr22$CAT=="CT"]<-"GC to AT"
chr22$Category[chr22$CAT=="GC" | chr22$CAT=="CG"]<-"GC to CG"
chr22$Category[chr22$CAT=="GT" | chr22$CAT=="CA"]<-"GC to TA"

#Read in bins file from perl script
bins<-read.table("bin_out.txt", header=F, stringsAsFactors=F)
names(bins)<-c("AT", "CG")
bins$BIN<-c(1:nrow(bins))

#Read in perl script output containing local sequence and bin number
chr22a<-read.table("cpg_out.txt", header=F, stringsAsFactors=F)

if (cpg_flag=="on") {
	names(chr22a)<-c("POS","PAIR","BIN", "CPGI")
	chr22m<-merge(chr22, chr22a, by="POS")
	chr22m$CpG[chr22m$PAIR=="CG" | chr22m$PAIR=="GC"]<-1
	chr22m$CpG[chr22m$PAIR!="CG" & chr22m$PAIR!="GC"]<-0

	chr22cpg<-chr22m[(chr22m$CpG==1 & chr22m$CPGI==0),]
	chr22cpg$CAT<-paste(chr22cpg$REF, chr22cpg$ALT, sep="")

	chr22cpg$Category[chr22cpg$CAT=="GA" | chr22cpg$CAT=="CT"]<-"GC to AT"
	chr22cpg$Category[chr22cpg$CAT=="GC" | chr22cpg$CAT=="CG"]<-"GC to CG"
	chr22cpg$Category[chr22cpg$CAT=="GT" | chr22cpg$CAT=="CA"]<-"GC to TA"

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
	
	suppressMessages(ggplot(count5cpg, aes(x=factor(BIN), y=rel_prop, colour=Category, fill=Category, alpha=0.5))+geom_bar(position="stack", stat="identity")+ scale_x_discrete(breaks=NULL)+xlab("Bin")+ylab("Proportion")+ggtitle(title3))
	suppressMessages(ggsave(out3))
} else {
	names(chr22a)<-c("POS","PAIR","BIN")
	chr22m<-merge(chr22, chr22a, by="POS")
}	

#All--Merge bins + summary files and process
count<-count(chr22m, c("Category", "BIN", "CHR"))
count2<-merge(count, bins, by="BIN")
countAT<-count2[grep("^AT", count2$Category),]
countGC<-count2[grep("^GC", count2$Category),]
countAT$prop<-countAT$freq/countAT$AT
countGC$prop<-countGC$freq/countGC$CG
count3<-rbind(countAT, countGC)

#Plot proportion relative to all possible sites in each bin
suppressMessages(ggplot(count3, aes(x=factor(BIN), y=prop, colour=Category, fill=Category, group=Category, alpha=0.5))+geom_bar(position="dodge", stat="identity")+ scale_x_discrete(breaks=NULL)+xlab("Bin")+ylab("Proportion")+ggtitle(title1))
suppressMessages(ggsave(out1))

#Obtain proportion of each mutation category relative to total hits in each bin 
count4<-aggregate(freq~BIN, data=count3, sum)
names(count4)<-c("BIN", "total")
count5<-merge(count4, count3, by="BIN")
count5$rel_prop<-count5$freq/count5$total
count5<-count5[order(count5$BIN, count5$Category),]

#Plot proportions as stacked bars adding up to 1
suppressMessages(ggplot(count5, aes(x=factor(BIN), y=rel_prop, colour=Category, fill=Category, alpha=0.5))+geom_bar(position="stack", stat="identity")+ scale_x_discrete(breaks=NULL)+xlab("Bin")+ylab("Proportion")+ggtitle(title2))
suppressMessages(ggsave(out2))
