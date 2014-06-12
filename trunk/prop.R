#! /usr/bin/env Rscript

##############################################################################
# Process Options/Args
##############################################################################

options(warn=-1)
sink("R.log")

suppressMessages(require(ggplot2))
suppressMessages(require(plyr))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))

args<-commandArgs(TRUE)

chr<-as.character(args[1])
macl<-as.character(args[2])
binw<-as.numeric(args[3])
cpg_flag<-as.character(args[4])
summ<-as.character(args[5])
adj<-as.numeric(args[6])
hot_flag<-as.character(args[7])
imgdir<-args[8]

if (macl=="singletons") mac<-"Singleton"
if (macl=="doubletons") mac<-"Doubleton"

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

##############################################################################
#read in summary file and bins file and update columns
#
#summ<-print(paste0("/net/bipolar/jedidiah/bcftools/summaries/",macl,"/all/chr",chr,".",macl,".summary.txt"))
#chr22<-read.table("/net/bipolar/jedidiah/testpipe/summaries/chr22.summary", header=F, stringsAsFactors=F)
#names(chr22)<-c("CHR", "POS", "REF", "ALT", "DP", "AN", "ANNO")
##############################################################################

chr22<-read.table(summ, header=T, stringsAsFactors=F)
bins<-read.table("bin_out.txt", header=T, stringsAsFactors=F, check.names=F)

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

xmax<-floor(max(chr22$BIN)/100)*100

##############################################################################
# Annotation-specific analyses
##############################################################################

##############################################################################
# Default plots for 6 basic categories across chromosome:
# -Dodged/Stacked distribution (counts)
# -Stacked proportions summing to 1
# -Relative mutation rate per bin
##############################################################################

# Initialize titles/output strings
main_rel_title<-print(paste0("Chr",chr," ",mac, " Relative Mutation Rate"))
main_rel_out<-print(paste0(imgdir,"/chr",chr,"_",mac,"_mutation_prop.png"))

main_scale_title<-print(paste0("Chr",chr," ",mac, " Mutations--scaled proportions"))
main_scale_out<-print(paste0(imgdir,"/chr",chr,"_",mac,"_mutation_prop2.png"))

main_dist_title<-print(paste0("Chr",chr," ",mac, " Distribution by Mutation Type"))
main_ident_out<-print(paste0(imgdir,"/chr",chr,"_",mac,"_dist_ident.png"))
main_dodge_out<-print(paste0(imgdir,"/chr",chr,"_",mac,"_dist_dodge.png"))

#Plot dodged and stacked distributions of 6 main categories
ggplot(chr22, aes(x=POS, colour=Category, fill=Category, group=Category, alpha=0.5))+
	geom_histogram(binwidth=binw, position="identity")+
	ggtitle(main_dist_title)
suppressMessages(ggsave(main_ident_out))

ggplot(chr22, aes(x=POS, colour=Category, fill=Category, group=Category, alpha=0.5))+
	geom_histogram(binwidth=binw, position="dodge")+
	ggtitle(main_dist_title)
suppressMessages(ggsave(main_dodge_out))


# Plot proportions as stacked bars adding up to 1
# out4<-print(paste0("chr",chr,"_",mac,"_mutation_prop2.txt"))
# write.table(count, out4, col.names=T, row.names=F, quote=F, sep="\t")
chr22b<-merge(chr22, bins, by="BIN")
count<-count(chr22b, c("Category", "BIN", "AT", "CG"))
count<-merge(count, aggregate(freq~BIN, data=count, sum), by="BIN")
count$rel_prop<-count$freq.x/count$freq.y
count<-count[order(count$BIN, count$Category),]

ggplot(count, aes(x=factor(BIN), y=rel_prop, colour=Category, fill=Category, alpha=0.5))+
	geom_bar(position="stack", stat="identity")+
	scale_x_discrete(breaks=seq(0,xmax,50))+
	xlab("Bin")+
	ylab("Proportion")+
	ggtitle(main_scale_title)
suppressMessages(ggsave(main_scale_out))

# Plot relative mutation rate per bin for 6 main categories 
countAT<-count[grep("^AT", count$Category),]
countGC<-count[grep("^GC", count$Category),]
countAT$prop<-countAT$freq.x/countAT$AT
countGC$prop<-countGC$freq.x/countGC$CG
count2<-rbind(countAT, countGC)

ggplot(count2, aes(x=factor(BIN), y=prop, colour=Category, fill=Category, group=Category, alpha=0.5))+
	geom_bar(position="dodge", stat="identity")+
	scale_x_discrete(breaks=seq(0,xmax,50))+
	xlab("Bin")+
	ylab("Proportion")+
	ggtitle(main_rel_title)
suppressMessages(ggsave(main_rel_out))

# Kataegis plot (experimental--likely easier to create column in perl script)

# chr22$PREV<-c(chr22$POS[1], chr22$POS[-length(chr22$POS)])
# chr22$DIST_NEXT<-chr22$POS-chr22$PREV
# ggplot(chr22, aes(x=POS, y=DIST_NEXT, colour=Category))+geom_point()+scale_y_log10()
# ggsave(imgdir,"/chr22_kataegis.png", width=10, height=4)

##############################################################################
# Local Sequence Analysis
# -currently only for tri-nucleotide sequences
# -will expand for penta-nucleotide sequences
##############################################################################

if (adj==1) {
	bins2<-read.table("bin_out2.txt", header=T, stringsAsFactors=F)
	
	chr22$Sequence<-paste0(pmin(chr22$SEQ, chr22$ALTSEQ),"(",pmax(chr22$SEQ, chr22$ALTSEQ),")")
	cats<-factor(chr22$Category)

	for (i in 1:6) {
		cat<-levels(cats)[i]
		chr22s<-chr22[chr22$Category==cat,]
		title<-paste0("Chr22: ",cat," Mutations by Local Sequence")
		ggplot(chr22s, aes(x=POS))+
			geom_histogram(binwidth=binw, position="identity")+
			ggtitle(title)+
			facet_wrap(~Sequence)
		out<-paste0(imgdir,"/chr",chr," ",cat," by local seq.png")
		suppressMessages(ggsave(out))
	}

	chr22m<-merge(chr22, bins2, by="SEQ")
	chr22m2<-merge(chr22m, bins2, by.x="ALTSEQ", by.y="SEQ")
	chr22m2$COUNT=chr22m2$COUNT.x+chr22m2$COUNT.y

	aggseq<-count(chr22m2, c("Sequence", "Category", "COUNT"))
	aggseq$rel_prop<-aggseq$freq/aggseq$COUNT

	aggseq_a<-aggseq[grep("^A", aggseq$Category),]
	aggseq_g<-aggseq[grep("^G", aggseq$Category),]

	at_seq_title<-paste0("Chr",chr,": AT Mutations by Local Sequence")
	gc_seq_title<-paste0("Chr",chr,": GC Mutations by Local Sequence")
	at_rel_prop_title<-paste0("Chr",chr,": AT Relative Mutation Rate by Local Sequence")
	gc_rel_prop_title<-paste0("Chr",chr,": GC Relative Mutation Rate by Local Sequence")

	at_seq_out<-paste0(imgdir,"/chr",chr,"_AT_seq.png")
	gc_seq_out<-paste0(imgdir,"/chr",chr,"_GC_seq.png")
	at_rel_out<-paste0(imgdir,"/chr",chr,"_AT_rel.png")
	gc_rel_out<-paste0(imgdir,"/chr",chr,"_GC_rel.png")

	ggplot(aggseq_a, aes(x=Category, y=freq, fill=Sequence))+
		geom_bar(position="dodge", stat="identity")+
		ggtitle(at_seq_title)
	suppressMessages(ggsave(at_seq_out))

	ggplot(aggseq_g, aes(x=Category, y=freq, fill=Sequence))+
		geom_bar(position="dodge", stat="identity")+
		ggtitle(gc_seq_title)
	suppressMessages(ggsave(gc_seq_out))

	ggplot(aggseq_a, aes(x=Category, y=rel_prop, fill=Sequence))+
		geom_bar(position="dodge", stat="identity")+
		ggtitle(at_rel_prop_title)+
		ylab("Relative Mutation Rate")
	suppressMessages(ggsave(at_rel_out))

	ggplot(aggseq_g, aes(x=Category, y=rel_prop, fill=Sequence))+
		geom_bar(position="dodge", stat="identity")+
		ggtitle(gc_rel_prop_title)+
		ylab("Relative Mutation Rate")
	suppressMessages(ggsave(gc_rel_out))
	

	#Output datasets and heatmaps

	chr22r<-chr22[,c('BIN', 'Category', 'Sequence')]
	

	#Counts per bin
	pc<-dcast(chr22r, BIN~Category+Sequence)
	bins_r<-bins[pc$BIN,4:ncol(bins)]
	#bins_r<-bins_r[pc$BIN,]
	write.csv(pc, "tri_counts_100kb.csv", row.names=F)
	
	log.pc<-as.matrix(log(pc[-1]+1,2))
	ybreaks<-names(pc[-1])[seq(8,96,16)]
	ylabs<-unique(substr(names(pc[-1]),1,8))
	xbreaks=seq(50-(min(pc[1]) %% 50),nrow(pc),50)
	xlabs=xbreaks+min(pc[1])
	
	count_heat_title<-paste0("Chr",chr,": Mutation Counts Heatmap")
	count_out<-paste0(imgdir,"/chr",chr,"_count_heatmap.png")
	
	ggplot(melt(log.pc), aes(Var1,Var2,fill=value))+
		geom_raster()+
		scale_fill_gradientn(colours=myPalette(100))+
		theme_bw()+
		theme(panel.border=element_blank(),
			legend.position="none",
			axis.text.x = element_text(angle = 90, hjust = 0.5),
			axis.text.y = element_text(angle = 90, hjust = 0.5))+
		xlab(NULL)+
		ylab(NULL)+
		ggtitle(count_heat_title)+
		scale_y_discrete(breaks=ybreaks, labels=ylabs)+
		scale_x_discrete(breaks=xbreaks, labels=xlabs)
	ggsave(count_out)
	

	#Relative mutation rate per bin
	pc1<-pc

	for (i in 2:ncol(pc))  {
		pc1[,i]<-pc[,i]/bins_r[,substr(names(pc)[i], 10,nchar(names(pc)[i]))]
	}
	
	catnames<-substr(names(pc),1,8)
	
	write.csv(pc1, "tri_rel_mut_rate_100kb.csv", row.names=F)
	

	
	log.pc1<-as.matrix(log(pc1[-1]*10000+1,2))
	
	rel_heat_title<-paste0("Chr",chr,": Relative Mutation Rate Heatmap")
	rel_rate_out<-paste0(imgdir,"/chr",chr,"_rel_rate_heatmap.png")
	
	ggplot(melt(log.pc1), aes(Var1,Var2,fill=value))+
		geom_raster()+
		scale_fill_gradientn(colours=myPalette(100))+
		theme_bw()+
		theme(panel.border=element_blank(),
			legend.position="none",
			axis.text.x = element_text(angle = 90, hjust = 0.5),
			axis.text.y = element_text(angle = 90, hjust = 0.5))+
		xlab(NULL)+
		ylab(NULL)+
		ggtitle(rel_heat_title)+
		scale_y_discrete(breaks=ybreaks, labels=ylabs)+
		scale_x_discrete(breaks=xbreaks, labels=xlabs)
	ggsave(rel_rate_out)
	
}


##############################################################################
# CpG Analysis
# -same plots as main analysis, but for the 3 CpG categories
# -relative mutation rate currently excluded
##############################################################################

if (cpg_flag=="on") {
	#names(chr22a)<-c("POS","PAIR","CPGI")
	cpg_scale_title<-print(paste0("Chr",chr," ",mac, " CpG Mutations--scaled proportions"))
	cpg_scale_out<-print(paste0(imgdir,"/chr",chr,"_",mac,"_cpg_mutation_prop2.png"))

	cpg_dist_title<-print(paste0("Chr",chr," ",mac, " CpG Distribution by Mutation Type"))
	cpg_ident_out<-print(paste0(imgdir,"/chr",chr,"_",mac,"_cpg_dist_ident.png"))
	cpg_dodge_out<-print(paste0(imgdir,"/chr",chr,"_",mac,"_cpg_dist_dodge.png"))
	
	chr22b$CpG[chr22b$PAIR=="CG" | chr22b$PAIR=="GC"]<-1
	chr22b$CpG[chr22b$PAIR!="CG" & chr22b$PAIR!="GC"]<-0

	chr22cpg<-chr22b[(chr22b$CpG==1 & chr22b$CPGI==0),]
	chr22cpg$CAT<-paste(chr22cpg$REF, chr22cpg$ALT, sep="")

	chr22cpg$Category[chr22cpg$CAT=="GA" | chr22cpg$CAT=="CT"]<-"GC to AT"
	chr22cpg$Category[chr22cpg$CAT=="GC" | chr22cpg$CAT=="CG"]<-"GC to CG"
	chr22cpg$Category[chr22cpg$CAT=="GT" | chr22cpg$CAT=="CA"]<-"GC to TA"
	
	#Plot distribution of 3 CpG categories
	ggplot(chr22cpg, aes(x=POS, colour=Category, fill=Category, group=Category, alpha=0.5))+
		geom_histogram(binwidth=binw, position="identity")+
		ggtitle(cpg_dist_title)
	ggsave(cpg_ident_out)
	
	ggplot(chr22cpg, aes(x=POS, colour=Category, fill=Category, group=Category, alpha=0.5))+
		geom_histogram(binwidth=binw, position="dodge")+
		ggtitle(cpg_dist_title)
	ggsave(cpg_dodge_out)

	#CpG--Merge bins + summary files and process
	countcpg<-count(chr22cpg, c("Category", "BIN"))
	countcpg<-merge(countcpg, aggregate(freq~BIN, data=countcpg, sum), by="BIN")
	countcpg$rel_prop<-countcpg$freq.x/countcpg$freq.y
	countcpg<-countcpg[order(countcpg$BIN, countcpg$Category),]
	
	ggplot(countcpg, aes(x=factor(BIN), y=rel_prop, colour=Category, fill=Category, alpha=0.5))+
		geom_bar(position="stack", stat="identity")+
		scale_x_discrete(breaks=seq(0,xmax,50))+
		xlab("Bin")+
		ylab("Proportion")+
		ggtitle(cpg_scale_title)
	suppressMessages(ggsave(cpg_scale_out))
} 

##############################################################################
# Heatmaps
# flag is for distance to nearest hotspot heatmap
# also including code for depth heatmap--needs work
##############################################################################

if (hot_flag=="on") {
	#RECOMBINATION HOTSPOTS
	hotspot_heat_out<-print(paste0(imgdir,"/chr",chr,"_",mac,"_mutation_vs_hotspot_heatmap.png"))

	hotspot_agg<-aggregate(DIST ~ BIN+Category, data=chr22, mean)

	ggplot(hotspot_agg, aes(x=BIN, y=Category, fill=DIST))+
		geom_raster()+
		scale_fill_gradientn(colours=myPalette(4))+
		scale_x_continuous(breaks=seq(0,xmax,50))
	suppressMessages(ggsave(hotspot_heat_out))

	#DEPTH	
	depth_heat_out<-print(paste0(imgdir,"/chr",chr,"_",mac,"_mutation_vs_depth_heatmap.png"))

	aggdata<-aggregate(AVGDP ~ BIN+Category, data=chr22, mean)

	#subset to exclude bins with very high avg. depth
	agg2<-aggdata[aggdata$BIN>250,]

	#Graphics--heatmap of avg. depth per bin for the 6 mutation categories
	ggplot(agg2, aes(x=BIN, y=Category, fill=AVGDP))+
		geom_raster()+
		scale_fill_gradientn(colours=myPalette(4))+
		scale_x_continuous(breaks=seq(0,xmax,50))
	suppressMessages(ggsave(depth_heat_out))
}
