#! /usr/bin/env Rscript

##############################################################################
# Process Options/Args
# Define color palettes
##############################################################################

options(warn=-1)
sink("R.log")

{
	# Load packages
	suppressMessages(require(ggplot2))
	suppressMessages(require(plyr))
	suppressMessages(require(reshape2))
	suppressMessages(require(RColorBrewer))
	suppressMessages(require(grid))

	args<-commandArgs(TRUE)

	# Read args
	chr<-as.character(args[1])
	macl<-as.character(args[2])
	binw<-as.numeric(args[3])
	cpg_flag<-as.character(args[4])
	summ<-as.character(args[5])
	adj<-as.numeric(args[6])
	hot_flag<-as.character(args[7])
	imgdir<-args[8]
	bin1<-as.character(args[9])
	# bin2<-as.character(args[10])

	if (macl=="singletons") mac<-"Singleton"
	if (macl=="doubletons") mac<-"Doubleton"

	# Define palettes
	myPaletteCat <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
	myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
	myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
	myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
	myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
	rb<-c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
	g<-myPaletteG(6)[1:3]
}
##############################################################################
# Read in summary file and bins file and update columns
#
#summ<-print(paste0("/net/bipolar/jedidiah/bcftools/summaries/",macl,"/all/chr",chr,".",macl,".summary.txt"))
#chr22<-read.table("expanded.summary", header=T, stringsAsFactors=F)
#names(chr22)<-c("CHR", "POS", "REF", "ALT", "DP", "AN", "ANNO")
##############################################################################

{
	chr22<-read.table(summ, header=T, stringsAsFactors=F)
	bins<-read.table(bin1, header=T, stringsAsFactors=F, check.names=F)

	chr22$BIN<-ceiling(chr22$POS/binw)
	chr22$CAT<-paste(chr22$REF, chr22$ALT, sep="")

	# chr22<-chr22[ which(chr22$BIN<260 | chr22$BIN>300),]

	chr22$Category[chr22$CAT=="AC" | chr22$CAT=="TG"]<-"AT_CG"
	chr22$Category[chr22$CAT=="AG" | chr22$CAT=="TC"]<-"AT_GC"
	chr22$Category[chr22$CAT=="AT" | chr22$CAT=="TA"]<-"AT_TA"
	chr22$Category[chr22$CAT=="GA" | chr22$CAT=="CT"]<-"GC_AT"
	chr22$Category[chr22$CAT=="GC" | chr22$CAT=="CG"]<-"GC_CG"
	chr22$Category[chr22$CAT=="GT" | chr22$CAT=="CA"]<-"GC_TA"
	
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
		names(bins)[i]<-paste0(names(bins)[i], "(", revcomp(names(bins)[i]), ")" )
	}
	
	bins2<-bins[,names(bins)%in%unique(chr22$Sequence)]
	bins<-cbind(bins[,1:4],bins2)

	xmax<-floor(max(chr22$BIN)/100)*100

################
# spec$BIN<-ceiling(spec$POS/binw)
# spec$CAT<-paste(spec$REF, spec$ALT, sep="")

# spec$Category[spec$CAT=="AC" | spec$CAT=="TG"]<-"AT_CG"
# spec$Category[spec$CAT=="AG" | spec$CAT=="TC"]<-"AT_GC"
# spec$Category[spec$CAT=="AT" | spec$CAT=="TA"]<-"AT_TA"
# spec$Category[spec$CAT=="GA" | spec$CAT=="CT"]<-"GC_AT"
# spec$Category[spec$CAT=="GC" | spec$CAT=="CG"]<-"GC_CG"
# spec$Category[spec$CAT=="GT" | spec$CAT=="CA"]<-"GC_TA"
}



##############################################################################
# Annotation-specific analyses
##############################################################################



##############################################################################
# Default plots for 6 basic categories across chromosome:
# -Dodged/Stacked distribution (counts)
# -Stacked proportions summing to 1
# -Relative mutation rate per bin
##############################################################################

############# INITIALIZE TITLES/OUTPUT STRINGS
{
	main_rel_title<-paste0("Chr",chr," ",mac, " Relative Mutation Rate")
	main_rel_out<-paste0(imgdir,"/chr",chr,"_",mac,"_mutation_prop.png")

	main_scale_title<-paste0("Chr",chr," ",mac, " Mutations--scaled proportions")
	main_scale_out<-paste0(imgdir,"/chr",chr,"_",mac,"_mutation_prop2.png")

	main_dist_title<-paste0("Chr",chr," ",mac, " Distribution by Mutation Type")
	main_dist_out<-paste0(imgdir,"/chr",chr,"_",mac,"_dist_count.png")
	
	main_dist_title2<-paste0("Distribution of Chr",chr," singletons (low exon density bins)")
	main_dist_out2<-paste0(imgdir,"/chr",chr,"_",mac,"_dist_count_all.png")
}

############## PLOT DISTRIBUTION OF COUNTS FOR 6 MAIN CATEGORIES
{
	ggplot(chr22, aes(x=POS, colour=Category, fill=Category))+
		geom_histogram(binwidth=binw, position="identity", alpha=0.5)+
		facet_wrap(~Category)+
		scale_fill_manual(values=rb)+
		scale_colour_manual(values=rb)+
		ggtitle(main_dist_title)+
		theme_bw()+
		# scale_y_log10()+
		theme(panel.border=element_blank(),
			legend.position="none",
			axis.text.x = element_text(angle = 90, hjust = 0.5),
			axis.text.y = element_text(angle = 90, hjust = 0.5))
	suppressMessages(ggsave(main_dist_out))
}

{
	ggplot(chr22, aes(x=POS, colour=myPaletteG(6)[1], fill=myPaletteG(6)[1]))+
		geom_histogram(binwidth=binw, position="identity")+
		ggtitle(main_dist_title2)+
		# scale_fill_gradient(low="blue", high="red")+
		# scale_colour_gradient(low="blue", high="red")+
		theme_bw()+
		theme(panel.border=element_blank(),
			legend.position="none",
			axis.text.x = element_text(angle = 90, hjust = 0.5),
			axis.text.y = element_text(angle = 90, hjust = 0.5))
	suppressMessages(ggsave(main_dist_out2, width=12, height=4))
}

############# PLOT PROPORTIONS AS STACKED BARS ADDING UP TO 1
{
	# out4<-print(paste0("chr",chr,"_",mac,"_mutation_prop2.txt"))
	# write.table(count, out4, col.names=T, row.names=F, quote=F, sep="\t")
	# chr22b<-merge(chr22, bins, by="BIN", all=TRUE)
	chr22b<-merge(chr22, bins, by="BIN")
	# count<-count(chr22b, c("Category", "BIN", "AT", "CG", "prop_GC", "pct")) ## Use this version if considering %exonic, etc.
	count<-count(chr22b, c("Category", "BIN", "AT", "CG", "prop_GC"))
	rm(chr22b)
	count<-merge(count, aggregate(freq~BIN, data=count, sum), by="BIN", all=TRUE)
	count$rel_prop<-count$freq.x/count$freq.y
	count<-count[order(count$BIN, count$Category),]
	
	#plot dist across chr; colored by pct exonic
			# scale_fill_manual(values=rb)+
		# scale_colour_manual(values=rb)+
		# scale_x_discrete(breaks=seq(0,xmax,50))+
	# ggplot(count, aes(x=BIN, y=freq.y, colour=pct, fill=pct))+
		# geom_bar(position="identity", stat="identity", alpha=0.5)+
		# scale_fill_gradient("%exonic", low="blue", high="green")+
		# scale_colour_gradient(low="blue", high="green", guide=FALSE)+
		# xlab("Bin")+
		# ylab("Count")+
		# ggtitle(main_dist_title2)+
		# theme_bw()+
		# theme(panel.border=element_blank(),
			# axis.text.x = element_text(angle = 90, hjust = 0.5),
			# axis.text.y = element_text(angle = 90, hjust = 0.5))
	# suppressMessages(ggsave(main_dist_out2, width=12, height=4))

	ggplot(count, aes(x=factor(BIN), y=rel_prop, colour=Category, fill=Category))+
		geom_bar(position="stack", stat="identity", alpha=0.5)+
		# scale_x_discrete(breaks=seq(0,xmax,50))+
		scale_fill_manual(values=rb)+
		scale_colour_manual(values=rb)+
		xlab("Bin")+
		ylab("Proportion")+
		ggtitle(main_scale_title)+
		theme_bw()+
		theme(panel.border=element_blank(),
			axis.text.x = element_text(angle = 90, hjust = 0.5),
			axis.text.y = element_text(angle = 90, hjust = 0.5))
	suppressMessages(ggsave(main_scale_out))


	countAT<-count[grep("^AT", count$Category),]
	countGC<-count[grep("^GC", count$Category),]
	countNULL<-count[count$prop_GC==0,]
	countAT$prop<-countAT$freq.x/countAT$AT
	countGC$prop<-countGC$freq.x/countGC$CG
	countNULL$prop<-countNULL$AT/countNULL$freq.x
	count2<-rbind(countAT, countGC, countNULL)
	count2[is.na(count2)]<-"AT_CG"
	count2$prop[count2$prop>0.5]<-0
}

############# AVERAGE RELATIVE RATES ACROSS WHOLE CHROMOSOME
{
	AT_s<-aggregate(freq.x~Category, data=countAT, sum)
	AT_t<-aggregate(AT~Category, data=countAT, sum)
	AT_s$prop<-AT_s$freq.x/AT_t$AT

	GC_s<-aggregate(freq.x~Category, data=countGC, sum)
	GC_t<-aggregate(CG~Category, data=countGC, sum)
	GC_s$prop<-GC_s$freq.x/GC_t$CG

	cwa<-rbind(AT_s, GC_s)
	names(cwa)<-c("variable", "freq.x", "avg")
}
	
############# PLOT RELATIVE MUTATION RATE PER BIN FOR 6 MAIN CATEGORIES 
{

	p1<-ggplot(count2, aes(x=BIN, y=prop, colour=Category, fill=Category))+
		geom_bar(position="identity", stat="identity", alpha=0.5)+
		facet_wrap(~Category)+
		# scale_x_discrete(breaks=seq(0,xmax,50))+
		scale_fill_manual(values=rb)+
		scale_colour_manual(values=rb)+
		#ggtitle(main_rel_title)+
		theme_bw()+
		theme(panel.border=element_blank(),
			legend.position="none",
			axis.text.x = element_text(angle = 90, hjust = 0.5))
	
	mainplot<-p1+
			coord_cartesian(ylim = c(0, 0.015))+
			xlab("Bin")+
			ylab("Relative Mutation Rate")
			
	subplot<-p1+ggtitle("full scale")+theme(strip.background = element_blank(), strip.text = element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
	
	vp <- viewport(width = 0.4, height = 0.4, x = 0.8, 0.8)

	vpb_ <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # the larger map
    vpa_ <- viewport(width = 0.25, height = 0.25, x = 0.27, y = 0.8)  # the inset in upper right
	
	# print(mainplot, vp=vpb_)
	# print(p1, vp=vpa_)
	
	png(main_rel_out)
	
	full <- function() {
		 print(mainplot)
		 theme_set(theme_bw(base_size = 4))
		 print(subplot, vp = vpa_)
		 theme_set(theme_bw())
	 }

	 full()
	 dev.off()
	 
    # ggsave("test_window.png")
	# suppressMessages(ggsave(main_rel_out))
}

############# GC CONTENT HEATMAP	
{
	gc_heat_title<-paste0("Chr",chr," ",mac," GC Content Heatmap")
	gc_heat_out<-paste0(imgdir,"/chr",chr,"_",mac,"_mutation_vs_gc_heatmap.png")

	aggdata<-aggregate(GC ~ BIN+Category, data=chr22, mean)

	ggplot(aggdata, aes(x=BIN, y=Category, fill=GC))+
		geom_raster()+
		scale_fill_gradientn(colours=myPalette(4))+
		scale_x_continuous(breaks=seq(0,xmax,50))+
		scale_y_discrete(breaks=NULL)+
		ylab(NULL)+
		ggtitle(gc_heat_title)+
		theme_bw()+
		theme(panel.border=element_blank(),
			axis.text.x = element_text(angle = 90, hjust = 0.5),
			axis.text.y = element_text(angle = 90, hjust = 0.5))
	suppressMessages(ggsave(gc_heat_out))
}


# Kataegis plot (experimental--not informative yet)

# chr22$PREV<-c(chr22$POS[1], chr22$POS[-length(chr22$POS)])
# chr22$DIST_NEXT<-chr22$POS-chr22$PREV
# ggplot(chr22, aes(x=POS, y=DIST_NEXT, colour=Category))+geom_point()+scale_y_log10()
# ggsave(imgdir,"/chr22_kataegis.png", width=10, height=4)

##############################################################################
# Local Sequence Analysis
# -currently only for tri-nucleotide sequences
##############################################################################

if (adj>=1) {
	cats<-factor(chr22$Category)
	
	#Counts per bin
	chr22r<-chr22[,c('BIN', 'Category', 'Sequence')]
	pc<-dcast(chr22r, BIN~Category+Sequence)
	pcm<-merge(bins, pc, by="BIN", all=TRUE)
	pcm<-pcm[,names(pc)]
	
	write.csv(pcm, "tri_counts_100kb.csv", row.names=F)
	
	log.pcm<-as.matrix(log(pcm[-1]+1,2))
	ybreaks<-names(pcm[-1])[seq(8,96,16)]
	ylabs<-unique(substr(names(pcm[-1]),1,5))
	xbreaks=seq(50-(min(pcm[1]) %% 50),nrow(pcm),50)
	xlabs=xbreaks+min(pcm[1])
	
	count_heat_title<-paste0("Chr",chr," ",mac, " Mutation Counts Heatmap")
	count_out<-paste0(imgdir,"/chr",chr,"_",mac,"_count_heatmap.png")
	
	ggplot(melt(log.pcm), aes(Var1,Var2,fill=value))+
		geom_raster()+
		scale_fill_gradientn(colours=myPalette(100), 
			limits=c(min(log.pcm, na.rm=T),max(log.pcm, na.rm=T)),
			breaks=c(min(log.pcm, na.rm=T),max(log.pcm, na.rm=T)),
			labels=c(min(pcm[2:ncol(pcm)], na.rm=T), max(pcm[2:ncol(pcm)], na.rm=T)))+
		theme_bw()+
		theme(panel.border=element_blank(),
			axis.text.x = element_text(angle = 90, hjust = 0.5),
			axis.text.y = element_text(angle = 90, hjust = 0.5))+
		xlab(NULL)+
		ylab(NULL)+
		ggtitle(count_heat_title)+
		guides(fill = guide_colorbar(title="Count", title.position = "bottom"))+
		scale_y_discrete(breaks=ybreaks, labels=ylabs)+
		scale_x_discrete(breaks=xbreaks, labels=xlabs)
	suppressMessages(ggsave(count_out))
	
	#Relative mutation rate per bin
	pc1<-pcm
	bins[bins==0]<-1
	bins_r<-bins[pcm$BIN,4:ncol(bins)]
	
	for (i in 2:ncol(pc1))  {
		pc1[,i]<-round(pc1[,i]/bins_r[,substr(names(pc1)[i], 7, nchar(names(pc1)[i]))], 4)
	}
	
	catnames<-substr(names(pc),1,8)
	
	write.csv(pc1, "tri_rel_mut_rate_100kb.csv", row.names=F)
	
	log.pc1<-as.matrix(log(pc1[-1]*10000+1,2))
	
	rel_heat_title<-paste0("Chr",chr,": Relative Mutation Rate Heatmap")
	rel_rate_out<-paste0(imgdir,"/chr",chr,"_",mac,"_rel_rate_heatmap.png")
	
	ggplot(melt(log.pc1), aes(Var1,Var2,fill=value))+
		geom_raster()+
		scale_fill_gradientn(colours=myPalette(100),
			limits=c(min(log.pc1, na.rm=T),max(log.pc1, na.rm=T)),
			breaks=c(min(log.pc1, na.rm=T),max(log.pc1, na.rm=T)),
			labels=c(min(pc1[2:ncol(pc1)], na.rm=T), max(pc1[2:ncol(pc1)], na.rm=T)))+
		theme_bw()+
		theme(panel.border=element_blank(),
			axis.text.x = element_text(angle = 90, hjust = 0.5),
			axis.text.y = element_text(angle = 90, hjust = 0.5))+
		xlab(NULL)+
		ylab(NULL)+
		ggtitle(rel_heat_title)+
		guides(fill = guide_colorbar(title="Relative Mutation Rate", title.position = "bottom"))+
		scale_y_discrete(breaks=ybreaks, labels=ylabs)+
		scale_x_discrete(breaks=xbreaks, labels=xlabs)
	ggsave(rel_rate_out)

	############# PANELED BARCHARTS (same info as heatmap, may want to delete)
	# for (i in 1:6) {
		# cat<-levels(cats)[i]
		# chr22s<-chr22[chr22$Category==cat,]
		# title<-paste0("Chr22: ",cat," Mutations by Local Sequence")
		# ggplot(chr22s, aes(x=POS))+
			# geom_histogram(binwidth=binw, position="identity")+
			# ggtitle(title)+
			# theme_bw()+
			# facet_wrap(~Sequence)
		# out<-paste0(imgdir,"/chr",chr,"_",mac," ",cat," by local seq.png")
		# suppressMessages(ggsave(out))
	# }
	
	# Get counts of each subsequence across whole chromosome
	bins2<-melt(bins[,4:((4^(adj*2+1))/2+4)], id="BIN")
	bins2<-aggregate(data=bins2, value ~ variable, sum)
	names(bins2)<-c("Sequence", "COUNT")
	bins2$Sequence<-sub("[.]", "(", bins2$Sequence)
	bins2$Sequence<-sub("[.]", ")", bins2$Sequence)
	
	bins2$SEQ1<-substr(bins2$Sequence, 0, adj*2+1)
	bins2$SEQ2<-substr(bins2$Sequence, (adj*2+1)+2, (adj*2+2)+(adj*2+1))
	bins2$SEQMIN<-pmin(bins2$SEQ1, bins2$SEQ2)
	bins2<-data.frame(bins2$COUNT, bins2$SEQMIN)
	names(bins2)<-c("COUNT", "SEQMIN")

	# Merge with summary data using min of seq vs altseq
	chr22$SEQMIN<-pmin(chr22$SEQ, chr22$ALTSEQ)
	chr22<-merge(chr22, bins2, by="SEQMIN")
	
	# chr22m<-merge(chr22, bins2, by="SEQ")
	# chr22m2<-merge(chr22m, bins2, by.x="ALTSEQ", by.y="SEQ")
	# chr22m2$COUNT=chr22m2$COUNT.x+chr22m2$COUNT.y

	aggseq<-count(chr22, c("Sequence", "Category", "COUNT"))
	aggseq$rel_prop<-aggseq$freq/aggseq$COUNT
	
	std <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
	y<-apply(pc, 2, std)
	z<-apply(pc1, 2, std)
	# aggseq$std_count<-y[2:ncol(pc1)]
	# aggseq$std_prop<-z[2:ncol(pc1)]
	
	# limits_count<-aes(ymax=freq+std_count, ymin=freq-std_count)
	# limits_prop<-aes(ymax=rel_prop+std_prop, ymin=rel_prop-std_prop)

	aggseq_a<-aggseq[grep("^A", aggseq$Category),]
	aggseq_g<-aggseq[grep("^G", aggseq$Category),]
	
	############# TEST PROPORTIONS--GC
	print("Test proportions--GC")
	aggseq_cg<-aggseq_g[grep("G$", aggseq_g$Category),]
	print(prop.test(aggseq_cg$freq, aggseq_cg$COUNT))

	aggseq_ta<-aggseq_g[grep("A$", aggseq_g$Category),]
	print(prop.test(aggseq_ta$freq, aggseq_ta$COUNT))
	
	aggseq_at<-aggseq_g[grep("T$", aggseq_g$Category),]
	print(prop.test(aggseq_at$freq, aggseq_at$COUNT))
	
	############# TEST PROPORTIONS--AT
	print("Test proportions--AT")
	aggseq_cg<-aggseq_a[grep("G$", aggseq_a$Category),]
	print(prop.test(aggseq_cg$freq, aggseq_cg$COUNT))

	aggseq_gc<-aggseq_a[grep("C$", aggseq_a$Category),]
	print(prop.test(aggseq_gc$freq, aggseq_gc$COUNT))
	
	aggseq_ta<-aggseq_a[grep("A$", aggseq_a$Category),]
	print(prop.test(aggseq_ta$freq, aggseq_ta$COUNT))
	
	############# INITIALIZE PLOT TITLES/OUTPUT
	at_seq_title<-paste0("Chr",chr,": AT Mutations by Local Sequence")
	gc_seq_title<-paste0("Chr",chr,": GC Mutations by Local Sequence")
	at_rel_prop_title<-paste0("Chr",chr,": AT Relative Mutation Rate by Local Sequence")
	gc_rel_prop_title<-paste0("Chr",chr,": GC Relative Mutation Rate by Local Sequence")

	at_seq_out<-paste0(imgdir,"/chr",chr,"_",mac,"_AT_seq.png")
	gc_seq_out<-paste0(imgdir,"/chr",chr,"_",mac,"_GC_seq.png")
	at_rel_out<-paste0(imgdir,"/chr",chr,"_",mac,"_AT_rel.png")
	gc_rel_out<-paste0(imgdir,"/chr",chr,"_",mac,"_GC_rel.png")
	
	# Function to reverse sequence--used to correctly plot right flank in subsequence heatmaps
	reverse_chars <- function(string){
		string_split = strsplit(as.character(string), split = "")
		reversed_split = string_split[[1]][nchar(string):1]
		paste(reversed_split, collapse="")
	}
	
	
	a_seqs<-aggseq_a$Sequence
	map_a<-data.frame(v1=a_seqs)
	map_a$v2<-substr(map_a$v1,1,adj)
	map_a$v2a<-as.character(lapply(as.vector(map_a$v2), reverse_chars))
	map_a$v2a<-factor(map_a$v2a)
	map_a$v3<-substr(map_a$v1,adj+2,adj*2+1)
	map_a$v4<-aggseq_a$rel_prop
	map_a$v5<-aggseq_a$Category

	g_seqs<-aggseq_g$Sequence
	map_g<-data.frame(v1=g_seqs)
	map_g$v2<-substr(map_g$v1,1,adj)
	map_g$v2a<-as.character(lapply(as.vector(map_g$v2), reverse_chars))
	map_g$v2a<-factor(map_g$v2a)
	map_g$v3<-substr(map_g$v1,adj+2,adj*2+1)
	map_g$v4<-aggseq_g$rel_prop
	map_g$v5<-aggseq_g$Category
	
	levs_a<-as.character(lapply(as.vector(levels(map_a$v2a)), reverse_chars))
	levs_g<-as.character(lapply(as.vector(levels(map_g$v2a)), reverse_chars))
	
	# map_a<-data.frame(v1=aggseq_a$Sequence, v2=substr(aggseq_a$Sequence, 1, adj), v3=substr(aggseq_a$Sequence, adj+2, ncar(aggseq_a$Sequence))
	# map_g<-data.frame(v1=aggseq_g$Sequence, v2=substr(aggseq_g$Sequence, 1, adj), v3=substr(aggseq_g$Sequence, adj+2, ncar(aggseq_g$Sequence))
	
	at_map_out<-paste0(imgdir,"/chr",chr,"_",mac,"_AT_map.png")
	gc_map_out<-paste0(imgdir,"/chr",chr,"_",mac,"_GC_map.png")
		
	nbox<-length(unique(map_g$v2a))
	nint<-nbox/4
	
	xhi<-rep(1:4,4)*nint+0.5
	xlo<-xhi-nint
	
	yhi<-rep(1:4,each=4)*nint+0.5
	ylo<-yhi-nint
	
	f<-data.frame(xlo,xhi,ylo,yhi)
	
	# Plot relative rate heatmaps
	ggplot()+
		geom_tile(data=map_a, aes(x=v2a, y=v3, fill=log(v4*10000+1,2)))+
		geom_text(data=map_a, aes(x=v2a, y=v3, label=round(v4,3), family="Courier", size=0.1))+
		geom_rect(data=f, size=1.4, colour="grey30", aes(xmin=xlo, xmax=xhi, ymin=ylo, ymax=yhi), fill=NA)+
		scale_fill_gradientn(colours=myPalette((ncol(pc1)-1)/6))+
		xlab("Left flank")+
		ylab("Right flank")+
		theme(legend.position="none")+
		scale_x_discrete(labels=levs_a)+
		facet_wrap(~v5, ncol=1)
	suppressMessages(ggsave(at_map_out, width=12, height=24))
	
	
	ggplot()+
		geom_tile(data=map_g, aes(x=v2a, y=v3, fill=log(v4*10000+1,2)))+
		# geom_rect(data=frames_g, size=1, colour="black", aes(xmin=v2a1-0.5, xmax=v2a1+0.5, ymin=v3a1-0.5, ymax=v3a1+0.5))+
		geom_text(data=map_g, aes(x=v2a, y=v3, label=round(v4,3), family="Courier", size=0.1))+
		geom_rect(data=f, size=1.4, colour="grey30", aes(xmin=xlo, xmax=xhi, ymin=ylo, ymax=yhi), fill=NA)+
		scale_fill_gradientn(colours=myPalette((ncol(pc1)-1)/6))+
		xlab("Left flank")+
		ylab("Right flank")+
		theme(legend.position="none")+
		scale_x_discrete(labels=levs_g)+
		# scale_x_discrete(labels=unique(as.character(lapply(as.vector(map_a$v2a), reverse_chars))))+
		facet_wrap(~v5, ncol=1)
	suppressMessages(ggsave(gc_map_out, width=12, height=24))
	
	# Plot count heatmaps
	ggplot(aggseq_a, aes(x=Category, y=freq, fill=Sequence))+
		geom_bar(position="dodge", stat="identity")+
		# geom_errorbar(limits_count, position="dodge")+
		scale_fill_manual(values=myPalette((ncol(pc1)-1)/6))+
		ggtitle(at_seq_title)+
		theme_bw()+
		theme(panel.border=element_blank(),
			axis.ticks.x=element_blank(),
			legend.position="none",
			panel.grid.major.y=element_blank())
		# guides(col = guide_legend(nrow = 40, byrow = TRUE))
	suppressMessages(ggsave(at_seq_out, width=12, height=7))

	ggplot(aggseq_g, aes(x=Category, y=freq, fill=Sequence))+
		geom_bar(position="dodge", stat="identity")+
		# geom_errorbar(limits_count, position="dodge")+
		scale_fill_manual(values=myPalette((ncol(pc1)-1)/6))+
		ggtitle(gc_seq_title)+
		theme_bw()+
		theme(panel.border=element_blank(),
			axis.ticks.x=element_blank(),
			legend.position="none",
			panel.grid.major.y=element_blank())
		# guides(col = guide_legend(nrow = 40, byrow = TRUE))
	suppressMessages(ggsave(gc_seq_out, width=12, height=7))

	# Plot relative rate barcharts
	ggplot(aggseq_a, aes(x=Sequence, y=rel_prop, fill=Sequence))+
		geom_bar(position="dodge", stat="identity")+
		# geom_errorbar(limits_prop, position="dodge")+
		scale_fill_manual(values=myPalette((ncol(pc1)-1)/6))+
		#ggtitle(at_rel_prop_title)+
		ylab("Relative Mutation Rate")+
		theme_bw()+
		theme(panel.border=element_blank(),
			axis.ticks.x=element_blank(),
			axis.text.x = element_text(angle = 90, hjust = 0.5),
			legend.position="none",
			panel.grid.major.y=element_blank())+
		facet_wrap(~Category, ncol=1)
		# guides(col = guide_legend(nrow = 40, byrow = TRUE))
	suppressMessages(ggsave(at_rel_out, width=24, height=12))
	
	ggplot(aggseq_g, aes(x=Category, y=rel_prop, fill=Sequence))+
		geom_bar(position="dodge", stat="identity")+
		# geom_errorbar(limits_prop, position="dodge")+
		scale_fill_manual(values=myPalette((ncol(pc1)-1)/6))+
		#ggtitle(gc_rel_prop_title)+
		ylab("Relative Mutation Rate")+
		theme_bw()+
		theme(panel.border=element_blank(),
			axis.ticks.x=element_blank(),
			legend.position="none",
			panel.grid.major.y=element_blank())
		# guides(col = guide_legend(nrow = 40, byrow = TRUE))
	suppressMessages(ggsave(gc_rel_out, width=12, height=7))
}

##############################################################################
# CpG Analysis
# -same plots as main analysis, but for the 3 CpG categories
# -relative mutation rate currently excluded
##############################################################################

if (cpg_flag=="on") {
	#names(chr22a)<-c("POS","PAIR","CPGI")
	cpg_scale_title<-paste0("Chr",chr," ",mac, " CpG Mutations--scaled proportions")
	cpg_scale_out<-paste0(imgdir,"/chr",chr,"_",mac,"_cpg_mutation_prop2.png")

	cpg_dist_title<-paste0("Chr",chr," ",mac, " CpG Distribution by Mutation Type")
	cpg_ident_out<-paste0(imgdir,"/chr",chr,"_",mac,"_cpg_dist_ident.png")
	cpg_dodge_out<-paste0(imgdir,"/chr",chr,"_",mac,"_cpg_dist_dodge.png")
	
	chr22b$CpG[chr22b$PAIR=="CG" | chr22b$PAIR=="GC"]<-1
	chr22b$CpG[chr22b$PAIR!="CG" & chr22b$PAIR!="GC"]<-0

	chr22cpg<-chr22b[(chr22b$CpG==1 & chr22b$CPGI==0),]
	chr22cpg$CAT<-paste(chr22cpg$REF, chr22cpg$ALT, sep="")

	chr22cpg$Category[chr22cpg$CAT=="GA" | chr22cpg$CAT=="CT"]<-"GC to AT"
	chr22cpg$Category[chr22cpg$CAT=="GC" | chr22cpg$CAT=="CG"]<-"GC to CG"
	chr22cpg$Category[chr22cpg$CAT=="GT" | chr22cpg$CAT=="CA"]<-"GC to TA"
	
	############# PLOT DODGED AND STACKED DISTRIBUTIONS OF 3 CPG CATEGORIES
	ggplot(chr22cpg, aes(x=POS, colour=Category, fill=Category, group=Category))+
		geom_histogram(binwidth=binw, position="identity", alpha=0.5)+
		ggtitle(cpg_dist_title)+
		scale_colour_manual(values=g)+
		theme_bw()+
		theme(panel.border=element_blank(),
			axis.text.x = element_text(angle = 90, hjust = 0.5),
			axis.text.y = element_text(angle = 90, hjust = 0.5))
	ggsave(cpg_ident_out)
	
	ggplot(chr22cpg, aes(x=POS, colour=Category, fill=Category, group=Category))+
		geom_histogram(binwidth=binw, position="dodge", alpha=0.5)+
		ggtitle(cpg_dist_title)+
		scale_colour_manual(values=g)+
		theme_bw()+
		theme(panel.border=element_blank(),
			axis.text.x = element_text(angle = 90, hjust = 0.5),
			axis.text.y = element_text(angle = 90, hjust = 0.5))
	ggsave(cpg_dodge_out)

	#CpG--Merge bins + summary files and process
	countcpg<-count(chr22cpg, c("Category", "BIN"))
	countcpg<-merge(countcpg, aggregate(freq~BIN, data=countcpg, sum), by="BIN")
	countcpg$rel_prop<-countcpg$freq.x/countcpg$freq.y
	countcpg<-countcpg[order(countcpg$BIN, countcpg$Category),]
	
	############# PLOT RELATIVE MUTATION RATES FOR CPG CATEGORIES
	ggplot(countcpg, aes(x=factor(BIN), y=rel_prop, colour=Category, fill=Category))+
		geom_bar(position="stack", stat="identity", alpha=0.5)+
		scale_x_discrete(breaks=seq(0,xmax,50))+
		xlab("Bin")+
		ylab("Proportion")+
		scale_colour_manual(values=g)+
		scale_fill_manual(values=g)+
		ggtitle(cpg_scale_title)+
		theme_bw()+
		theme(panel.border=element_blank(),
			axis.text.x = element_text(angle = 90, hjust = 0.5),
			axis.text.y = element_text(angle = 90, hjust = 0.5))
	suppressMessages(ggsave(cpg_scale_out))
} 

##############################################################################
# Recombination Hotspot Analysis
# flag is for distance to nearest hotspot heatmap
# also including code for depth heatmap--needs work
##############################################################################

if (hot_flag=="on") {

	recomb_title<-paste0("Chr",chr," ",mac, " Relative Mutation Rates in Recombination Hotspots")
	recomb_out<-paste0(imgdir,"/chr",chr,"_",mac,"_hotspot_rel_mut_rate.png")

	sites<-read.table("new_sites.txt", header=T, stringsAsFactors=F)
	aggsites<-aggregate(Rate~Start, data=sites, mean)
	
	hot<-read.table("hotspot_counts.txt", header=T, stringsAsFactors=F)
	hotm<-melt(hot, id=names(hot)[1:10])
	AT<-hotm[grep("^AT", hotm$variable),]
	GC<-hotm[grep("^GC", hotm$variable),]
	AT$prop<-AT$value/AT$AT
	GC$prop<-GC$value/GC$GC
	hotm2<-rbind(AT,GC)
	
	hotm2<-merge(hotm2, aggsites, by="Start")
	
	#average relative rates across recombination hotspots
	AT_s<-aggregate(value~variable, data=AT, sum)
	AT_t<-aggregate(AT~variable, data=AT, sum)
	AT_s$prop<-AT_s$value/AT_t$AT
	
	GC_s<-aggregate(value~variable, data=GC, sum)
	GC_t<-aggregate(GC~variable, data=GC, sum)
	GC_s$prop<-GC_s$value/GC_t$GC
	
	cwa2<-rbind(AT_s, GC_s)
	#cwa2<-aggregate(prop~variable, data=hotm2, mean)
	
	############# PLOT RELATIVE MUTATION RATES ACROSS HOTSPOTS
	ggplot(hotm2, aes(x=Centre, y=prop, colour=Rate, width=End-Start))+
		geom_bar(position="identity", stat="identity")+
		facet_wrap(~variable)+
		scale_colour_gradientn("Recombination Rate (cM/Mb)", colours=myPalette(4))+
		geom_hline(data=cwa, aes(yintercept=avg))+
		geom_hline(data=cwa2, aes(yintercept=prop), linetype="dashed")+
		theme_bw()+
		theme(panel.border=element_blank())+
		ggtitle(recomb_title)+
		xlab("Position")+
		ylab("Relative Mutation Rate")
	suppressMessages(ggsave(recomb_out, width=12, height=7))

	
	############# PLOT RELATIVE MUTATION RATES FOR 96 SUBTYPES--RECOMBINATION HOTSPOTS ONLY
	bins3<-read.table("bin_out3.txt", header=T, stringsAsFactors=F)
	chr22s<-chr22[chr22$DIST==1,]
	chr22sm<-merge(chr22s, bins3, by="SEQ")
	chr22sm2<-merge(chr22sm, bins3, by.x="ALTSEQ", by.y="SEQ")
	chr22sm2$COUNT=chr22sm2$COUNT.x+chr22sm2$COUNT.y

	aggseq<-count(chr22sm2, c("Sequence", "Category", "COUNT"))
	aggseq$rel_prop<-aggseq$freq/aggseq$COUNT

	aggseq_a<-aggseq[grep("^A", aggseq$Category),]
	aggseq_g<-aggseq[grep("^G", aggseq$Category),]

	at_rel_rec_prop_title<-paste0("Chr",chr,": AT Relative Mutation Rate by Local Sequence\n(recombination hotspots only)")
	gc_rel_rec_prop_title<-paste0("Chr",chr,": GC Relative Mutation Rate by Local Sequence\n(recombination hotspots only)")
	
	at_rel_rec_out<-paste0(imgdir,"/chr",chr,"_",mac,"_AT_rel_rec.png")
	gc_rel_rec_out<-paste0(imgdir,"/chr",chr,"_",mac,"_GC_rel_rec.png")

	ggplot(aggseq_a, aes(x=Category, y=rel_prop, fill=Sequence))+
		geom_bar(position="dodge", stat="identity")+
		#geom_errorbar(limits_prop, position="dodge")+
		scale_fill_manual(values=myPalette(16))+
		ggtitle(at_rel_rec_prop_title)+
		ylab("Relative Mutation Rate")+
		theme_bw()+
		theme(panel.border=element_blank(),
		axis.ticks.x=element_blank(),
		panel.grid.major.y=element_blank())
	ggsave(at_rel_rec_out)

	ggplot(aggseq_g, aes(x=Category, y=rel_prop, fill=Sequence))+
		geom_bar(position="dodge", stat="identity")+
		#geom_errorbar(limits_prop, position="dodge")+
		scale_fill_manual(values=myPalette(16))+
		ggtitle(gc_rel_rec_prop_title)+
		ylab("Relative Mutation Rate")+
		theme_bw()+
		theme(panel.border=element_blank(),
		axis.ticks.x=element_blank(),
		panel.grid.major.y=element_blank())
	ggsave(gc_rel_rec_out)
	
	
	############# PLOT HEATMAP OF RECOMBINATION HOTSPOTS (INDICATOR METHOD)
	# hotspot_heat_out<-paste0(imgdir,"/chr",chr,"_",mac,"_mutation_vs_hotspot_heatmap.png")

	# chr22s<-chr22[chr22$DIST==1,]
	# hotspot_agg<-count(chr22s, c("Category", "BIN"))
	
	# ggplot(hotspot_agg, aes(x=BIN, y=Category, fill=freq))+
		# geom_raster()+
		# scale_fill_gradientn(colours=myPalette(4))+
		# theme_bw()+
		# theme(panel.border=element_blank(),
			# axis.text.x = element_text(angle = 90, hjust = 0.5),
			# axis.text.y = element_text(angle = 90, hjust = 0.5))
	# suppressMessages(ggsave(hotspot_heat_out))
	
	#RECOMBINATION HOTSPOTS--OLD METHOD
	# hotspot_agg<-aggregate(DIST ~ BIN+Category, data=chr22, mean)
	
	# ggplot(hotspot_agg, aes(x=BIN, y=Category, fill=DIST))+
		# geom_raster()+
		# scale_fill_gradientn(colours=myPalette(4))+
		# scale_x_continuous(breaks=seq(0,xmax,50))+
		# theme_bw()
	# suppressMessages(ggsave(hotspot_heat_out))

	#DEPTH	
	
	# chr22$DP<-as.numeric(chr22$DP)
	# chr22$DP[is.na(chr22$DP)]=mean(chr22$DP, na.rm=T)
	# chr22$AVGDP<-chr22$DP/(chr22$AN/2)
	
	# depth_heat_out<-paste0(imgdir,"/chr",chr,"_",mac,"_mutation_vs_depth_heatmap.png")

	# aggdata<-aggregate(AVGDP ~ BIN+Category, data=chr22, mean)

	#subset to exclude bins with very high avg. depth
	# agg2<-aggdata[aggdata$BIN>250,]

	#Graphics--heatmap of avg. depth per bin for the 6 mutation categories
	# ggplot(agg2, aes(x=BIN, y=Category, fill=AVGDP))+
		# geom_raster()+
		# scale_fill_gradientn(colours=myPalette(4))+
		# scale_x_continuous(breaks=seq(0,xmax,50))+
		# theme_bw()+
		# theme(panel.border=element_blank(),
			# axis.text.x = element_text(angle = 90, hjust = 0.5),
			# axis.text.y = element_text(angle = 90, hjust = 0.5))
	# suppressMessages(ggsave(depth_heat_out))
}