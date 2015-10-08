#! /usr/bin/env Rscript

##############################################################################
# Process Options/Args
# Define color palettes
# Define functions
# Initialize strings for titles and file output
##############################################################################

options(warn=-1)
sink("R.log")

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

nbp<-adj*2+1

# Define palettes
myPaletteCat <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
rb<-c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g<-myPaletteG(6)[1:3]


source("get_functions.R")
source("init_titles.R")

##############################################################################
# Read in summary file and bins file and update columns
##############################################################################

{
	chr22 <- read.table(summ, header=T, stringsAsFactors=F)
	bins <- read.table(bin1, header=T, stringsAsFactors=F, check.names=F)
	
	### Specify input files without using args
	# chr22 <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.expanded.summary", header=T, stringsAsFactors=F)
	# bins <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.bin_out.txt", header=T, stringsAsFactors=F)
	
	### Use full input files with all autosomes
	# chr22 <- read.table("/net/bipolar/jedidiah/mutation/output/full.summary", header=T, stringsAsFactors=F)
	# bins <- read.table("/net/bipolar/jedidiah/mutation/output/full_bin.txt", header=T, stringsAsFactors=F, check.names=F)
	
	# chr22 <- chr22[-grep(",", chr22$ALT),] #<-remove multiallelic sites, if not already filtered

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
		substr(chr22$SEQ,adj+1,adj+1) < substr(chr22$ALTSEQ,adj+1,adj+1),
		paste0(chr22$SEQ,"(",chr22$ALTSEQ,")"),
		paste0(chr22$ALTSEQ,"(",chr22$SEQ,")")
	)

	# get complement of sequence columns in bin file and remove duplicates
	for(i in 6:ncol(bins)){
		names(bins)[i] <- paste0(names(bins)[i], "(", revcomp(names(bins)[i]), ")" )
	}
	
	bins2 <- bins[,names(bins)%in%unique(chr22$Sequence)]
	bins <- cbind(bins[,1:5],bins2)

	xmax <- floor(max(chr22$BIN)/100)*100
}

##############################################################################
# Plot distribution of counts faceted by 6 main categories
##############################################################################
	ggplot(chr22, aes(x=POS, colour=Category, fill=Category))+
		geom_histogram(binwidth=binw, position="identity", alpha=0.5)+
		facet_wrap(~Category, scales="free")+
		scale_fill_manual(values=rb)+
		scale_colour_manual(values=rb)+
		# ggtitle(main_dist_title)+
		xlab("Position")+
		theme_bw()+
		# scale_y_log10()+
		theme(panel.border=element_blank(),
			legend.position="none",
			strip.text.x = element_text(size=14),
			axis.text.x = element_text(size=14, angle = 90, hjust = 0.5),
			axis.text.y = element_text(size=14, angle = 90, hjust = 0.5),
			axis.title.y = element_text(size=16),
			axis.title.x = element_text(size=16))
	suppressMessages(ggsave(main_dist_out, width=15.25, height=8.75))

##############################################################################
# Plot distribution of counts for all singletons together
##############################################################################
	ggplot(chr22, aes(x=POS, colour=myPaletteG(6)[1], fill=myPaletteG(6)[1]))+
		geom_histogram(binwidth=binw, position="identity")+
		ggtitle(main_dist_title2)+
		theme_bw()+
		theme(panel.border=element_blank(),
			legend.position="none",
			axis.text.x = element_text(angle = 90, hjust = 0.5),
			axis.text.y = element_text(angle = 90, hjust = 0.5))
	suppressMessages(ggsave(main_dist_out2, width=12, height=4))

##############################################################################
# Plot stacked barchart showing % of each basic mutation subtype
##############################################################################
	chr22b <- merge(chr22, bins, by="BIN")
	# count <- count(chr22b, c("Category", "BIN", "AT", "CG", "prop_GC", "pct")) ## Use this version if considering %exonic, etc.
	count <- count(chr22b, c("Category", "BIN", "AT", "CG", "prop_GC"))
	rm(chr22b)
	count <- merge(count, aggregate(freq~BIN, data=count, sum), by="BIN", all=TRUE)
	count$rel_prop <- count$freq.x/count$freq.y
	count <- count[order(count$BIN, count$Category),]

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

##############################################################################
#Plot relative mutation rate per bin for 6 main categories
##############################################################################
	countAT <- count[grep("^AT", count$Category),]
	countGC <- count[grep("^GC", count$Category),]
	countNULL <- count[count$prop_GC==0,]
	countAT$prop <- countAT$freq.x/countAT$AT
	countGC$prop <- countGC$freq.x/countGC$CG
	countNULL$prop <- countNULL$AT/countNULL$freq.x
	count2 <- rbind(countAT, countGC, countNULL)
	count2[is.na(count2)] <- "AT_CG"
	count2$prop[count2$prop>0.5] <- 0

	p1 <- ggplot(count2, aes(x=BIN, y=prop, colour=Category, fill=Category))+
			geom_bar(position="identity", stat="identity", alpha=0.5)+
			facet_wrap(~Category, scales="free")+
			scale_fill_manual(values=rb)+
			scale_colour_manual(values=rb)+
			#ggtitle(main_rel_title)+
			theme_bw()+
			theme(panel.border=element_blank(),
				legend.position="none",
				axis.text.x = element_text(angle = 90, hjust = 0.5))
	
	mainplot <- p1+
				coord_cartesian(ylim = c(0, 0.03))+
				xlab("Bin")+
				ylab("Relative Mutation Rate")
			
	subplot <- p1+
			   ggtitle("full scale")+
			   theme(strip.background = element_blank(), 
			         strip.text = element_blank(), 
					 axis.text.x=element_blank(), 
					 axis.text.y=element_blank(), 
					 axis.ticks=element_blank(), 
					 axis.title.x=element_blank(), 
					 axis.title.y=element_blank())
	
	vp <- viewport(width = 0.4, height = 0.4, x = 0.8, 0.8)
	vpb_ <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # the larger map
    vpa_ <- viewport(width = 0.25, height = 0.25, x = 0.27, y = 0.8)  # the inset in upper right
	
	png(main_rel_out)
	insetPlot(mainplot, subplot, vpa_)
	dev.off()

##############################################################################
# Local Sequence Analysis
##############################################################################

if (adj>=1) {
	cats <- factor(chr22$Category)
	
	#Counts per bin
	chr22r <- chr22[,c('BIN', 'Category', 'Sequence')]
	pc <- dcast(chr22r, BIN~Category+Sequence)
	pcm <- merge(bins, pc, by="BIN", all=TRUE)
	pcm <- pcm[,names(pc)]
	
	# write.csv(pcm, "tri_counts_100kb.csv", row.names=F)
	
	log.pcm <- as.matrix(log(pcm[-1]+1,2))
	ybreaks <- names(pcm[-1])[seq(8,96,16)]
	ylabs <- unique(substr(names(pcm[-1]),1,5))
	xbreaks <- seq(50-(min(pcm[1]) %% 50),nrow(pcm),50)
	xlabs <- xbreaks+min(pcm[1])

	ggplot(melt(log.pcm), aes(Var1,Var2,fill=value))+
		geom_raster()+
		scale_fill_gradientn(colours=myPalette(100), 
			limits=c(min(log.pcm, na.rm=T),max(log.pcm, na.rm=T)),
			breaks=c(min(log.pcm, na.rm=T),max(log.pcm, na.rm=T)),
			labels=c(min(pcm[2:ncol(pcm)], na.rm=T), max(pcm[2:ncol(pcm)], na.rm=T)))+
		theme_bw()+
		theme(panel.border = element_blank(),
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
	pc1 <- pcm
	bins[bins==0] <- 1
	bins_r <- bins[pcm$BIN,5:ncol(bins)]
	
	for (i in 2:ncol(pc1))  {
		pc1[,i] <- round(pc1[,i]/bins_r[,substr(names(pc1)[i], 7, nchar(names(pc1)[i]))], 4)
	}
	
	catnames <- substr(names(pc),1,8)
	
	# write.csv(pc1, "tri_rel_mut_rate_100kb.csv", row.names=F)
	
	log.pc1 <- as.matrix(log(pc1[-1]*10000+1,2))
	
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
	
	# Get counts of each subsequence across whole chromosome
	bins2 <- melt(bins[,5:ncol(bins)], id="BIN")
	bins2 <- aggregate(data=bins2, value ~ variable, sum)
	names(bins2) <- c("Sequence", "COUNT")
	bins2$Sequence <- sub("[.]", "(", bins2$Sequence)
	bins2$Sequence <- sub("[.]", ")", bins2$Sequence)
	
	bins2$SEQ1 <- substr(bins2$Sequence, 0, adj*2+1)
	bins2$SEQ2 <- substr(bins2$Sequence, (adj*2+1)+2, (adj*2+2)+(adj*2+1))
	bins2$SEQMIN <- pmin(bins2$SEQ1, bins2$SEQ2)
	bins2 <- data.frame(bins2$COUNT, bins2$SEQMIN)
	names(bins2) <- c("COUNT", "SEQMIN")

	# Merge with summary data using min of seq vs altseq
	chr22$SEQMIN <- pmin(chr22$SEQ, chr22$ALTSEQ)
	chr22 <- merge(chr22, bins2, by="SEQMIN")

	aggseq <- count(chr22, c("Sequence", "Category", "CAT", "COUNT", "SEQ"))
	aggseq$rel_prop <- aggseq$freq/aggseq$COUNT
	
	# Test for uniformity of 5bp motifs that share a 3bp motif
	b<-c("A", "C", "G", "T")
	cats<-unique(aggseq$Category)
	
	if(adj==2){
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
	}

	aggseq_a <- aggseq[grep("^A", aggseq$Category),]
	aggseq_g <- aggseq[grep("^G", aggseq$Category),]

	a_seqs <- aggseq_a$Sequence
	map_a <- data.frame(v1=a_seqs)
	map_a$v2 <- substr(map_a$v1,1,adj)
	map_a$v2a <- as.character(lapply(as.vector(map_a$v2), reverse_chars))
	map_a$v2a <- factor(map_a$v2a)
	map_a$v3 <- substr(map_a$v1,adj+2,adj*2+1)
	map_a$v4 <- aggseq_a$rel_prop
	map_a$v5 <- factor(aggseq_a$Category)
	map_a$v6 <- aggseq_a$CAT

	g_seqs <- aggseq_g$Sequence
	map_g <- data.frame(v1=g_seqs)
	map_g$v2 <- substr(map_g$v1,1,adj)
	map_g$v2a <- as.character(lapply(as.vector(map_g$v2), reverse_chars))
	map_g$v2a <- factor(map_g$v2a)
	map_g$v3 <- substr(map_g$v1,adj+2,adj*2+1)
	map_g$v4 <- aggseq_g$rel_prop
	map_g$v5 <- factor(aggseq_g$Category)
	map_g$v6 <- aggseq_g$CAT
	
	levs_a <- as.character(lapply(as.vector(levels(map_a$v2a)), reverse_chars))
	levs_g <- as.character(lapply(as.vector(levels(map_g$v2a)), reverse_chars))
	
	levels(map_a$v5) <- c("A>C", "A>G", "A>T")
	levels(map_g$v5) <- c("C>T", "C>G", "C>A")
	
	map_a1<-aggregate(v4~v1+v2+v2a+v3+v5, map_a, mean)
	map_a1$v4a <- round(map_a1$v4, 3)
	map_a1$v4a[map_a1$v4a<0.001]<-"<0.001"
	
	map_g1<-aggregate(v4~v1+v2+v2a+v3+v5, map_g, mean)
	map_g1$v4a <- round(map_g1$v4, 3)
	map_g1$v4a[map_g1$v4a<0.001]<-"<0.001"
	
	# Define parameters for grouping 3bp motifs
	nbox<-length(unique(map_g$v2a))
	nint<-nbox/4
	xhi <- rep(1:4,4)*nint+0.5
	xlo <- xhi-nint
	yhi <- rep(1:4,each=4)*nint+0.5
	ylo <- yhi-nint
	f <- data.frame(xlo,xhi,ylo,yhi)
	
	##############################################################################
	# Plot relative rate heatmaps
	##############################################################################
	at_heat <- rrheat(map_a1, levs_a, "v5")
	gc_heat <- rrheat(map_g1, levs_g, "v5")
	
	png(at_map_out, width=24, height=24, units="in", res=300)
	multiplot(at_heat, gc_heat, cols=2)
	dev.off()
	

}

##############################################################################
# Run various statistical tests of interest and print to output
##############################################################################
# source("stat_tests.R")

##############################################################################
# CpG Analysis
# -same plots as main analysis, but for the 3 CpG categories
# -relative mutation rate currently excluded
##############################################################################

if (cpg_flag=="on") {
	source("cpg_plots.R")
} 

##############################################################################
# Recombination hotspot analysis
##############################################################################

if (hot_flag=="on") {
	source("recomb_plots.R")
}

##############################################################################
# Annotation-specific analyses
##############################################################################
# source("anno_plots.R")