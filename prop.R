#! /usr/bin/env Rscript

##############################################################################
# Process Options/Args
# Define color palettes
# Define functions
##############################################################################

source("get_opts.R")
source("get_functions.R")

##############################################################################
# Initialize strings for titles and file output
##############################################################################
source("init_titles.R")

##############################################################################
# Read in summary file and bins file and update columns
##############################################################################

{
	chr22 <- read.table(summ, header=T, stringsAsFactors=F)
	bins <- read.table(bin1, header=T, stringsAsFactors=F, check.names=F)

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
}

##############################################################################
# Plot distribution of counts faceted by 6 main categories
##############################################################################
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
		facet_wrap(~Category)+
		scale_fill_manual(values=rb)+
		scale_colour_manual(values=rb)+
		#ggtitle(main_rel_title)+
		theme_bw()+
		theme(panel.border=element_blank(),
			  legend.position="none",
			  axis.text.x = element_text(angle = 90, hjust = 0.5))
	
	mainplot <- p1+
			coord_cartesian(ylim = c(0, 0.015))+
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
# -currently only for tri-nucleotide sequences
##############################################################################

if (adj>=1) {
	cats <- factor(chr22$Category)
	
	#Counts per bin
	chr22r <- chr22[,c('BIN', 'Category', 'Sequence')]
	pc <- dcast(chr22r, BIN~Category+Sequence)
	pcm <- merge(bins, pc, by="BIN", all=TRUE)
	pcm <- pcm[,names(pc)]
	
	write.csv(pcm, "tri_counts_100kb.csv", row.names=F)
	
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
	pc1 <- pcm
	bins[bins==0] <- 1
	bins_r <- bins[pcm$BIN,4:ncol(bins)]
	
	for (i in 2:ncol(pc1))  {
		pc1[,i] <- round(pc1[,i]/bins_r[,substr(names(pc1)[i], 7, nchar(names(pc1)[i]))], 4)
	}
	
	catnames <- substr(names(pc),1,8)
	
	write.csv(pc1, "tri_rel_mut_rate_100kb.csv", row.names=F)
	
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

	# Merge with summary data using min of seq vs altseq
	chr22$SEQMIN <- pmin(chr22$SEQ, chr22$ALTSEQ)
	chr22 <- merge(chr22, bins2, by="SEQMIN")

	aggseq <- count(chr22, c("Sequence", "Category", "CAT", "COUNT", "SEQ"))
	aggseq$rel_prop <- aggseq$freq/aggseq$COUNT

	aggseq_a <- aggseq[grep("^A", aggseq$Category),]
	aggseq_g <- aggseq[grep("^G", aggseq$Category),]

	a_seqs <- aggseq_a$Sequence
	map_a <- data.frame(v1=a_seqs)
	map_a$v2 <- substr(map_a$v1,1,adj)
	map_a$v2a <- as.character(lapply(as.vector(map_a$v2), reverse_chars))
	map_a$v2a <- factor(map_a$v2a)
	map_a$v3 <- substr(map_a$v1,adj+2,adj*2+1)
	map_a$v4 <- aggseq_a$rel_prop
	map_a$v5 <- aggseq_a$Category
	map_a$v6 <- aggseq_a$CAT

	g_seqs <- aggseq_g$Sequence
	map_g <- data.frame(v1=g_seqs)
	map_g$v2 <- substr(map_g$v1,1,adj)
	map_g$v2a <- as.character(lapply(as.vector(map_g$v2), reverse_chars))
	map_g$v2a <- factor(map_g$v2a)
	map_g$v3 <- substr(map_g$v1,adj+2,adj*2+1)
	map_g$v4 <- aggseq_g$rel_prop
	map_g$v5 <- aggseq_g$Category
	map_g$v6 <- aggseq_g$CAT
	
	levs_a <- as.character(lapply(as.vector(levels(map_a$v2a)), reverse_chars))
	levs_g <- as.character(lapply(as.vector(levels(map_g$v2a)), reverse_chars))
	
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
	
	# rrheat(map_a, v5)
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
	
	##############################################################################
	# Plot relative rate heatmaps for uncombined categories
	##############################################################################
	##### Redo data subsets for 12 individual mutations
	a_seqs <- aggseq_a$SEQ
	map_a$v1 <- a_seqs
	map_a$v2 <- substr(map_a$v1,1,adj)
	map_a$v2a <- as.character(lapply(as.vector(map_a$v2), reverse_chars))
	map_a$v2a <- factor(map_a$v2a)
	map_a$v3 <- substr(map_a$v1,adj+2,adj*2+1)
	
	g_seqs <- aggseq_g$SEQ
	map_g$v1 <- g_seqs
	map_g$v2 <- substr(map_g$v1,1,adj)
	map_g$v2a <- as.character(lapply(as.vector(map_g$v2), reverse_chars))
	map_g$v2a <- factor(map_g$v2a)
	map_g$v3 <- substr(map_g$v1,adj+2,adj*2+1)
	
	map_t <- map_a[grep("^T", map_a$v6),]
	map_a <- map_a[grep("^A", map_a$v6),]
	map_c <- map_g[grep("^C", map_g$v6),]
	map_g <- map_g[grep("^G", map_g$v6),]
	
	levs_a <- as.character(lapply(as.vector(levels(map_a$v2a)), reverse_chars))
	levs_g <- as.character(lapply(as.vector(levels(map_g$v2a)), reverse_chars))
	levs_c <- as.character(lapply(as.vector(levels(map_c$v2a)), reverse_chars))
	levs_t <- as.character(lapply(as.vector(levels(map_t$v2a)), reverse_chars))
	
	ggplot()+
		geom_tile(data=map_a, aes(x=v2a, y=v3, fill=log(v4*10000+1,2)))+
		geom_text(data=map_a, aes(x=v2a, y=v3, label=round(v4,3), family="Courier", size=0.1))+
		geom_rect(data=f, size=1.4, colour="grey30", aes(xmin=xlo, xmax=xhi, ymin=ylo, ymax=yhi), fill=NA)+
		scale_fill_gradientn(colours=myPalette((ncol(pc1)-1)/6))+
		xlab("Left flank")+
		ylab("Right flank")+
		theme(legend.position="none")+
		scale_x_discrete(labels=levs_a)+
		facet_wrap(~v6, ncol=1)
	suppressMessages(ggsave(a_map_out, width=12, height=24))
	
	
	ggplot()+
		geom_tile(data=map_g, aes(x=v2a, y=v3, fill=log(v4*10000+1,2)))+
		geom_text(data=map_g, aes(x=v2a, y=v3, label=round(v4,3), family="Courier", size=0.1))+
		geom_rect(data=f, size=1.4, colour="grey30", aes(xmin=xlo, xmax=xhi, ymin=ylo, ymax=yhi), fill=NA)+
		scale_fill_gradientn(colours=myPalette((ncol(pc1)-1)/6))+
		xlab("Left flank")+
		ylab("Right flank")+
		theme(legend.position="none")+
		scale_x_discrete(labels=levs_g)+
		facet_wrap(~v6, ncol=1)
	suppressMessages(ggsave(g_map_out, width=12, height=24))
	
	ggplot()+
		geom_tile(data=map_c, aes(x=v2a, y=v3, fill=log(v4*10000+1,2)))+
		geom_text(data=map_c, aes(x=v2a, y=v3, label=round(v4,3), family="Courier", size=0.1))+
		geom_rect(data=f, size=1.4, colour="grey30", aes(xmin=xlo, xmax=xhi, ymin=ylo, ymax=yhi), fill=NA)+
		scale_fill_gradientn(colours=myPalette((ncol(pc1)-1)/6))+
		xlab("Left flank")+
		ylab("Right flank")+
		theme(legend.position="none")+
		scale_x_discrete(labels=levs_c)+
		facet_wrap(~v6, ncol=1)
	suppressMessages(ggsave(c_map_out, width=12, height=24))
	
	ggplot()+
		geom_tile(data=map_t, aes(x=v2a, y=v3, fill=log(v4*10000+1,2)))+
		geom_text(data=map_t, aes(x=v2a, y=v3, label=round(v4,3), family="Courier", size=0.1))+
		geom_rect(data=f, size=1.4, colour="grey30", aes(xmin=xlo, xmax=xhi, ymin=ylo, ymax=yhi), fill=NA)+
		scale_fill_gradientn(colours=myPalette((ncol(pc1)-1)/6))+
		xlab("Left flank")+
		ylab("Right flank")+
		theme(legend.position="none")+
		scale_x_discrete(labels=levs_t)+
		facet_wrap(~v6, ncol=1)
	suppressMessages(ggsave(t_map_out, width=12, height=24))
	
	##############################################################################
	# Plot count barcharts
	# TODO: Fix error bars?
	##############################################################################
	# y <- apply(pc, 2, std)
	# z <- apply(pc1, 2, std)
	# aggseq$std_count <- y[2:ncol(pc1)]
	# aggseq$std_prop <- z[2:ncol(pc1)]
	
	# limits_count <- aes(ymax=freq+std_count, ymin=freq-std_count)
	# limits_prop <- aes(ymax=rel_prop+std_prop, ymin=rel_prop-std_prop)
	
	
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

	##############################################################################
	# Plot relative rate barcharts
	##############################################################################
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
# Run various statistical tests of interest and print to output
##############################################################################
source("stat_tests.R")

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