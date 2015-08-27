##############################################################################
# Script for running genome-wide models
##############################################################################

##############################################################################
# Process Options/Args
# Define color palettes
# Define functions
##############################################################################
suppressMessages(require(ggplot2))
suppressMessages(require(scales))
suppressMessages(require(plyr))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))
suppressMessages(require(MASS))
suppressMessages(require(speedglm))
suppressMessages(require(grid))
suppressMessages(require(ggbio))
data(hg19IdeogramCyto, package = "biovizBase")

source("/net/bipolar/jedidiah/mutation/smaug-genetics/get_functions.R")

binw <- 100000

myPaletteCat <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
rb <- c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g <- myPaletteG(6)[1:3]

# chr22 <- chr22[-grep(",", chr22$ALT),]

# chr22 <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.expanded.summary", header=T, stringsAsFactors=F)
# bins <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.bin_out.txt", header=T, stringsAsFactors=F)

summ_5bp_100k <- read.table("/net/bipolar/jedidiah/mutation/output/5bp_100k/full.summary", header=T, stringsAsFactors=F)

summ_5bp_100k$BIN <- ceiling(summ_5bp_100k$POS/binw)

bins_5bp_100k <- read.table("/net/bipolar/jedidiah/mutation/output/5bp_100k/full_bin.txt", header=T, stringsAsFactors=F, check.names=F)

# summ_5bp_10k <- read.table("/net/bipolar/jedidiah/mutation/output/5bp_10k/full.summary", header=T, stringsAsFactors=F)
# bins_5bp_10k <- read.table("/net/bipolar/jedidiah/mutation/output/5bp_10k/full_bin.txt", header=T, stringsAsFactors=F, check.names=F)

# summ_3bp_100k <- read.table("/net/bipolar/jedidiah/mutation/output/3bp_100k/full.summary", header=T, stringsAsFactors=F)
# bins_3bp_100k <- read.table("/net/bipolar/jedidiah/mutation/output/3bp_100k/full_bin.txt", header=T, stringsAsFactors=F, check.names=F)

# summ_3bp_10k <- read.table("/net/bipolar/jedidiah/mutation/output/3bp_10k/full.summary", header=T, stringsAsFactors=F)
# bins_3bp_10k <- read.table("/net/bipolar/jedidiah/mutation/output/3bp_10k/full_bin.txt", header=T, stringsAsFactors=F, check.names=F)

##############################################################################
# Add Bin/Category/Sequence columns to summary file, update colnames of bin
# file, and merge overall count for each sequence motif into summary file
##############################################################################
updateData <- function(summfile, binfile, adj){

	nbp <- adj*2+1

	# summfile$BIN <- ceiling(summfile$POS/binw)
	summfile$CAT <- paste(summfile$REF, summfile$ALT, sep="")

	# Manually remove bins near chr20 centromere
	# chr22 <- chr22[ which(chr22$BIN<260 | chr22$BIN>300),]
	summfile$Category[summfile$CAT=="AC" | summfile$CAT=="TG"] <- "AT_CG"
	summfile$Category[summfile$CAT=="AG" | summfile$CAT=="TC"] <- "AT_GC"
	summfile$Category[summfile$CAT=="AT" | summfile$CAT=="TA"] <- "AT_TA"
	summfile$Category[summfile$CAT=="GA" | summfile$CAT=="CT"] <- "GC_AT"
	summfile$Category[summfile$CAT=="GC" | summfile$CAT=="CG"] <- "GC_CG"
	summfile$Category[summfile$CAT=="GT" | summfile$CAT=="CA"] <- "GC_TA"
	

	summfile$Sequence <- ifelse(
		substr(summfile$SEQ,adj+1,adj+1)<substr(summfile$ALTSEQ,adj+1,adj+1),
		paste0(summfile$SEQ,"(",summfile$ALTSEQ,")"),
		paste0(summfile$ALTSEQ,"(",summfile$SEQ,")")
	)
	
	# Second category column to include +3 CpG categories
	summfile$Category2 <- ifelse(substr(summfile$Sequence,adj+1,adj+2)=="CG", 
								paste0("cpg_",summfile$Category), 
								summfile$Category)

	# get complement of sequence columns in bin file and remove duplicates
	for(i in 6:ncol(binfile)){
		names(binfile)[i] <- paste0(names(binfile)[i], "(", revcomp(names(binfile)[i]), ")" )
	}

	bins2 <- binfile[,names(binfile)%in%unique(summfile$Sequence)]
	binfile <- cbind(binfile[,1:5],bins2)
	xmax <- floor(max(summfile$BIN)/100)*100

	bins2 <- melt(binfile[,5:ncol(binfile)], id="BIN")
	bins2 <- aggregate(data=bins2, value ~ variable, sum)
	names(bins2) <- c("Sequence", "COUNT")
	bins2$Sequence <- sub("[.]", "(", bins2$Sequence)
	bins2$Sequence <- sub("[.]", ")", bins2$Sequence)
	bins2$SEQ1 <- substr(bins2$Sequence, 0, adj*2+1)
	bins2$SEQ2 <- substr(bins2$Sequence, (adj*2+1)+2, (adj*2+2)+(adj*2+1))
	bins2$SEQMIN <- pmin(bins2$SEQ1, bins2$SEQ2)
	bins2 <- data.frame(bins2$COUNT, bins2$SEQMIN)
	names(bins2) <- c("COUNT", "SEQMIN")

	summfile$SEQMIN <- pmin(summfile$SEQ, summfile$ALTSEQ)
	summfile <- merge(summfile, bins2, by="SEQMIN")
	
	datalist<- list("summ"=summfile, "bin"=binfile)
	return(datalist)
}

dat_5bp_100k <- updateData(summ_5bp_100k, bins_5bp_100k, 2)
# dat_5bp_2 <- updateData(summ_5bp_100k[summ_5bp_100k$CHR==2 & summ_5bp_100k$BIN>=300 & summ_5bp_100k$BIN<600,], bins_5bp_100k[bins_5bp_100k$CHR=="chr2" & bins_5bp_100k$BIN>=300 & bins_5bp_100k$BIN<600,], 2)
rm(summ_5bp_100k)
rm(bins_5bp_100k)
# gc()
# dat_5bp_10k<-updateData(summ_5bp_10k, bins_5bp_10k, 2)
# dat_3bp_100k<-updateData(summ_3bp_100k, bins_3bp_100k, 1)
# dat_3bp_10k<-updateData(summ_3bp_10k, bins_3bp_10k, 1)

##############################################################################
# Aggregate data, plot motif heatmap, output predicted counts into each bin
##############################################################################
aggData <- function(datfile, adj){
	
	summfile <- datfile$summ
	binfile <- datfile$bin
	nbp <- adj*2+1
	
	# summfile <- dat_5bp_100k$summ
	# binfile <- dat_5bp_100k$bin
	
	# summfile <- dat_5bp_2$summ
	# binfile <- dat_5bp_2$bin
	
	# summfile <-summfile[summfile$CHR==2 & summfile$BIN<=600 & summfile$BIN>=300,]
	# binfile <- binfile[binfile$CHR=="chr2" & binfile$BIN<=600 & binfile$BIN>=300,]
	
	# Plot genome-wide motif heatmaps
	plot_heatmap <- 0
	if(plot_heatmap==1){
	
		aggseq <- count(summfile, c("Sequence", "Category", "CAT", "COUNT", "SEQ"))
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
		
		# Plot relative rate heatmaps
		at_heat <- rrheat(map_a1, f, levs_a, "v5", nbp)
		gc_heat <- rrheat(map_g1, f, levs_g, "v5", nbp)
		
		png("/net/bipolar/jedidiah/mutation/smaug-genetics/gw_map.png", width=18, height=24, units="in", res=300)
		multiplot(at_heat, gc_heat, cols=2)
		dev.off()
		
		# Plot relative rate heatmaps for uncombined categories
		# Redo data subsets for 12 individual mutations
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
		
		map_a$v4a <- round(map_a$v4, 3)
		map_a$v4a[map_a$v4a<0.001]<-"<0.001"
		map_t$v4a <- round(map_t$v4, 3)
		map_t$v4a[map_t$v4a<0.001]<-"<0.001"
		map_c$v4a <- round(map_c$v4, 3)
		map_c$v4a[map_c$v4a<0.001]<-"<0.001"
		map_g$v4a <- round(map_g$v4, 3)
		map_g$v4a[map_g$v4a<0.001]<-"<0.001"
		
		levs_a <- as.character(lapply(as.vector(levels(map_a$v2a)), reverse_chars))
		levs_g <- as.character(lapply(as.vector(levels(map_g$v2a)), reverse_chars))
		levs_c <- as.character(lapply(as.vector(levels(map_c$v2a)), reverse_chars))
		levs_t <- as.character(lapply(as.vector(levels(map_t$v2a)), reverse_chars))
			
		a_heat <- rrheat(map_a, f, levs_a, "v6", nbp)
		t_heat <- rrheat(map_t, f, levs_t, "v6", nbp)
		c_heat <- rrheat(map_c, f, levs_c, "v6", nbp)
		g_heat <- rrheat(map_g, f, levs_g, "v6", nbp)
		
		png("/net/bipolar/jedidiah/mutation/images/gw_map_uncollapsed.png", width=48, height=24, units="in", res=300)
		multiplot(a_heat, t_heat, c_heat, g_heat, cols=4)
		dev.off()
			
		cats<-unique(aggseq$Category)
		seqs<-unique(aggseq$Sequence)
		
		# Test for differences in proportions on opposite strands
		if(adj==2){
			for(i in 1:6){
				for(j in 1:512){
					cat<-cats[i]
					seq<-seqs[j]
					dat<-aggseq[(aggseq$Sequence==seqs[j] & aggseq$Category==cats[i]),]
					if(nrow(dat)==2){
						test<-prop.test(dat$freq, dat$COUNT)
						if(test$p.value<0.05/1536){
							print(cat)
							print(seq)
							print(dat)
							print(test)
							print(test$p.value)
							# print(dat)
						}
					}
				}
			}
		}
	}
	
	# Get dataframe of observed and predicted counts
	{
		aggseq <- count(summfile, c("Sequence", "Category2", "COUNT"))
		aggseq$rel_prop <- aggseq$freq/aggseq$COUNT
		
		ct <- count(summfile,c("CHR","BIN","Category2"))
		ct.ord <- ct[order(ct$Category2),]
		
		binfile$BIN <- paste0(binfile$CHR,".",binfile$BIN)
		
		binsT <- setNames(data.frame(t(binfile[,-c(1:5)])), binfile$BIN)
		binsT$Sequence <- rownames(binsT)

		# merge counts per sequence/category with counts of mutable motifs
		aggseq.m <- merge(aggseq, binsT, by="Sequence")

		# get expected counts for each row
		aggseq.m[,6:ncol(aggseq.m)] <- aggseq.m$rel_prop*aggseq.m[,6:ncol(aggseq.m)]

		# sum all sequence combinations for each category/bin
		atcg <- colSums(aggseq.m[aggseq.m$Category2=="AT_CG",-c(1:5)])
		atgc <- colSums(aggseq.m[aggseq.m$Category2=="AT_GC",-c(1:5)])
		atta <- colSums(aggseq.m[aggseq.m$Category2=="AT_TA",-c(1:5)])
		gcat <- colSums(aggseq.m[aggseq.m$Category2=="GC_AT",-c(1:5)])
		gccg <- colSums(aggseq.m[aggseq.m$Category2=="GC_CG",-c(1:5)])
		gcta <- colSums(aggseq.m[aggseq.m$Category2=="GC_TA",-c(1:5)])
		cpg_gcat <- colSums(aggseq.m[aggseq.m$Category2=="cpg_GC_AT",-c(1:5)])
		cpg_gccg <- colSums(aggseq.m[aggseq.m$Category2=="cpg_GC_CG",-c(1:5)])
		cpg_gcta <- colSums(aggseq.m[aggseq.m$Category2=="cpg_GC_TA",-c(1:5)])

		# vector of expected counts
		exp <- c(atcg, atgc, atta, gcat, gccg, gcta, cpg_gcat, cpg_gccg, cpg_gcta)

		# vector of categories
		catrep <- c(rep("AT_CG",ncol(aggseq.m)-5), rep("AT_GC",ncol(aggseq.m)-5), rep("AT_TA",ncol(aggseq.m)-5), 
				  rep("GC_AT",ncol(aggseq.m)-5), rep("GC_CG",ncol(aggseq.m)-5), rep("GC_TA",ncol(aggseq.m)-5),
				  rep("cpg_GC_AT",ncol(aggseq.m)-5), rep("cpg_GC_CG",ncol(aggseq.m)-5), rep("cpg_GC_TA",ncol(aggseq.m)-5))
				  
		# create data frame with observed, expected, and bin	
		BIN <- as.integer(gsub(".*\\.", "", names(exp)))
		CHR <- as.integer(substring(gsub("\\..*", "", names(exp)), 4))
		
		oe2 <- data.frame(CHR, BIN, Category2=catrep, exp)
		oe2 <- merge(oe2,ct.ord, by=c("CHR","Category2","BIN"))
		names(oe2)[5] <- "obs"
		oe2$res <- paste0(nbp,"bp")

		# get odds ratio for each bin/category
		oe2$odds <- oe2$obs/oe2$exp

		# order by OR and add column of ranks
		# oe2.ord <- oe2[order(oe2$Category2, oe2$odds),]
		# oe2.ord <- ddply(oe2.ord, .(Category2), transform, rk=seq_along(Category2))
		# oe2.ord$res <- paste0(nbp,"bp")

		# return(oe2.ord)
		
		datalist<- list("agg"=aggseq, "oe"=oe2)
		return(datalist)
	}
}

aggV <- aggData(dat_5bp_100k, 2)
aggV2 <- aggData(dat_5bp_2, 2)

agg_5bp_100k <- aggV$oe
rates1 <- aggV$agg
agg_5bp_2 <- aggV2$oe
rates2 <- aggV2$agg

agg_5bp_100k2 <- rbind(agg_5bp_100k[!(agg_5bp_100k$CHR %in% agg_5bp_2$CHR & agg_5bp_100k$BIN %in% agg_5bp_2$BIN),], agg_5bp_2)
# agg_5bp_10k<-aggData(dat_5bp_10k, 2)
# agg_3bp_100k<-aggData(dat_3bp_100k, 1)
# agg_3bp_10k<-aggData(dat_3bp_10k, 1)

agg_5bp_100k <- agg_5bp_100k[order(agg_5bp_100k$BIN),]
agg_5bp_100k$diff <- agg_5bp_100k$obs-agg_5bp_100k$exp
agg_5bp_100k$Category2 <- as.character(agg_5bp_100k$Category2)


##############################################################################
# Test if rates on chr2p are significantly different than genomewide
##############################################################################
mut_cats <- unique(rates1$Category2)

sigmotifs<-data.frame()
for(i in 1:6){
	mut_seqs <- unique(rates1[rates1$Category2==mut_cats[i],]$Sequence)
	for(j in 1:256){

		mut <- c(rates1[rates1$Category2==mut_cats[i] & rates1$Sequence==mut_seqs[j],]$freq, rates2[rates2$Category2==mut_cats[i] & rates2$Sequence==mut_seqs[j],]$freq)
		sites <- c(rates1[rates1$Category2==mut_cats[i] & rates1$Sequence==mut_seqs[j],]$COUNT, rates2[rates2$Category2==mut_cats[i] & rates2$Sequence==mut_seqs[j],]$COUNT)
		
		if(!is.na(mut)){
			z<-prop.test(mut, sites)
			
			if(z$p.value<0.05/1536){
			
				mut_cat<-mut_cats[i]
				mut_seq<-mut_seqs[j]
				
				motifrow<-data.frame(mut_cat, mut_seq, z$p.value)
				
				sigmotifs<-rbind(sigmotifs, motifrow)
			}
		}
	}
}

##############################################################################
# Compare 3bp/5bp motifs at 100kb resolution
##############################################################################
{
	# oe.full <- rbind(agg_3bp_100k, agg_5bp_100k)
	# oe.full2 <- oe.full[(oe.full$odds>0.5 & oe.full$odds<2),]
	# ggplot(oe.full2, aes(x=rk, y=odds, group=res, colour=res))+
		# geom_point()+
		# facet_wrap(~Category, scales="free")
	# ggsave("/net/bipolar/jedidiah/mutation/images/3bp_vs_5bp.png")

	# ggplot(oe.full, aes(factor(res), exp))+
		# geom_boxplot()+
		# facet_wrap(~Category, scales="free")+
		# theme_bw()+
		# theme(strip.text.x=element_text(size=14))
	# ggsave("/net/bipolar/jedidiah/mutation/images/3bp_vs_5bp_box.png", width=13.5, height=8.5)
}

##############################################################################
# Get 1bp predictions
##############################################################################
singlebp<-0
if(singlebp==1){
	aggseq6<-ddply(agg_5bp_100k, .(Category2), summarize, COUNT=sum(obs), freq=sum(freq))
	aggseq6$rel_prop<-aggseq6$freq/aggseq6$obs

	atcg<-aggseq6[1,4]*bins[,1]
	atgc<-aggseq6[2,4]*bins[,1]
	atta<-aggseq6[3,4]*bins[,1]
	gcat<-aggseq6[4,4]*bins[,2]
	gccg<-aggseq6[5,4]*bins[,2]
	gcta<-aggseq6[6,4]*bins[,2]

	exp<-c(atcg, atgc, atta, gcat, gccg, gcta)

	# vector of categories
	catrep<-c(rep("AT_CG",length(atcg)), 
			  rep("AT_GC",length(atgc)), 
			  rep("AT_TA",length(atta)), 
			  rep("GC_AT",length(gcat)), 
			  rep("GC_CG",length(gccg)), 
			  rep("GC_TA",length(gcta)))
			  
	oe1bp<-data.frame(exp, Category=catrep, BIN=rep(1:length(atcg),6))
	oe1bp<-merge(oe1bp,ct.ord, by=c("Category","BIN"))
	names(oe1bp)[4]<-"obs"

	oe1bp$odds<-oe1bp$obs/oe1bp$exp

	# order by OR and add column of ranks
	oe1bp.ord<-oe1bp[order(oe1bp$Category, oe1bp$odds),]
	oe1bp.ord <- ddply(oe1bp.ord, .(Category), transform, rk=seq_along(Category))
	oe1bp.ord$res<-paste0(1,"bp")
}

##############################################################################
# Build matrix of genomic features
##############################################################################
{
	getHistBins <- function(file, mark){

		data <- read.table(file, header=F, stringsAsFactors=F, sep="\t")

		names(data) <- c("CHR", "Start", "End", "mark")
		names(data)[4] <- mark
		
		data$CHR <- as.integer(substring(data$CHR, 4))
		data$BIN <- ceiling(data$End/binw)

		# Old version--parse default output from bedops with "|" as delimiter
		# data$count <- as.integer(gsub(".*\\|", "", data$End))
		# data$End <- as.integer(gsub("\\|.*", "", data$End))

		return(data[,c(1,4,5)])
	}

	# Number of reads per histone mark
	H3K4me1 <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K4me1.100kb.bed", "H3K4me1")
	H3K4me3 <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K4me3.100kb.bed", "H3K4me3")
	H3K9ac <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K9ac.100kb.bed", "H3K9ac")
	H3K9me3 <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K9me3.100kb.bed", "H3K9me3")
	H3K27ac <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K27ac.100kb.bed", "H3K27ac")
	H3K27me3 <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K27me3.100kb.bed", "H3K27me3")
	H3K36me3 <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K36me3.100kb.bed", "H3K36me3")

	# Pct in CpG island
	CPGI <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/cpg_islands_100kb.bed", "CPGI")
	CPGI$CPGI <- CPGI$CPGI/binw

	# Pct exonic 
	EXON <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/coding_exon_100kb.bed", "EXON")
	EXON$EXON <- EXON$EXON/binw

	# Replication timing
	reptime <- read.table("/net/bipolar/jedidiah/mutation/reference_data/lymph_rep_time.txt", header=F, stringsAsFactors=F, sep="\t")
	names(reptime) <- c("CHR", "POS", "TIME")
	reptime$BIN <- ceiling(reptime$POS/binw)
	rtagg <- aggregate(TIME~CHR+BIN, reptime, mean)

	# Recombination rate
	rcrate <- read.table("/net/bipolar/jedidiah/mutation/reference_data/recomb_rate.bed", header=T, stringsAsFactors=F, sep="\t")
	rcrate$CHR <- as.integer(substring(rcrate$CHR, 4))
	rcrate$POS <- (rcrate$START+rcrate$END)/2
	rcrate$BIN <- ceiling(rcrate$POS/binw)
	rcagg <- aggregate(RATE~CHR+BIN, rcrate, mean)

	# GC content
	pctgc <- dat_5bp_100k$bin[,c(1,5,4)]
	pctgc$CHR <- as.integer(substring(pctgc$CHR, 4))

	# Merge all data and sort

	# Old Version--reduce-merge
	# mut_cov <- Reduce(function(x, y) merge(x, y, all=TRUE),
						  # list(H3K4me1, H3K4me3, H3K9ac, H3K9me3, H3K27ac, H3K27me3, H3K36me3, CPGI, EXON))

	# New version--cbind histone marks, cpgi, exon, then merge with rep timing and rc rate
	mut_cov <- cbind(H3K4me1, H3K4me3[,2], H3K9ac[,2], H3K9me3[,2], 
				   H3K27ac[,2], H3K27me3[,2], H3K36me3[,2], CPGI[,2], EXON[,2])
	names(mut_cov) <- c("CHR", "H3K4me1", "BIN", "H3K4me3", "H3K9ac", "H3K9me3", 
					  "H3K27ac", "H3K27me3", "H3K36me3", "CPGI", "EXON")
	mut_cov <- merge(mut_cov, rtagg, by=c("CHR", "BIN"))
	mut_cov <- merge(mut_cov, rcagg, by=c("CHR", "BIN"))
	mut_cov <- merge(mut_cov, pctgc, by=c("CHR", "BIN"))
	
	pc.dat <- prcomp(mut_cov[,3:ncol(mut_cov)], center=T, scale=T)
	scores <- as.data.frame(round(pc.dat$x, 3))
		
	mut_cov2<-cbind(mut_cov[,1:2], scores)
	write.table(mut_cov2, "/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/mut_cov2.txt", sep="\t", col.names=F, row.names=F, quote=F) 
}

##############################################################################
# Compare 5bp predictions vs genomic features model across 100kb bins
# TO-DO:
# -add category of CpG mutations and run model on 9 categories instead of 6
##############################################################################
agg_cov <- merge(agg_5bp_100k[,c(1,2,3,4,5)], mut_cov, by=c("CHR", "BIN"))

# Motifs + genomic features
compare.all <- data.frame()
mut_cats <- unique(agg_5bp_100k$Category2)
for(i in 1:length(mut_cats)){
	cat1 <- mut_cats[i]
	aggcat <- agg_cov[agg_cov$Category==mut_cats[i],]
	pc.dat <- prcomp(aggcat[,6:17], center=T, scale=T)
	scores <- as.data.frame(pc.dat$x)
	
	aggcat <- cbind(aggcat, scores)

	mut_lm_feat <- glm.nb(obs~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12, data=aggcat)
	mut_lm_motif <- glm.nb(obs~exp, data=aggcat)
	mut_lm_full <- glm.nb(obs~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+exp, data=aggcat)
	
	fits_feat <- mut_lm_feat$fitted.values
	names(fits_feat) <- paste0(aggcat$CHR,".",aggcat$BIN)
	
	fits_motif <- mut_lm_motif$fitted.values
	names(fits_motif) <- paste0(aggcat$CHR,".",aggcat$BIN)
	
	fits_full <- mut_lm_full$fitted.values
	names(fits_full) <- paste0(aggcat$CHR,".",aggcat$BIN)
	
	BIN <- as.integer(gsub(".*\\.", "", names(fits_feat)))
	CHR <- as.integer(gsub("\\..*", "", names(fits_feat)))
	
	model_dat_feat <- data.frame(CHR,Category2=cat1, BIN, 
				   exp=fits_feat, 
				   obs=aggcat$obs,
				   stringsAsFactors=F)
				   
	model_dat_feat$res <- "features"
	
	model_dat_motif <- data.frame(CHR,Category2=cat1, BIN, 
				   exp=fits_motif, 
				   obs=aggcat$obs,
				   stringsAsFactors=F)
				   
	model_dat_motif$res <- "5bp"
	
	model_dat_full <- data.frame(CHR,Category2=cat1, BIN, 
				   exp=fits_full, 
				   obs=aggcat$obs,
				   stringsAsFactors=F)
				   
	model_dat_full$res <- "5bp+features"
	
	compare.all <- rbind(compare.all,model_dat_feat,model_dat_motif,model_dat_full)
	
	# z <- summary(mut.lm)
	# print(cat1)
	# print(z)
}


# compare.all<-rbind(agg_5bp_100k[,c(1,2,3,4,5,8)], d, d1) #<- this version uses direct estimates from motif rates; not fair comparison
compare.all$res <- factor(compare.all$res, levels = c("features", "5bp", "5bp+features"))
compare.all$diff <- compare.all$obs-compare.all$exp

compare.all$Category2 <- factor(compare.all$Category2)
levels(compare.all$Category2) <- c("AT>CG", "AT>GC", "AT>TA", "(CpG)GC>AT", "(CpG)GC>CG", "(CpG)GC>TA",  "GC>AT", "GC>CG", "GC>TA")

# Function to compute standard error for correlations
corSE<-function(corval, ct){
	sqrt((1-corval^2)/(ct-2))
}

# Plot scatterplot across chr of model errors (using negbin motif predictions)
p2 <- ggplot(compare.all[compare.all$CHR==2 & compare.all$diff>-150,], aes(x=BIN, y=diff, colour=res))+
	# geom_bar(stat="identity", position="identity")+
	geom_point(alpha=0.1, size=4)+
	scale_colour_manual("Model", values=myPaletteCat(8)[6:8], labels=c("exp"="5bp+features", "obs"="5bp only"))+
	facet_wrap(~Category2, scales="free", ncol=1)+
	# scale_x_discrete(labels=c("exp"="5bp+features", "obs"="5bp only"))+
	ylab("Error")+
	xlab(NULL)+
	theme_bw()+
	theme(axis.title.y=element_text(size=16), 
	      axis.text.y=element_text(size=14),
		  legend.title=element_text(size=16),
		  legend.text=element_text(size=14),
		  axis.text.x=element_blank(), 
		  axis.ticks.x=element_blank())
ggsave("/net/bipolar/jedidiah/mutation/images/hier_diffs_pt5.png", width=9, height=18)

# Plot barcharts comparing obs/exp correlation for different models
mod.corr <- ddply(compare.all, .(Category2, res), summarize, num=length(exp), cor=cor(exp, obs, method="spearman"))
mod.corr$SE <- corSE(mod.corr$cor, mod.corr$num)
limits <- aes(ymax = mod.corr$cor + mod.corr$SE, ymin=mod.corr$cor - mod.corr$SE)
dodge <- position_dodge(width=0.9)

ggplot(mod.corr, aes(x=Category2, y=cor, fill=res))+
	geom_bar(stat="identity", position=dodge)+
	scale_colour_brewer("Predictor", palette="Dark2")+
	scale_fill_brewer("Predictor", palette="Dark2")+
	xlab("Category")+
	ylab("Correlation with observed count")+
	geom_errorbar(limits, position=dodge, width=0.25)+
	theme_bw()+
	theme(legend.title = element_text(size=18),
		  legend.text = element_text(size=16),
		  axis.title.x = element_text(size=20),
		  axis.title.y = element_text(size=20),
		  axis.text.y = element_text(size=16), 
		  axis.text.x = element_text(size=16, angle = 45,  vjust=1, hjust=1.01))
ggsave("/net/bipolar/jedidiah/mutation/images/gw_5bp_vs_mod.png", width=7, height=7)


mut.diff<-merge(agg_5bp_100k, mut_cov, by=c("CHR", "BIN"))

d.diff <- data.frame()
for(i in 1:length(mut_cats)){
	cat1 <- mut_cats[i]
	aggcat <- mut.diff[mut.diff$Category2==cat1,]
	
	pc.dat <- prcomp(aggcat[,9:20], center=T, scale=T)
	scores <- as.data.frame(pc.dat$x)
	
	aggcat <- cbind(aggcat, scores)

	mut.lm <- lm(diff~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12, data=aggcat)

	# mut.lm <- lm(diff~H3K4me1+H3K4me3+H3K9ac+H3K9me3+H3K27ac+H3K27me3+H3K36me3+
					  # CPGI+EXON+TIME+RATE+prop_GC, data=aggcat)
	
	fits <- mut.lm$fitted.values
	names(fits) <- paste0(aggcat$CHR,".",aggcat$BIN)
	
	BIN <- as.integer(gsub(".*\\.", "", names(fits)))
	CHR <- as.integer(gsub("\\..*", "", names(fits)))
	
	df <- data.frame(CHR, Category2=cat1, BIN,
				   exp=fits, 
				   obs=aggcat$diff,
				   stringsAsFactors=F)
				   
	df$diff <- df$obs-df$exp
	df$pred <- "model"
	
	d.diff <- rbind(d.diff,df)
	
	z <- summary(mut.lm)
	
	print(cat1)
	print(z)
}	



##############################################################################
# Plot error from 5bp motif model predictions vs 5bp+features model predictions
# with ideogram track
##############################################################################
diffm <- melt(d.diff[,1:5], id.vars=c("CHR", "Category2", "BIN"), value.var=c("exp", "obs"))
diffm <- diffm[diffm$CHR==2,]
diffm$BIN2 <- diffm$BIN*100000

p <- plotIdeogram(hg19IdeogramCyto, "chr2", xlabel=TRUE, alpha=0, 
				  zoom.region=c(min(diffm[diffm$Category=="AT_GC",]$BIN2),max(diffm[diffm$Category=="AT_GC",]$BIN2)))

# p2 <- ggplot(diffm[diffm$Category=="GC_AT",], aes(x=BIN2, y=value, colour=variable))+
p2 <- ggplot(diffm, aes(x=BIN, y=value, colour=variable))+
	# geom_bar(stat="identity", position="identity")+
	geom_point(alpha=0.2, size=3)+
	scale_colour_manual("Model", values=myPaletteCat(8)[6:7], labels=c("exp"="5bp+features", "obs"="5bp only"))+
	facet_wrap(~Category2, scales="free", ncol=1)+
	scale_x_continuous(breaks=pretty_breaks(n=25))+
	ylab("Error")+
	xlab(NULL)+
	theme_bw()+
	theme(axis.title.y=element_text(size=16), 
	      axis.text.y=element_text(size=14),
		  legend.title=element_text(size=16),
		  legend.text=element_text(size=14))
		  
		  # ,
		  # axis.text.x=element_blank(), 
		  # axis.ticks.x=element_blank()
	
# tracks(p,p2, heights=c(1.5,7))
ggsave("/net/bipolar/jedidiah/mutation/images/hier_diffs_pt3.png", width=9, height=18)