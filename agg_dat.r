##############################################################################
# Defines function for aggregating singleton data and obtaining genome-wide
# and per-bin summaries
##############################################################################

aggData <- function(datfile, adj){

	summfile <- datfile$summ
	binfile <- datfile$bin
	mct <- datfile$mct
	nbp <- adj*2+1

	# summfile <- dat_5bp_100k$summ
	# binfile <- dat_5bp_100k$bin

	# summfile <- dat_5bp_2$summ
	# binfile <- dat_5bp_2$bin

	# summfile <-summfile[summfile$CHR==2 & summfile$BIN<=600 & summfile$BIN>=300,]
	# binfile <- binfile[binfile$CHR=="chr2" & binfile$BIN<=600 & binfile$BIN>=300,]

	# Compute Ts/Tv by bin--can use as filter set to exclude windows with Ts/Tv
	# below a certain threshold
	cat("Computing Ts/Tv per bin...\n")
	summfile$TS <- ifelse(
		((summfile$REF=="A" | summfile$REF=="G") &
			(summfile$ALT=="A" | summfile$ALT=="G")) |
		((summfile$REF=="C" | summfile$REF=="T") &
			(summfile$ALT=="C" | summfile$ALT=="T")),
		"TS", "TV")

	tstv <- summfile %>%
		group_by(CHR, BIN, TS) %>%
		summarise(freq=n())
	tstv <- dcast(tstv, CHR+BIN~TS, value.var="freq")
	tstv$TOT <- tstv$TS+tstv$TV
	tstv$TSTV <- tstv$TS/tstv$TV

	# Cutoff of 1.5 removes 44Mb and ~500,000 variants
	# 1.2 removes 3Mb and ~60,000 variants (mostly Chr8)

	filter<-0
	if(filter){
		cat("Applying Ts/Tv filter...\n")
		cutoff <- 1.2
		filterset <- tstv[tstv$TSTV<cutoff,]

		summfile$IND <- paste0(summfile$CHR, "_", summfile$BIN)
		filterset$IND <- paste0(filterset$CHR, "_", filterset$BIN)

		summfile2 <- summfile[!(summfile$IND %in% filterset$IND),]
	}

	# Plot genome-wide motif heatmaps
	plot_heatmap <- 1
	if(plot_heatmap==1){
		cat("Generating data for relative rate heatmap...\n")
		# aggseq <- count(summfile, c("Sequence", "Category", "CAT", "COUNT", "SEQ"))
		aggseq <- count(summfile, Sequence, Category, CAT, SEQMIN, SEQ)
		aggseq <- merge(aggseq, mct, by="SEQMIN")
		aggseq$rel_prop <- aggseq$n/aggseq$COUNT

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
		cat("Plotting heatmaps...\n")
		at_heat <- rrheat(map_a1, f, levs_a, "v5", nbp)
		gc_heat <- rrheat(map_g1, f, levs_g, "v5", nbp)

		gwmapfile<-paste0(parentdir, "/images/gw_map.png")
		png(gwmapfile, width=18, height=24, units="in", res=300)
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

		gwmapfile2<-paste0(parentdir, "/images/gw_map_uncollapsed.png")
		png(gwmapfile2, width=48, height=24, units="in", res=300)
		multiplot(a_heat, t_heat, c_heat, g_heat, cols=4)
		dev.off()

		cats<-unique(aggseq$Category)
		seqs<-unique(aggseq$Sequence)

		# Test for differences in proportions on opposite strands
		pvals<-data.frame()
		if(adj==2){
			for(i in 1:6){
				for(j in 1:512){
					cat<-cats[i]
					seq<-seqs[j]
					dat<-aggseq[(aggseq$Sequence==seqs[j] & aggseq$Category==cats[i]),]
					if(nrow(dat)==2){
						test<-prop.test(dat$n, dat$COUNT)
						row<-data.frame(cat, seq, test$p.value)
						pvals<-rbind(pvals, row)
					}
				}
			}
		}

		pvals$adj<-p.adjust(pvals$test.p.value, "BH")
		head(pvals)

	}

	# Get dataframe of observed and predicted counts
	{
		# aggseq <- count(summfile, c("Sequence", "Category2", "COUNT")) #<-plyr
		cat("Generating motif relative rates...\n")
		aggseq <- summfile %>%
			group_by(Sequence, Category2, SEQMIN, BIN) %>%
			summarise(n=n()) %>%
			summarise(num=sum(n), mean=mean(n), sd=sd(n))
		aggseq <- merge(aggseq, mct, by="SEQMIN")
		# aggseq2 <- count(summfile, Sequence, Category2, COUNT) #<-dplyr
		aggseq$rel_prop <- aggseq$num/aggseq$COUNT

		# ct <- count(summfile,c("CHR","BIN","Category2")) #<-plyr
		# ct <- count(summfile, CHR, BIN, Category2) #<-dplyr
		# ct.ord <- arrange(ct, Category2)

		# NEW VERSION -- need to update so nmotifs is consistent for all 3 categories
		# currently, motifs without a present singleton in a given bin are not
		# accounted for
		cat("Counting mutations per motif per bin...\n")
		summagg <- summfile %>%
				dplyr::select(CHR, POS, BIN, Sequence, Category2) %>%
				group_by(CHR, BIN, Sequence, Category2) %>%
				summarise(obs=n())

		# summagg <- merge(summagg, aggseq[,c("Sequence", "Category2", "COUNT")],
		# 	by=c("Sequence", "Category2"), all=TRUE)

		# binfile <- gather(binfile, Sequence, Count, 6:ncol(binfile))
		# binfile$CHR <- as.integer(substring(binfile$CHR, 4))
		# binfile <- binfile %>% arrange(substring(Sequence, adj+1, adj+1))
		#
		# b2 <- binfile[rep(seq_len(nrow(binfile)), each=3),]
		# b2$Category <- c(rep(c("AT_CG", "AT_GC", "AT_TA"), nrow(b2)/6),
		# 				  rep(c("GC_AT", "GC_CG", "GC_TA"), nrow(b2)/6))
		#
		# b2$Category2 <- ifelse(substr(b2$Sequence,adj+1,adj+2)=="CG",
		# 						paste0("cpg_",b2$Category),
		# 						b2$Category)
		# b2 <- b2 %>% arrange(CHR, BIN)
		# {
		# row.names(b2) <- paste0(b2$CHR, "_", b2$BIN, "_", b2$Sequence, "_", b2$Category2)
		# row.names(summagg) <- paste0(summagg$CHR, "_",
			# summagg$BIN, "_", summagg$Sequence, "_", summagg$Category2)

		# summagg2 <- cbind(summagg,
			# count=b2[,"Count"][match(rownames(summagg), rownames(b2))])
		# }
		# summagg2 <- merge(b2, summagg, by=c("CHR", "BIN", "Category2", "Sequence"), all=T)
		# summagg2$exp <- summagg2$rel_prop*summagg2$Count

		# Count observed singletons per bin per basic category
		spbc <- summagg %>%
			group_by(CHR, BIN, Category2) %>%
			summarise(obs=sum(obs, na.rm=T))

		# s2 <- summagg2 %>%
		# 	group_by(CHR, BIN, Category2) %>%
		# 	summarise(exp=sum(exp, na.rm=T),
		# 		obs=sum(obs, na.rm=T),
		# 		nmotifs=sum(Count, na.rm=T),
		# 		motifvar=var(Count),
		# 		maxn=max(Count),
		# 		minn=min(Count))

		# summcor <- summagg2 %>%
					# group_by(Category2, Sequence) %>%
					# summarise(cor=cor(exp, obs, use="complete.obs"))

		datalist<- list("agg"=aggseq, "oe"=spbc, "summagg2"=summagg)
		return(datalist)
	}
}
