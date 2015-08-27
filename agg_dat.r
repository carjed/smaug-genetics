##############################################################################
# Defines function for aggregating singleton data and obtaining genome-wide
# and per-bin summaries
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
	plot_heatmap <- 1
	if(plot_heatmap==1){
	
		# aggseq <- count(summfile, c("Sequence", "Category", "CAT", "COUNT", "SEQ"))
		aggseq <- count(summfile, Sequence, Category, CAT, COUNT, SEQ)
		aggseq$rel_prop <- aggseq$n/aggseq$COUNT
		
		b<-c("A", "C", "G", "T")
		cats<-unique(aggseq$Category)
		for(i in 1:6){
			aggcat <- aggseq[aggseq$Category==cats[i],]
			# aggcat <- aggcat[aggcat$freq>25,]

			# test <- chisq.test(aggcat$obs, p=aggcat$exp/sum(aggcat$exp))

			# cat <- cats[i]
			# msg <- paste0("checking category ",cat)
			# print(msg)


			for(j in 1:4){

				for(k in 1:4){
					ref <- unique(substr(aggcat$Sequence,3,3))
					motif <- paste0(b[j],ref,b[k])
					dat <- aggcat[substr(aggcat$Sequence,2,4)==motif,]
					dat$exp <- dat$COUNT*(sum(dat$n)/sum(dat$COUNT))
					# print(head(dat))
					cat <- cats[i]
					test <- chisq.test(dat$n, p=dat$exp/sum(dat$exp))
					
					if(test$p.value>0.05/96){
						print(cat)
						print(motif)
						print(test$p.value)
						# print(dat)
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
		
		# Plot relative rate heatmaps
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
		if(adj==2){
			for(i in 1:6){
				for(j in 1:512){
					cat<-cats[i]
					seq<-seqs[j]
					dat<-aggseq[(aggseq$Sequence==seqs[j] & aggseq$Category==cats[i]),]
					if(nrow(dat)==2){
						test<-prop.test(dat$n, dat$COUNT)
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
		# aggseq <- count(summfile, c("Sequence", "Category2", "COUNT")) #<-plyr
		agg.gp <- group_by(summfile, Sequence, Category2, COUNT, BIN)
		aggseq <- summarise(agg.gp, n=n()) %>% summarise(num=sum(n), mean=mean(n), sd=sd(n))
		# aggseq2 <- count(summfile, Sequence, Category2, COUNT) #<-dplyr
		aggseq$rel_prop <- aggseq$num/aggseq$COUNT
		
		# ct <- count(summfile,c("CHR","BIN","Category2")) #<-plyr
		ct <- count(summfile, CHR, BIN, Category2) #<-dplyr
		# ct.ord <- arrange(ct, Category2)
		
		binfile$BIN <- paste0(binfile$CHR,".",binfile$BIN)
		
		binsT <- setNames(data.frame(t(binfile[,-c(1:5)])), binfile$BIN)
		binsT$Sequence <- rownames(binsT)

		# merge counts per sequence/category with counts of mutable motifs
		aggseq.m <- merge(aggseq, binsT, by="Sequence")

		# get number of mutable motifs per bin
		atcg.s <- colSums(aggseq.m[aggseq.m$Category2=="AT_CG",-c(1:7)])
		atgc.s <- colSums(aggseq.m[aggseq.m$Category2=="AT_GC",-c(1:7)])
		atta.s <- colSums(aggseq.m[aggseq.m$Category2=="AT_TA",-c(1:7)])
		gcat.s <- colSums(aggseq.m[aggseq.m$Category2=="GC_AT",-c(1:7)])
		gccg.s <- colSums(aggseq.m[aggseq.m$Category2=="GC_CG",-c(1:7)])
		gcta.s <- colSums(aggseq.m[aggseq.m$Category2=="GC_TA",-c(1:7)])
		cpg_gcat.s <- colSums(aggseq.m[aggseq.m$Category2=="cpg_GC_AT",-c(1:7)])
		cpg_gccg.s <- colSums(aggseq.m[aggseq.m$Category2=="cpg_GC_CG",-c(1:7)])
		cpg_gcta.s <- colSums(aggseq.m[aggseq.m$Category2=="cpg_GC_TA",-c(1:7)])
		
		# get expected counts for each row
		aggseq.m[,8:ncol(aggseq.m)] <- aggseq.m$rel_prop*aggseq.m[,8:ncol(aggseq.m)]

		# sum all sequence combinations for each category/bin
		atcg <- colSums(aggseq.m[aggseq.m$Category2=="AT_CG",-c(1:7)])
		atgc <- colSums(aggseq.m[aggseq.m$Category2=="AT_GC",-c(1:7)])
		atta <- colSums(aggseq.m[aggseq.m$Category2=="AT_TA",-c(1:7)])
		gcat <- colSums(aggseq.m[aggseq.m$Category2=="GC_AT",-c(1:7)])
		gccg <- colSums(aggseq.m[aggseq.m$Category2=="GC_CG",-c(1:7)])
		gcta <- colSums(aggseq.m[aggseq.m$Category2=="GC_TA",-c(1:7)])
		cpg_gcat <- colSums(aggseq.m[aggseq.m$Category2=="cpg_GC_AT",-c(1:7)])
		cpg_gccg <- colSums(aggseq.m[aggseq.m$Category2=="cpg_GC_CG",-c(1:7)])
		cpg_gcta <- colSums(aggseq.m[aggseq.m$Category2=="cpg_GC_TA",-c(1:7)])
		
		# aggseq.m<-aggseq.m[(aggseq.m$Category2=="GC_TA" & (substr(aggseq.m$Sequence,2,3)=="TC" | substr(aggseq.m$Sequence,2,3)=="AC")),]


		# vector of expected counts
		exp <- c(atcg, atgc, atta, gcat, gccg, gcta, cpg_gcat, cpg_gccg, cpg_gcta)
		nmotifs <- c(atcg.s, atgc.s, atta.s, gcat.s, gccg.s, gcta.s, cpg_gcat.s, cpg_gccg.s, cpg_gcta.s)

		# vector of categories
		catrep <- c(rep("AT_CG",ncol(aggseq.m)-7), rep("AT_GC",ncol(aggseq.m)-7), rep("AT_TA",ncol(aggseq.m)-7), 
				  rep("GC_AT",ncol(aggseq.m)-7), rep("GC_CG",ncol(aggseq.m)-7), rep("GC_TA",ncol(aggseq.m)-7),
				  rep("cpg_GC_AT",ncol(aggseq.m)-7), rep("cpg_GC_CG",ncol(aggseq.m)-7), rep("cpg_GC_TA",ncol(aggseq.m)-7))
				  
		z<-unlist(strsplit(names(exp), "[.]"))
		CHR<-substring(z[seq(1, length(z), 2)], 4)
		BIN<-z[seq(2, length(z), 2)]
		
		# create data frame with observed, expected, and bin	
		# BIN <- as.integer(gsub(".*\\.", "", names(exp)))
		# CHR <- as.integer(substring(gsub("\\..*", "", names(exp)), 4))
		
		oe2 <- data.frame(CHR, BIN, Category2=catrep, exp, nmotifs)
		oe2 <- merge(oe2,ct, by=c("CHR","Category2","BIN"))
		names(oe2)[6] <- "obs"
		oe2$res <- paste0(nbp,"bp")

		# get odds ratio for each bin/category
		oe2$odds <- oe2$obs/oe2$exp
		oe2$diff <- oe2$obs-oe2$exp
		oe2$Category2 <- as.character(oe2$Category2)

		# order by OR and add column of ranks
		# oe2.ord <- oe2[order(oe2$Category2, oe2$odds),]
		# oe2.ord <- ddply(oe2.ord, .(Category2), transform, rk=seq_along(Category2))
		# oe2.ord$res <- paste0(nbp,"bp")

		# return(oe2.ord)
		
		datalist<- list("agg"=aggseq, "oe"=oe2)
		return(datalist)
	}
}