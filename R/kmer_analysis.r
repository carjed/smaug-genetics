# Summarise aggseq rates for shorter motifs
ptm <- proc.time()

cbp <- adj+1

ratelist <- list()
testlist <- list()
modlist <- list()
for(j in 1:5){
	i <- j-1

	# using nMotifs in the full file results in incorrect counts when collapsing
	# shorter motifs. Here we merge with binned outputs obtained from running
	# the data_pipeline/augment_summary.pl script with shorter motifs to ensure
	# proper counting of mutable sites
	nbptmp <- i*2+1
	bindir <- paste0(parentdir, "/motif_counts/", nbptmp, "-mers/full")

	gpdat <- full_data$aggseq %>%
		mutate(Type=gsub("cpg_", "", Category2),
			SEQA=substr(Motif, cbp-i, cbp+i),
			SEQB=substr(Motif, cbp*3-i, cbp*3+i),
			Motif=paste0(SEQA, "(", SEQB, ")"))

	if(i==0){
		gpdat <- gpdat %>%
			dplyr::select(Type, Motif, nERVs) %>%
			group_by(Type, Motif) %>%
			summarise(nERVs=sum(nERVs))

		mcfile <- paste0(parentdir, "/output/", nbptmp, "bp_final_rates.txt")
		mcount <- read.table(mcfile, header=T, stringsAsFactors=F)
		mcount <- mcount %>%
			mutate(Motif=ifelse(grepl("^A", Type),
				"A(T)",
				"C(G)")) %>%
			dplyr::select(Type, Motif, nMotifs)

		gpdat <- merge(gpdat, mcount, by=c("Type", "Motif")) %>%
			mutate(ERV_rel_rate=nERVs/nMotifs)
	}
	else if(i>0 & i<4){
 		gpdat<- gpdat %>%
			dplyr::select(Type, Motif, nERVs) %>%
			group_by(Type, Motif) %>%
			summarise(nERVs=sum(nERVs))

		bins <- get_bins(bindir)
		mcount <- get_mct(bins)
		# mcfile <- paste0(parentdir, "/output/", nbptmp, "bp_final_rates.txt")
		# mcount <- read.table(mcfile, header=T, stringsAsFactors=F)
		# mcount <- mcount %>%
		# 	dplyr::select(Type, Motif, nMotifs)

		gpdat <- merge(gpdat, mcount, by=c("Motif")) %>%
			mutate(ERV_rel_rate=nERVs/nMotifs)

		# Plot heatmap panels
		plotdat <- gpdat %>%
			mutate(v2=substr(Motif,1,i),
				v2a=factor(as.character(lapply(as.vector(v2), reverse_chars))),
				v3=substr(Motif, i+2, i*2+1),
				v4=ERV_rel_rate,
				Category=Type,
				v5=factor(gsub("_", ">", Type)))

		nbox <- length(unique(plotdat$v2a))
	  nint <- nbox/(4^(i-1))
	  xhi <- rep(1:(4^(i-1)),4^(i-1))*nint+0.5
	  xlo <- xhi-nint
	  yhi <- rep(1:(4^(i-1)),each=4^(i-1))*nint+0.5
	  ylo <- yhi-nint
	  f <- data.frame(xlo,xhi,ylo,yhi)

	  levs_a <- as.character(lapply(as.vector(levels(plotdat$v2a)),
			reverse_chars))

		for(k in 1:6){
			categ <- orderedcats[k]
			p1 <- rrheat2(plotdat[plotdat$Category==categ,], f, levs_a, "v5", nbptmp)
			p1a <- p1+theme(legend.position="none")

			png(paste0(parentdir, "/images/", categ, "_", nbptmp, "bp_heatmap.png"),
			height=5, width=5, units="in", res=300)
			pushViewport(viewport(width=unit(5, "in"), height=unit(5, "in")))
			grid.draw(ggplotGrob(p1a))
			dev.off()
		}

		# trim whitespace on panels with imagemagick mogrify
		trimcmd <- paste0("mogrify -trim ",
			parentdir, "/images/*", nbptmp, "bp_heatmap.png")
		system(trimcmd)

		# extract legend
		legend <- get_legend(p1)
		png(paste0(parentdir, "/images/heatmap_legend.png"),
			height=8, width=3, units="in", res=300)
		grid.draw(legend)
		dev.off()

	}
	else { # don't draw heatmap for 9-mers; use original motif counts
		gpdat <- gpdat %>%
			dplyr::select(Type, Motif, nERVs, nMotifs) %>%
			group_by(Type, Motif) %>%
			summarise(nERVs=sum(nERVs), nMotifs=sum(nMotifs)) %>%
			mutate(ERV_rel_rate=nERVs/nMotifs)
	}

	ratelist[[j]] <- gpdat

	# Test for heterogeneity among subtypes sharing same (K-2)-mer parent
	if(i>0){
		parentdat <- gpdat %>%
			# dplyr::select(Type, Motif, nERVs, nMotifs, rel_prop) %>%
			filter(nERVs > 20) %>%
			mutate(SEQA=substr(Motif, 2, nbptmp-1),
				SEQB=substr(Motif, nbptmp+3, nchar(Motif)-2),
				MotifP=paste0(SEQA, "(", SEQB, ")")) %>%
			group_by(Type, MotifP) %>%
			arrange(Type, MotifP) %>%
			mutate(exp=sum(nERVs)/sum(nMotifs)*nMotifs,
				p=exp/sum(exp),
				n=n(),
				b1=substr(Motif,1,1),
				b2=substr(Motif,nbptmp,nbptmp)) %>%
			filter(n==16)

		moddat <- parentdat %>%
			do(tidy(glance(lm(ERV_rel_rate ~ b1+b2, data=.))))
		modlist[[i]] <- moddat

		hettests <- parentdat %>%
			summarise(pval=chisq.test(nERVs, p=p)$p.value) %>%
			ungroup() %>%
			mutate(fdr=p.adjust(pval, method="fdr"))
		testlist[[i]] <- hettests
	}

	write.table(gpdat,
		paste0(parentdir, "/output/rates/", nbptmp, "bp_final_rates2.txt"),
		col.names=T, row.names=F, quote=F, sep="\t")
}

tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
