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

		# NEW VERSION -- need to update so nmotifs is consistent for all 3 categories
		# currently, motifs without a present singleton in a given bin are not
		# accounted for
		cat("Counting mutations per motif per bin...\n")
		summagg <- summfile %>%
				dplyr::select(CHR, POS, BIN, Sequence, Category2) %>%
				group_by(CHR, BIN, Sequence, Category2) %>%
				summarise(obs=n())

		# Count observed singletons per bin per basic category
		spbc <- summagg %>%
			group_by(CHR, BIN, Category2) %>%
			summarise(obs=sum(obs, na.rm=T))

		datalist<- list("agg"=aggseq, "oe"=spbc, "summagg2"=summagg)
		return(datalist)
	}
}
