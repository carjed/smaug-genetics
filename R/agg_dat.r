##############################################################################
# Defines function for aggregating singleton data and obtaining genome-wide
# and per-bin summaries
##############################################################################

aggData <- function(datfile, adj){

	summfile <- datfile$summ
	binfile <- datfile$bin
	mct <- datfile$mct
	nbp <- adj*2+1

	# Get dataframe of observed and predicted counts
	{
		# aggseq <- count(summfile, c("Sequence", "Category2", "COUNT")) #<-plyr
		cat("Generating motif relative rates...\n")
		aggseq <- summfile %>%
			group_by(Sequence, Category2, BIN) %>%
			summarise(n=n()) %>%
			summarise(num=sum(n), mean=mean(n), sd=sd(n))
		aggseq <- merge(aggseq, mct, by="Sequence")
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
