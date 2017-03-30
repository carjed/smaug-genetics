##############################################################################
# Add Bin/Category/Sequence columns to summary file, update colnames of bin
# file, and merge overall count for each sequence motif into summary file
##############################################################################
updateData <- function(summfile, binfile, adj){

	nbp <- adj*2+1

	# summfile$BIN <- ceiling(summfile$POS/binw)


	# get complement of sequence columns in bin file and remove duplicates
	# change this to an apply()-based function?

	# cat("Updating bin file columns...\n")
	# for(i in 6:ncol(binfile)){
	# 	names(binfile)[i] <- paste0(names(binfile)[i], "(", revcomp(names(binfile)[i]), ")" )
	# }
	#
	# cat("Removing redundant columns from bin file...\n")
	# bins2 <- binfile[,names(binfile)%in%unique(summfile$Sequence)]
	# binfile <- cbind(binfile[,1:5],bins2)
	# xmax <- floor(max(summfile$BIN)/100)*100

	# cat("Counting total motifs in genome...\n")
	# mct <- melt(binfile[,5:ncol(binfile)], id="BIN") %>%
	# 	group_by(variable) %>%
	# 	summarise(value=sum(value))
	# # bins2 <- aggregate(data=bins2, value ~ variable, sum)
	# names(mct) <- c("Sequence", "COUNT")
	# mct$Sequence <- sub("[.]", "(", mct$Sequence)
	# mct$Sequence <- sub("[.]", ")", mct$Sequence)
	# mct$SEQ1 <- substr(mct$Sequence, 0, adj*2+1)
	# mct$SEQ2 <- substr(mct$Sequence, (adj*2+1)+2, (adj*2+2)+(adj*2+1))
	# mct$SEQMIN <- pmin(mct$SEQ1, mct$SEQ2)
	# mct <- mct %>%
	# 	dplyr::select(SEQMIN, COUNT)

	mct <- bins_5bp_100k %>% 
		dplyr::select(CHR, Sequence=MOTIF, COUNT) %>%
		group_by(Sequence) %>%
		summarise(COUNT=sum(COUNT))


	datalist<- list("summ"=summfile, "bin"=binfile, "mct"=mct)
	return(datalist)
}
