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
	cat("Updating bin file...\n")
	for(i in 6:ncol(binfile)){
		names(binfile)[i] <- paste0(names(binfile)[i], "(", revcomp(names(binfile)[i]), ")" )
	}

	bins2 <- binfile[,names(binfile)%in%unique(summfile$Sequence)]
	binfile <- cbind(binfile[,1:5],bins2)
	xmax <- floor(max(summfile$BIN)/100)*100

	bins2 <- melt(binfile[,5:ncol(binfile)], id="BIN")
	# bins2$variable<-substr(as.character(bins2$variable),1,5)
	bin.gp <- group_by(bins2, variable)
	bins2 <- summarise(bin.gp, value=sum(value))
	# bins2 <- aggregate(data=bins2, value ~ variable, sum)
	names(bins2) <- c("Sequence", "COUNT")
	bins2$Sequence <- sub("[.]", "(", bins2$Sequence)
	bins2$Sequence <- sub("[.]", ")", bins2$Sequence)
	bins2$SEQ1 <- substr(bins2$Sequence, 0, adj*2+1)
	bins2$SEQ2 <- substr(bins2$Sequence, (adj*2+1)+2, (adj*2+2)+(adj*2+1))
	bins2$SEQMIN <- pmin(bins2$SEQ1, bins2$SEQ2)
	bins2 <- data.frame(bins2$COUNT, bins2$SEQMIN)
	names(bins2) <- c("COUNT", "SEQMIN")

	cat("Updating summary file...\n")
	summfile$SEQMIN <- pmin(summfile$SEQ, summfile$ALTSEQ)
	summfile <- left_join(summfile, bins2, by="SEQMIN")
	
	datalist<- list("summ"=summfile, "bin"=binfile)
	return(datalist)
}
