##############################################################################
# Script for running genome-wide models
##############################################################################

##############################################################################
# Function to install and load packages
##############################################################################
usePackage <- function(p) {
	
	p <- as.character(substitute(p))

    if (!is.element(p, installed.packages()[,1])){
        install.packages(p, dep = TRUE)
		require(p, character.only = TRUE)
	} else {
		require(p, character.only = TRUE)
	}
}

##############################################################################
# Function to read file from disk
##############################################################################
make.data<-function(filename, chunksize,...){
	conn<-NULL
	function(reset=FALSE){
		if(reset){
			if(!is.null(conn)) close(conn)
			conn<<-file(filename,open="r")
		} else{
			rval<-read.table(conn, nrows=chunksize,...)
			if ((nrow(rval)==0)) {
				close(conn)
				conn<<-NULL
				rval<-NULL
			}
			return(rval)
		}
	}
}

suppressMessages(usePackage(ggplot2))
suppressMessages(usePackage(scales))
suppressMessages(usePackage(plyr))
suppressMessages(usePackage(reshape2))
suppressMessages(usePackage(RColorBrewer))
suppressMessages(usePackage(MASS))
suppressMessages(usePackage(speedglm))
suppressMessages(usePackage(grid))

source("/net/bipolar/jedidiah/mutation/smaug-genetics/get_functions.R")

binw <- 100000

myPaletteCat <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
rb <- c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g <- myPaletteG(6)[1:3]

summ_5bp_100k <- read.table("/net/bipolar/jedidiah/mutation/output/5bp_100k/full.summary", header=T, stringsAsFactors=F)

summ_5bp_100k$BIN <- ceiling(summ_5bp_100k$POS/binw)

bins_5bp_100k <- read.table("/net/bipolar/jedidiah/mutation/output/5bp_100k/full_bin.txt", header=T, stringsAsFactors=F, check.names=F)

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
rm(summ_5bp_100k)
rm(bins_5bp_100k)
# gc()
# dat_5bp_10k<-updateData(summ_5bp_10k, bins_5bp_10k, 2)
# dat_3bp_100k<-updateData(summ_3bp_100k, bins_3bp_100k, 1)
# dat_3bp_10k<-updateData(summ_3bp_10k, bins_3bp_10k, 1)

##############################################################################
# Build matrix of genomic features
##############################################################################

# Function reads in .bed file and returns chr, bin, and value for each feature
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

# Number of reads of each histone mark per 100kb
H3K4me1 <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K4me1.100kb.bed", "H3K4me1")
H3K4me3 <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K4me3.100kb.bed", "H3K4me3")
H3K9ac <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K9ac.100kb.bed", "H3K9ac")
H3K9me3 <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K9me3.100kb.bed", "H3K9me3")
H3K27ac <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K27ac.100kb.bed", "H3K27ac")
H3K27me3 <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K27me3.100kb.bed", "H3K27me3")
H3K36me3 <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/H3K36me3.100kb.bed", "H3K36me3")

# Pct in CpG island per 100kb
CPGI <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/cpg_islands_100kb.bed", "CPGI")
CPGI$CPGI <- CPGI$CPGI/binw

# Pct exonic per 100kb
EXON <- getHistBins("/net/bipolar/jedidiah/mutation/reference_data/coding_exon_100kb.bed", "EXON")
EXON$EXON <- EXON$EXON/binw

# Replication timing per 100kb
reptime <- read.table("/net/bipolar/jedidiah/mutation/reference_data/lymph_rep_time.txt", header=F, stringsAsFactors=F, sep="\t")
names(reptime) <- c("CHR", "POS", "TIME")
reptime$BIN <- ceiling(reptime$POS/binw)
rtagg <- aggregate(TIME~CHR+BIN, reptime, mean)

# Recombination rate per 100kb
rcrate <- read.table("/net/bipolar/jedidiah/mutation/reference_data/recomb_rate.bed", header=T, stringsAsFactors=F, sep="\t")
rcrate$CHR <- as.integer(substring(rcrate$CHR, 4))
rcrate$POS <- (rcrate$START+rcrate$END)/2
rcrate$BIN <- ceiling(rcrate$POS/binw)
rcagg <- aggregate(RATE~CHR+BIN, rcrate, mean)

# GC content per 100kb
pctgc <- dat_5bp_100k$bin[,c(1,5,4)]
pctgc$CHR <- as.integer(substring(pctgc$CHR, 4))

# Merge all data
mut_cov <- cbind(H3K4me1, H3K4me3[,2], H3K9ac[,2], H3K9me3[,2], 
			   H3K27ac[,2], H3K27me3[,2], H3K36me3[,2], CPGI[,2], EXON[,2])
names(mut_cov) <- c("CHR", "H3K4me1", "BIN", "H3K4me3", "H3K9ac", "H3K9me3", 
				  "H3K27ac", "H3K27me3", "H3K36me3", "CPGI", "EXON")
mut_cov <- merge(mut_cov, rtagg, by=c("CHR", "BIN"))
mut_cov <- merge(mut_cov, rcagg, by=c("CHR", "BIN"))
mut_cov <- merge(mut_cov, pctgc, by=c("CHR", "BIN"))

# Run PCA on features, write table of chr, bin, and all PCs (rounded to 3 sig. figs)
pc.dat <- prcomp(mut_cov[,3:ncol(mut_cov)], center=T, scale=T)
scores <- as.data.frame(round(pc.dat$x, 3))
mut_cov2<-cbind(mut_cov[,1:2], scores)
write.table(mut_cov2, "/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/mut_cov2.txt", sep="\t", col.names=F, row.names=F, quote=F) 

##############################################################################
# Logistic regression model
##############################################################################

log_model <- 1
plots<-data.frame()
mut_cats<-unique(dat_5bp_100k$summ$Category)
if(log_model==1){

	for(i in 1:1){
		categ<-mut_cats[i]
		catopt<-substr(categ,0,2)
		# summfile1 <- dat_5bp_100k$summ[(dat_5bp_100k$summ$POS>=6000000 & dat_5bp_100k$summ$POS<=7000000 & dat_5bp_100k$summ$CHR=="20" & dat_5bp_100k$summ$Category==categ), c("CHR", "POS", "BIN", "Sequence")]
		
		# Get summary file for category i; merge with covariates
		summfile1 <- dat_5bp_100k$summ[dat_5bp_100k$summ$Category==categ, c("CHR", "BIN", "POS", "Sequence")]
		summfile1$mut <- 1
		summfile1 <- merge(summfile1, mut_cov2, by=c("CHR", "BIN"))
		
		posfile<-paste0("/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/",categ,"_pos_examples.txt")
		write.table(summfile1, posfile, col.names=F, row.names=F, quote=F, sep="\t")

		# Run perl script to get negative examples, merging with covariates from mut_cov2 file
		# Extracts positions from each chromosome to exclude from negative examples
		for(chr in 1:22){	
			z <- summfile1[summfile1$CHR==chr,]$POS
			exclfile<-paste0("/net/bipolar/jedidiah/mutation/smaug-genetics/chr",chr,"_",categ,"_exclusion_list.txt")
			write.table(z, exclfile, col.names=F, row.names=F, quote=F)
			
			perlcmd <- paste0("perl /net/bipolar/jedidiah/mutation/smaug-genetics/getNonMut.pl --b ", catopt, " --chr ", chr, " --categ ", categ)
			system(perlcmd)
		}
		
		# Run cat command to combine +/- data
		catcmd1 <- paste0("cat /net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/chr*_",categ,"_negative_examples.txt > /net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/",categ,"_negative_examples.txt")
		system(catcmd1)
		catcmd2 <- paste0("cat /net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/",categ,"_pos_examples.txt /net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/",categ,"_negative_examples.txt > /net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/",categ,"_full.txt")
		system(catcmd2)
		
		# Read in full dataset for category i and run logistic regression model
		danames<-c("CHR", "BIN", "POS", "Sequence", "mut", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12")
		fullfile<-paste0("/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/",categ,"_full.txt")
		da<-make.data(fullfile,chunksize=1000000,col.names=danames)

		log_mod<-shglm(mut~factor(Sequence)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12, datafun=da, family=binomial())
		
		# Read small subset for dummy model
		da2<-read.table("/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/full_sort.txt", header=F, stringsAsFactors=F, nrows=1000000)
		names(da2)<-danames

		# Further reduce data sample for glm, ensuring all 256 sequence factors are represented in dataset
		numseqs<-0
		
		while(numseqs!=256){
			da3<-da2[sample(nrow(da2), 20000),]
			numseqs<-length(unique(da3$Sequence))
		}
		
		# Run dummy model
		log_mod_null<-glm(mut~Sequence+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12, data=da3, family=binomial)
		
		# Create copy of dummy model
		log_mod_null2<-log_mod_null
		
		# Update coefficient estimates with those from full model
		log_mod_null2$coefficients<-log_mod$coefficients
		
		# Iterate over full data 1MB at a time and output predicted relative rates for each site
		chunksize<-1000000
		numsites<-1146613132 #<- need to update numsites for each of the 6 categories
		numchunks<-ceil(numsites/chunksize)

		# fullfile2<-"/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/AT_CG_full.txt"

		for(i in 1:numchunks){

			tmp<-read.table(fullfile, header=F, stringsAsFactors=F, nrows=chunksize, skip=((i-1)*chunksize))
			names(tmp)<-danames
			
			preds<-plogis(predict(log_mod_null3, tmp))
			
			out<-cbind(tmp[,1:3], preds)
			
			outpath<-paste0("/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/", categ, "_", i, ".txt")
			
			write.table(out, outpath, col.names=F, row.names=F, quote=F, sep="\t")
		}
		
		# OLD COMMANDS--used for plotting histogram of predicted values
		# z<-data.frame(V=logit(log_mod$linear.predictors))
		# z$cat<-categ
		# fithist<-ggplot(z, aes(x=V))+geom_histogram()+ggtitle(categ)+xlab("P(singleton)")
		# plots<-rbind(plots, z)
	}
	
	# multiplot(plots[1], plots[2], plots[3], plots[4], plots[5], plots[6], cols=3)
	# fithist<-ggplot(plots, aes(x=V))+geom_histogram()+ggtitle("mutation distributions")+xlab("P(singleton)")+facet_wrap(~cat, scales="free")
	# ggsave("/net/bipolar/jedidiah/mutation/images/fitted_hist.png")
}