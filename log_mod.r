##############################################################################
# Logistic regression model
##############################################################################

suppressMessages(require(data.table))
suppressMessages(require(foreach))
suppressMessages(require(doSNOW))

cluster <- makeCluster(8, type = "SOCK", outfile="")
registerDoSNOW(cluster)


binw <- 1000000
bink <- binw/1000

# Target mutation rate
mu <- 1e-8

# number of non-N bases
nbases <- 5.8e9

# number of individuals in sample
nind <- 3612

# estimated number of mutations per generation--
# 58 corresponds to a rate of ~1e-8 per Michaelson et al.
indmuts <- mu*nbases

# mean number of singletons per individual
nsing <- 10000

# estimate of maximum age (in generations) of a singleton
numgens <- ceiling(nsing/indmuts)

# estimate average age of singleton (~100 gens)
meangens <- mean(1:numgens)
denom <- nind*meangens

# Fast list of 6 basic categories from agg_5bp_100k data
mut_cats <- unique(agg_5bp_100k$Category2[nchar(agg_5bp_100k$Category2)==5])

# subset reference data to only AT or GC bases
catopt <- substr(categ,0,2)

# summfile1 <- dat_5bp_100k$summ[(dat_5bp_100k$summ$POS>=6000000 &
	# dat_5bp_100k$summ$POS<=7000000 &
	# dat_5bp_100k$summ$CHR=="20" &
	# dat_5bp_100k$summ$Category==categ), c("CHR", "POS", "BIN", "Sequence")]

##############################################################################
# Get summary file for category i; merge with covariates
##############################################################################
cat("Extracting", categ, "sites...\n")
summfile1 <- dat_5bp_100k$summ[dat_5bp_100k$summ$Category==categ,
	c("CHR", "BIN", "POS", "Sequence")]
summfile1$mut <- 1
summfile1$BIN <- ceiling(summfile1$POS/100000)
# summfile1 <- merge(summfile1, mut_cov, by=c("CHR", "BIN"))

##############################################################################
# Run perl script to build input data for model
# -outputs 1 row per site with reference base same as category
# -singletons have mut=1, otherwise 0
# -covariates are PCs from the 100kb mut_cov2 file
##############################################################################
# trainchr <- seq(1,10,2)
trainchr1 <- c(1:22)
nchr <- length(trainchr1)

# fullfile <- paste0(parentdir, "/output/logmod_data/",categ,"_full.txt")
# Only subset if specified temp file does not exist
testfile <- paste0(parentdir, "/output/logmod_data/chr22/chr22_",
	categ, "_TTTCTTG(CAAGAAA).txt")
if(!file.exists(testfile)){

	modtime <- proc.time()
	cat("Building data from training set...\n")

	foreach(chr=1:22) %dopar% {

		# Run with builddat=1 if *_m.txt.gz files do not already exist
		builddat<-0
		if(builddat){
			posfile <- paste0(parentdir,
					"/output/logmod_data/chr", chr, "_", categ,"_pos_examples.txt")
			dat <- summfile1[summfile1$CHR==chr,]
			dat <- dat[order(dat$POS),]

			write.table(dat, posfile, col.names=F, row.names=F, quote=F, sep="\t")

			perlcmd <- paste0("perl ",
				parentdir, "/smaug-genetics/getNonMut.pl --b ", catopt,
				" --chr ", chr,
				" --categ ", categ,
				" --bw ", 100,
				" --adj ", adj)
			system(perlcmd)

			# Gets position list for exons in specified chr
			chrcmd <- paste0("awk \'$1 == \"chr", chr, "\" { print }\' ",
				parentdir, "/reference_data/GRCh37_RefSeq_chop.bed | sort -k3,3 | cut -f 3,4 > ",
				parentdir, "/reference_data/chr", chr, "_RefSeq.bed")

			system(chrcmd)

			# Join exon positions with full data
			# cuts out regional exon info (column 14) and adds the binary column
			exoncmd <- paste0("join -a1 -1 3 -2 1 -e 0 -o auto -t $\'\\t\' ",
				parentdir, "/output/logmod_data/chr", chr, "_", categ, "_sites.txt ",
				parentdir, "/reference_data/chr", chr, "_RefSeq.bed | ",
				" cut -f 14 --complement > ",
				parentdir, "/output/logmod_data/chr", chr, "_", categ, "_m.txt")

			system(exoncmd)
		}

		# Subset chromosome file by motif
		subcmd <- paste0("zcat ", parentdir, "/output/logmod_data/chr", chr, "_", categ, "_m.txt.gz | ",
			"awk '{ print >> \"",
				parentdir, "/output/logmod_data/chr", chr, "/chr", chr, "_", categ, "_\" ", "$4 \".txt\" }' ")

		system(subcmd)
	}

	# Run Unix command to combine data, ordered by chromosome (1 to 22)
	# catcmd1 <- paste0("ls -v ", parentdir,
	# 	"/output/logmod_data/chr*_", categ, "_sites.txt | xargs cat >> ",
	# 	fullfile)
	# system(catcmd1)
	tottime <- (proc.time()-modtime)[3]
	cat("Done (", tottime, "s)\n")
}

##############################################################################
# Subset data by motif (using grep in system command) and run motif-specific
# model, creating data frame of parameter estimates
# (to be passed to prediction function)
##############################################################################
cat("Running model...\n")

coefdat <- data.frame(stringsAsFactors=F)
# int_only_rates <- data.frame(stringsAsFactors=F)
motifs <- sort(unique(summfile1$Sequence))

coefdat<-foreach(i=1:length(motifs), .combine=rbind) %dopar% {
# coefdat<-foreach(i=1:16, .combine=rbind) %dopar% {
	# motif <- substr(motifs[i], 0, nbp)
	motif <- motifs[i]
	cat("Running model", i, "on", motif, "sites...\n")

	require(speedglm)

	# Shortened motif
	escmotif <- substr(motif, 0, nbp)

	# Define name of temporary file for motif i
	tmpfile <- paste0(parentdir, "/output/logmod_data/motifs/",
		categ, "_", escmotif, ".txt")

	perchrtmp <- paste0(parentdir,
		"/output/logmod_data/chr*/chr*_", categ, "_", motif, ".txt")

	# Merge per-chromosome motif files to single file
	if(!(file.exists(tmpfile))){

		cat("Merging ", motif, " files...\n")

		catcmd1 <- paste0("find ", parentdir, "/output/logmod_data/chr* -name '*",
			escmotif, "*.txt' | sort -V | xargs cat >> ", tmpfile)
		system(catcmd1)
	}

	# Remove per-chromosome motif files once merged
	# unlink(perchrtmp)

	da1 <- read.table(tmpfile, header=F, stringsAsFactors=F)
	names(da1) <- c("CHR", "BIN", "POS", "Sequence", "mut",
		danames[-c(1:2,11)], "EXON")

	# Run logit model for categories with >10 singletons, return coefficients
	# Otherwise, returns single marginal rate
	if(sum(da1$mut)>10){
		log_mod_formula <- as.formula(paste("mut~",
			paste(names(da1)[-(1:5)], collapse="+")))
		log_mod <- speedglm(log_mod_formula, data=da1, family=binomial(), maxit=50)

		# Get coefficients from model summary, clean up data formats
		coefs<-data.frame(summary(log_mod)$coefficients, stringsAsFactors=F)
		coefs <- cbind(Cov = rownames(coefs), coefs)
		coefs$Cov<-as.character(coefs$Cov)
		rownames(coefs)<-NULL

		coefs[,-1] <- data.frame(apply(coefs[,-1], 2, function(x) as.numeric(as.character(x))))
		names(coefs)<-c("Cov", "Estimate", "SE", "Z", "pval")

		coefs$Sequence<-escmotif
		coefs
	} else {
		cat("Not enough data--using marginal rate only\n")
		alt <- c(sum(da1$mut)/nrow(da1), rep(0,13))
		# alt
	}

	# Remove motif file once model finished
	# unlink(tmpfile)
}

coeffile <- paste0(parentdir,
	"/output/logmod_data/", categ, "_", bink, "kb_coefs_p.txt")
write.table(coefdat, coeffile, col.names=F, row.names=F, quote=F, sep="\t")

##############################################################################
# If prediction is specified, creates and executes a slurm batch file
# to run predictions over all chromosomes specified in set
##############################################################################
if(run_predict){

	predictchr <- c(1:22)
	predictstr <- paste(predictchr, collapse=",")

	buildbatchcmd <- paste0("perl ",
		parentdir, "/smaug-genetics/build_batch.pl --trchr ", predictstr,
		" --cat ", categ,
		" --bink ", bink)
	system(buildbatchcmd)

	slurmcmd<-paste0("sbatch ", parentdir, "/smaug-genetics/slurm_predict.txt")
	system(slurmcmd)
}
