##############################################################################
# Builds per-chromosome datasets of all sites in category
# Splits each file and combines by motif to be passed to logit model
#
# To-Do:
# [-] check integration with log_mod.r to ensure parameters passed correctly
##############################################################################

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
# if(!file.exists(testfile)){
builddatouter <- FALSE;
if(builddatouter){

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
				parentdir, "/smaug-genetics/data_mgmt/logit_scripts/getNonMut.pl --b ", catopt,
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
		subcmd <- paste0("zcat ", parentdir, "/output/logmod_data/chr",
			chr, "_", categ, "_m.txt.gz | awk '{ print >> \"", parentdir,
				"/output/logmod_data/chr", chr, "/chr", chr, "_", categ, "_\" ",
				"$4 \".txt\" }' ")

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
