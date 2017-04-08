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
categ <- "AT_CG"

mut_cats <- c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA")

i<-3
for(categ in mut_cats){
	cat("Extracting", categ, "sites...\n")
	# summfile1 <- full_data$sites %>%
	summfile1 <- sites %>%
		filter(Category==categ) %>%
		mutate(Type=gsub("cpg_", "", Category2),
			SEQA=substr(Sequence, cbp-i, cbp+i),
			SEQB=substr(Sequence, cbp*3-i, cbp*3+i),
			Motif=paste0(SEQA, "(", SEQB, ")")) %>%
		dplyr::select(CHR, POS, Sequence=Motif) %>%
		mutate(mut=1)

	for(chr in 1:22){
		posfile <- paste0(parentdir,
				"/output/logmod_data/chr", chr, "_", categ,"_pos_examples.txt")
		dat <- summfile1[summfile1$CHR==chr,]
		dat <- dat[order(dat$POS),]

		write.table(dat, posfile, col.names=F, row.names=F, quote=F, sep="\t")
	}
}

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
# testfile <- paste0(parentdir, "/output/logmod_data/chr22/chr22_",
# 	categ, "_TTTCTTG(CAAGAAA).txt")
# if(!file.exists(testfile)){
builddatouter <- TRUE;
if(builddatouter){

	modtime <- proc.time()
	cat("Building data from training set...\n")

	foreach(chr=1:22) %dopar% {

		# Run with builddat=1 if *_m.txt.gz files do not already exist
		builddat<-0
		if(builddat){

			perlcmd <- paste0("perl ",
				parentdir, "/smaug-genetics/data_mgmt/logit_scripts/getNonMut.pl --b ", catopt,
				" --chr ", chr,
				" --categ ", categ,
				" --bw ", 100,
				" --adj ", adj)
			system(perlcmd)
		}

		fullfile <- paste0(parentdir, "/output/logmod_data/chr",
			chr, "_", categ, "_sites.txt.gz")

		# gzipcmd <- paste0("gzip ", fullfile)
		# system(gzipcmd)

		# Subset chromosome file by motif
		# subcmd <- paste0("zcat ", fullfile, ".gz | head -20000 | awk '{ print >> \"", parentdir,
		# 		"/output/logmod_data/chr", chr, "/chr", chr, "_", categ, "_\" ",
		# 		"$4 \".txt\" }' ")

		subcmd2 <- paste0("zcat ", fullfile, " | awk '{ print >> \"", parentdir,
				"/output/logmod_data/motifs/test/", categ, "_\" ",
				"$4 \".txt\" }' ")

		system(subcmd2)
	}

	tottime <- (proc.time()-modtime)[3]
	cat("Done (", tottime, "s)\n")
}
