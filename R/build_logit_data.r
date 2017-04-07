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

for(categ in mut_cats){
	cat("Extracting", categ, "sites...\n")
	summfile1 <- full_data$sites %>%
		filter(Category==categ) %>%
		dplyr::select("CHR", "BIN", "POS", "SEQUENCE") %>%
		mutate(mut=1, BIN=ceiling(POS/100000))

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
			chr, "_", categ, "_sites.txt")

		gzipcmd <- paste0("gzip ", fullfile)
		system(gzipcmd)

		# Subset chromosome file by motif
		subcmd <- paste0("zcat ", fullfile, ".gz | awk '{ print >> \"", parentdir,
				"/output/logmod_data/chr", chr, "/chr", chr, "_", categ, "_\" ",
				"$4 \".txt\" }' ")

		system(subcmd)
	}

	tottime <- (proc.time()-modtime)[3]
	cat("Done (", tottime, "s)\n")
}
