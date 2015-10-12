##############################################################################
# Logistic regression model
##############################################################################

suppressMessages(require(data.table))
suppressMessages(require(foreach))
suppressMessages(require(doSNOW))

# Target mutation rate
mu<-1e-8

# number of non-N bases
nbases<-5.8e9

# number of individuals in sample
nind<-3612

# estimated number of mutations per generation--
# 58 corresponds to a rate of ~1e-8 per Michaelson et al.
indmuts<-mu*nbases

# mean number of singletons per individual
nsing<-10000

# estimate of maximum age (in generations) of a singleton
numgens<-ceiling(nsing/indmuts)

# estimate average age of singleton (~100 gens)
meangens<-mean(1:numgens)
denom<-nind*meangens

# Fast list of 6 basic categories from agg_5bp_100k data
mut_cats<-unique(agg_5bp_100k$Category2[nchar(agg_5bp_100k$Category2)==5])

# subset reference data to only AT or GC bases
catopt<-substr(categ,0,2)

# summfile1 <- dat_5bp_100k$summ[(dat_5bp_100k$summ$POS>=6000000 &
	# dat_5bp_100k$summ$POS<=7000000 &
	# dat_5bp_100k$summ$CHR=="20" &
	# dat_5bp_100k$summ$Category==categ), c("CHR", "POS", "BIN", "Sequence")]

##############################################################################
# Get summary file for category i; merge with covariates
##############################################################################
if(!exists("summfile1")){
	cat("Extracting", categ, "sites...\n")
	summfile1 <- dat_5bp_100k$summ[dat_5bp_100k$summ$Category==categ,
		c("CHR", "BIN", "POS", "Sequence")]
	summfile1$mut <- 1
	# summfile1 <- merge(summfile1, mut_cov, by=c("CHR", "BIN"))
}

##############################################################################
# Run perl script to build input data for model
# -outputs 1 row per site with reference base same as category
# -singletons have mut=1, otherwise 0
# -covariates are PCs from the 100kb mut_cov2 file
##############################################################################
# trainchr <- seq(1,10,2)
trainchr <- c(20:22)
trainstr<-paste(trainchr, collapse=",")
nchr<-length(trainchr)

fullfile <- paste0(parentdir, "/output/logmod_data/",categ,"_full.txt")

if(!file.exists(fullfile)){

	modtime<-proc.time()
	cat("Building data from training set...\n")

	mutcov2file<-paste0(parentdir, "/output/logmod_data/100kb_mut_cov2.txt")

	for(chr in trainchr){

		posfile<-paste0(parentdir,
				"/output/logmod_data/chr", chr, "_", categ,"_pos_examples.txt")
		dat<-summfile1[summfile1$CHR==chr,]
		dat<-dat[order(dat$POS),]

		write.table(dat, posfile, col.names=F, row.names=F, quote=F, sep="\t")

		perlcmd <- paste0("perl ", parentdir,
			"/smaug-genetics/getNonMut.pl
				--b ", catopt,
				" --chr ", chr,
				" --categ ", categ,
				" --bw ", bink,
				" --covs ", mutcov2file)
		system(perlcmd)
	}

	# Run Unix command to combine data, ordered by chromosome (1 to 22)
	catcmd1 <- paste0("ls -v ", parentdir,
		"/output/logmod_data/chr*_", categ, "_sites.txt | xargs cat >> ", fullfile)
	system(catcmd1)
	tottime<-(proc.time()-modtime)[3]
	# cat("Done (", tottime, "s)\n")
}

##############################################################################
# Subset data by motif (using grep in system command) and run motif-specific
# model, creating data frame of parameter estimates
# (to be passed to prediction function)
##############################################################################
cat("Running model...\n")

coefdat<-data.frame(stringsAsFactors=F)
motifs<-sort(unique(summfile1$Sequence))
for(i in 1:length(motifs)){
	motif<-substr(motifs[i], 0, 5)
	# cat("Running model", i, "on", motif, "sites...\n")
	modtime <- proc.time()

	tmpfile <- paste0(parentdir, "/output/logmod_data/", categ, "_tmp.txt")
	grepcmd <- paste0("grep ", motif, " ", fullfile, " > ", tmpfile)
	system(grepcmd)

	da1<-read.table(tmpfile, header=F, stringsAsFactors=F)
	names(da1)<-danames

	log_mod_int <- speedglm(mut~., data=da1, family=binomial(), maxit=50)

	log_mod_formula<-as.formula(paste("mut~", paste(covnames, collapse="+")))
	log_mod<-speedglm(log_mod_formula, data=da1, family=binomial(), maxit=50)

	z<-as.numeric(log_mod$coefficients)
	coefdat<-rbind(coefdat, z)

	tottime<-(proc.time()-modtime)[3]
	cat("Finished category", i, "of 256", " (", tottime, "s)\n")
}

coefdat<-cbind(motifs, coefdat)
names(coefdat)<-c("Sequence", "(Intercept)", covnames)

coeffile <- paste0(parentdir,
	"/output/logmod_data/", categ, "_", bink, "kb_coefs.txt")
write.table(coefdat, coeffile, col.names=F, row.names=F, quote=F, sep="\t")

##############################################################################
# If prediction is specified, creates and executes a slurm batch file
# to run predictions over all chromosomes specified in set
##############################################################################
if(run_predict){
	buildbatchcmd<-paste0("perl ",
		parentdir, "/smaug-genetics/build_batch.pl --trchr ",
		trainstr, " --cat ", categ)
	system(buildbatchcmd)

	slurmcmd<-paste0("sbatch ", parentdir, "/smaug-genetics/slurm_predict.txt")
	system(slurmcmd)
}
