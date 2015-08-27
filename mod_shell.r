##############################################################################
# Script for running genome-wide models
##############################################################################

##############################################################################
# Process Options/Args
# Define color palettes
# Define functions
##############################################################################
ptm <- proc.time()

cat("Loading functions and packages...\n")
suppressMessages(require(ggplot2))
suppressMessages(require(scales))
# suppressMessages(require(plyr))
suppressMessages(require(dplyr))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))
suppressMessages(require(MASS))
suppressMessages(require(speedglm))
suppressMessages(require(grid))
suppressMessages(require(ggbio))
data(hg19IdeogramCyto, package = "biovizBase")

parentdir<-"/net/bipolar/jedidiah/mutation"

functionfile<-paste0(parentdir, "/smaug-genetics/get_functions.R")
source(functionfile)

binw <- 100000
bink <- binw/1000
adj <- 2
nbp <- adj*2+1

# options--will be updated to specify on cmd line
run_agg <- 1 # get genome-wide aggregate data for relative rates
pcs <- 0 # use PCs covariates
categ <- "AT_GC" # specify mutation category
negbin_model <- 1 # run negbin regression
log_model <- 0 # run logistic regression and predictions
run_predict <- 0

myPaletteCat <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
rb <- c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g <- myPaletteG(6)[1:3]
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

# chr22 <- chr22[-grep(",", chr22$ALT),]

# chr22 <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.expanded.summary", header=T, stringsAsFactors=F)
# bins <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.bin_out.txt", header=T, stringsAsFactors=F)

##############################################################################
# Read in data
##############################################################################
ptm <- proc.time()

summfile <- paste0(parentdir, "/output/", nbp, "bp_", bink, "k/full.summary")
binfile <- paste0(parentdir, "/output/", nbp, "bp_", bink, "k/full_bin.txt")

if(!file.exists(summfile)){
	cat("Merged summary/bin files do not exist---Merging now...\n")
	datadir<-paste0(parentdir, "/output/", nbp, "bp_", bink, "k")
	combinecmd <- paste0("awk 'FNR==1 && NR!=1{while(/^CHR/) getline; } 1 {print} ' ", datadir, "/chr*.expanded.summary > ", datadir, "/full.summary")
	combinecmd2 <- paste0("awk 'FNR==1 && NR!=1{while(/^CHR/) getline; } 1 {print} ' ", datadir, "/chr*.bin_out.txt > ", datadir, "/full_bin.txt")
	system(combinecmd)
	system(combinecmd2)
}

cat("Reading summary file:", summfile, "...\n")
summ_5bp_100k <- read.table(summfile, header=F, stringsAsFactors=F, skip=1)
names(summ_5bp_100k)<-c("CHR", "POS", "REF", "ALT", "DP", "AN", "SEQ", "ALTSEQ", "GC")
summ_5bp_100k$BIN <- ceiling(summ_5bp_100k$POS/binw)

cat("Reading bin file:", binfile, "...\n")
bins_5bp_100k <- read.table(binfile, header=T, stringsAsFactors=F, check.names=F)

tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

##############################################################################
# Update data
##############################################################################
ptm <- proc.time()
cat("Updating data...\n")
source("update_dat.r")
dat_5bp_100k <- updateData(summ_5bp_100k, bins_5bp_100k, adj)
rm(summ_5bp_100k)
rm(bins_5bp_100k)
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

##############################################################################
# Run aggregation script to get:
# -genome-wide heatmaps of relative rates
# -negative binomial model to summarize spatial distribution
# -table of relative rates
##############################################################################
if(run_agg==1){
	ptm <- proc.time()
	cat("Aggregating data...\n")
	source("agg_dat.r")
	aggV <- aggData(dat_5bp_100k, 2) #<-modify the adj value for 3bp data

	agg_5bp_100k <- aggV$oe
	rates1 <- aggV$agg
	
	ratefile <- paste0(parentdir, "/output/", nbp, "bp_", bink, "k_rates.txt")
	write.table(rates1, ratefile, col.names=T, row.names=F, quote=F, sep="\t")

	tottime<-(proc.time()-ptm)[3]
	cat("Done (", tottime, "s)\n")
}

##############################################################################
# Build file of covariates and run PCA
##############################################################################
ptm <- proc.time()
mutcov2file<-paste0(parentdir, "/output/logmod_data/", bink, "kb_mut_cov2.txt")
if(!file.exists(mutcov2file)){
	cat("Building covariate data...\n")
	source("get_covs.r")
} else {
	cat("Reading existing covariate datafile:", mutcov2file, "...\n")
	mut_cov<-read.table(mutcov2file, header=F, stringsAsFactors=F)
}

if(pcs==1){
	names(mut_cov)<-c("CHR", "BIN", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12")
	danames<-c("CHR", "BIN", "POS", "Sequence", "mut", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "cons")
} else {
	names(mut_cov)<-c("CHR", "BIN", "H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3", "H3K27ac", "H3K27me3", "H3K36me3", "CPGI", "EXON", "TIME", "RATE", "prop_GC")
	danames<-c("CHR", "BIN", "POS", "Sequence", "mut", "H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3", "H3K27ac", "H3K27me3", "H3K36me3", "CPGI", "EXON", "TIME", "RATE", "prop_GC", "cons")
}

covnames<-danames[-c(1:5)]

tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

##############################################################################
# Run negative binomial regression models and plot
##############################################################################
if(negbin_model==1){
	cat("Initializing negbin regression model...\n")
	ptm <- proc.time()
	source("negbin_mod.r")
	tottime<-(proc.time()-ptm)[3]
	cat("Done. Finished in", tottime, "seconds \n")
}

##############################################################################
# Run logistic regression model and obtain predicted probabilities
##############################################################################
if(log_model==1){
	cat("Initializing logistic regression model...\n")
	ptm <- proc.time()
	source("log_mod.r")
	tottime<-(proc.time()-ptm)[3]
	cat("Done. Model and predictions finished in", tottime, "seconds \n")
}

if(run_predict==1 & log_model==1){
	outfile<-paste0(parentdir, "/output/predicted/", categ, "_", bink, "kb_out.txt")
	predictcmd<-paste0("perl predict.pl --coefs ", coeffile, " --data ", fullfile, " --out ", outfile)
	system(predictcmd)
}