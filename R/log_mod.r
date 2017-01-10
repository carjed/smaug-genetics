#!/usr/bin/Rscript

##############################################################################
# Logistic regression model
##############################################################################

libpath <- "~/R/x86_64-pc-linux-gnu-library/2.13"

options(useHTTPS=FALSE)
suppressMessages(require(speedglm, lib.loc=libpath, quietly=T))
suppressMessages(require(bedr, lib.loc=libpath, quietly=T))
suppressMessages(require(dplyr, lib.loc=libpath, quietly=T))
suppressMessages(require(BSgenome.Hsapiens.UCSC.hg19, lib.loc=libpath, quietly=T))
suppressMessages(require(Repitools, lib.loc=libpath, quietly=T))
suppressMessages(require(boot, lib.loc=libpath, quietly=T))

source("./get_functions.r")

args <- getArgs(
	defaults=list(
		jobid=25,
		categ="AT_CG",
		nmotifs=4096,
		nodes=10))

nbp <- 7
# parentdir <- "/net/bipolar/jedidiah/mutation"

cat("Script will run with the following parameters:\n")
for(i in 1:length(args)){
	##first extract the object value
	tempobj=unlist(args[i])
	varname=names(args[i])

	# optional: print args
	cat(varname, ":", tempobj, "\n")

	##now create a new variable with the original name of the list item
	eval(parse(text=paste(names(args)[i],"= tempobj")))
}
cat("\n")

jobid <- as.numeric(jobid)

# cluster <- makeCluster(nodes, type = "SOCK", outfile="paste0(parentdir, "/snow.log"))
# registerDoSNOW(cluster)

# Target mutation rate
mu <- 1e-8
nbases <- 5.8e9
nind <- 3612
indmuts <- mu*nbases
nsing <- 10000
numgens <- ceiling(nsing/indmuts)
meangens <- mean(1:numgens)
denom <- nind*meangens

# Fast list of 6 basic categories from agg_5bp_100k data
mut_cats <- c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA")

# subset reference data to only AT or GC bases
catopt <- substr(categ,0,2)

motiffile <- paste0(parentdir, "/output/7bp_1000k_rates.txt")
incmd <- paste0("cut -f1-3 ", motiffile)
motifdat <- read.table(pipe(incmd), header=T, stringsAsFactors=F)
motifs <- motifdat %>%
	mutate(Category=gsub("cpg_", "", Category2)) %>%
	filter(Category==categ) %>%
	dplyr::select(Sequence) %>%
	unlist

runmotif <- motifs[jobid]
escmotif <- substr(runmotif, 0, nbp)
# comb <- function(x, ...) {
#       mapply(rbind,x,...,SIMPLIFY=FALSE)
# }

##############################################################################
# Run models
##############################################################################
cat("Running model on", runmotif, "sites...\n")
coefs <- logitMod(motif=runmotif, nbp=nbp, parentdir=parentdir, categ=categ)

# covlist <- clusterApply(cluster, motifs[1:nmotifs], logitMod, nbp=nbp, parentdir=parentdir, categ=categ)
# fullcoef <- rbind_all(covlist)

coefdir <- paste0(parentdir, "/output/logmod_data/coefs/", categ, "/")
dir.create(coefdir, recursive=T)
coeffile <- paste0(coefdir, categ, "_", escmotif, "_coefs.txt")
write.table(coefs, coeffile, col.names=F, row.names=F, quote=F, sep="\t")
