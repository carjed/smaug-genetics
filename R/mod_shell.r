#!/usr/bin/Rscript

##############################################################################
# Script for running genome-wide models
# Built under R version 3.3.2
##############################################################################
# Setup: load function helper scripts and packages
# The usePackage function loads packages if they already exist,
# otherwise installs from default CRAN repository
##############################################################################
cat("Loading functions and packages...\n")
source("./R/get_functions.r")

packages <- c("ggplot2", "dplyr", "tidyr", "broom", "RColorBrewer", "MASS",
"speedglm", "boot", "devtools", "psych", "lmtest", "fmsb", "stringr", "hexbin",
"cowplot", "grid", "gridExtra", "gtable", "yaml")

sapply(packages, function(x) suppressMessages(usePackage(x)))

# Load predefined color palettes, once RColorBrewer package is loaded
source("./R/palettes.r")

##############################################################################
# Setup: get arguments from _config.yaml file
##############################################################################
cat("Script will run with the following parameters:\n")

args <- yaml.load_file("./_config.yaml")
argParse(args)

# for(i in 1:length(args)){
#  	##first extract the object value
#  	tempobj=unlist(args[i])
#  	varname=names(args[i])
#
# 	# optional: print args
# 	cat(varname, ":", tempobj, "\n")
#
# 	##now create a new variable with the original name of the list item
# 	eval(parse(text=paste(names(args)[i],"= tempobj")))
# }

# Install the bedr package from github, if needed
# install_github('carjed/bedr')

# Load forked bedr package from custom library path specified in config file
suppressMessages(require(bedr, lib.loc=libpath, quietly=T))

# suppressMessages(usePackage(ggbio))

# Install ggrepel from github for advanced options
# install_github('slowkow/ggrepel')
# suppressMessages(require(ggrepel))

# Define additional variables for cleaner strings, etc.
bink <- binw/1000
nbp <- adj*2+1

datadir <- paste0(parentdir,
	"/output/", nbp, "bp_", bink, "k_singletons_", data)

summfile <- paste0(datadir, "/full.summary")
binfile <- paste0(datadir, "/full_motif_counts_", bin_scheme, ".txt")

##############################################################################
# Read and preprocess data
# returns list of 4 elements:
# - full_data$sites (summarized info for each singleton)
# - full_data$bins (raw bin file)
# - full_data$mct (motif counts genome-wide)
# -full_data$aggseq (initial relative mutation rates per K-mer subtype)
##############################################################################
full_data <- getData(datadir, bin_scheme)

##############################################################################
# Prepare singleton input for logit model
##############################################################################
if(build_logit){
	ptm <- proc.time()
	cat("Preparing data for logistic regression model...\n")

	mut_cats <- c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA")

	i<-3
	for(categ in mut_cats){
		cat("Extracting", categ, "sites...\n")
		summfile1 <- full_data$sites %>%
		# summfile1 <- sites %>%
			filter(Category==categ) %>%
			mutate(Type=gsub("cpg_", "", Category2),
				SEQA=substr(Sequence, cbp-i, cbp+i),
				SEQB=substr(Sequence, cbp*3-i, cbp*3+i),
				Motif=paste0(SEQA, "(", SEQB, ")")) %>%
			dplyr::select(CHR, POS, Sequence=Motif) %>%
			mutate(mut=1)

		for(chr in 1:22){
			posfile <- paste0(parentdir,
					"/output/logmod_data/chr", chr, "_", categ,"_sites.txt")
			dat <- summfile1 %>%
				filter(CHR==chr) %>%
				arrange(POS)

			write.table(dat, posfile, col.names=F, row.names=F, quote=F, sep="\t")
		}
	}
}

##############################################################################
# Get relative mutation rates per subtype; plot as heatmap
##############################################################################
# Summarise aggseq rates for shorter motifs
cbp <- adj+1

ratelist <- list()
testlist <- list()
modlist <- list()
for(j in 1:5){
	i <- j-1

	# using nMotifs in the full file results in incorrect counts when collapsing
	# shorter motifs. Here we merge with binned outputs obtained from running
	# the data_pipeline/augment_summary.pl script with shorter motifs to ensure
	# proper counting of mutable sites
	nbptmp <- i*2+1

	gpdat <- full_data$aggseq %>%
		mutate(Type=gsub("cpg_", "", Category2),
			SEQA=substr(Sequence, cbp-i, cbp+i),
			SEQB=substr(Sequence, cbp*3-i, cbp*3+i),
			Motif=paste0(SEQA, "(", SEQB, ")"))

	if(i==0){
		gpdat <- gpdat %>%
			dplyr::select(Type, Motif, nERVs) %>%
			group_by(Type, Motif) %>%
			summarise(nERVs=sum(nERVs))

		mcfile <- paste0(parentdir, "/output/", nbptmp, "bp_final_rates.txt")
		mcount <- read.table(mcfile, header=T, stringsAsFactors=F)
		mcount <- mcount %>%
			mutate(Motif=ifelse(grepl("^A", Type),
				"A(T)",
				"C(G)")) %>%
			dplyr::select(Type, Motif, nMotifs)

		gpdat <- merge(gpdat, mcount, by=c("Type", "Motif")) %>%
			mutate(ERV_rel_rate=nERVs/nMotifs)
	}
	else if(i>0 & i<4){
 		gpdat<- gpdat %>%
			dplyr::select(Type, Motif, nERVs) %>%
			group_by(Type, Motif) %>%
			summarise(nERVs=sum(nERVs))

		mcfile <- paste0(parentdir, "/output/", nbptmp, "bp_final_rates.txt")
		mcount <- read.table(mcfile, header=T, stringsAsFactors=F)
		mcount <- mcount %>%
			dplyr::select(Type, Motif, nMotifs)

		gpdat <- merge(gpdat, mcount, by=c("Type", "Motif")) %>%
			mutate(ERV_rel_rate=nERVs/nMotifs)

		# Plot heatmap panels
		plotdat <- gpdat %>%
			mutate(v2=substr(Motif,1,i),
				v2a=factor(as.character(lapply(as.vector(v2), reverse_chars))),
				v3=substr(Motif, i+2, i*2+1),
				v4=ERV_rel_rate,
				Category=Type,
				v5=factor(gsub("_", ">", Type)))

		nbox <- length(unique(plotdat$v2a))
	  nint <- nbox/(4^(i-1))
	  xhi <- rep(1:(4^(i-1)),4^(i-1))*nint+0.5
	  xlo <- xhi-nint
	  yhi <- rep(1:(4^(i-1)),each=4^(i-1))*nint+0.5
	  ylo <- yhi-nint
	  f <- data.frame(xlo,xhi,ylo,yhi)

	  levs_a <- as.character(lapply(as.vector(levels(plotdat$v2a)),
			reverse_chars))

		for(k in 1:6){
			categ <- orderedcats[k]
			p1 <- rrheat2(plotdat[plotdat$Category==categ,], f, levs_a, "v5", nbptmp)
			p1a <- p1+theme(legend.position="none")

			png(paste0(parentdir, "/images/", categ, "_", nbptmp, "bp_heatmap.png"),
			height=5, width=5, units="in", res=300)
			pushViewport(viewport(width=unit(5, "in"), height=unit(5, "in")))
			grid.draw(ggplotGrob(p1a))
			dev.off()
		}

		# trim whitespace on panels with imagemagick mogrify
		trimcmd <- paste0("mogrify -trim ",
			parentdir, "/images/*", nbptmp, "*bp_heatmap.png")
		system(trimcmd)

		# extract legend
		legend <- get_legend(p1)
		png(paste0(parentdir, "/images/heatmap_legend.png"),
			height=8, width=3, units="in", res=300)
		grid.draw(legend)
		dev.off()

	}
	else { # don't draw heatmap for 9-mers; use original motif counts
		gpdat <- gpdat %>%
			dplyr::select(Type, Motif, nERVs, nMotifs) %>%
			group_by(Type, Motif) %>%
			summarise(nERVs=sum(nERVs), nMotifs=sum(nMotifs)) %>%
			mutate(ERV_rel_rate=nERVs/nMotifs)
	}

	ratelist[[j]] <- gpdat

	# Test for heterogeneity among subtypes sharing same (K-2)-mer parent
	if(i>0){
		parentdat <- gpdat %>%
			# dplyr::select(Type, Motif, nERVs, nMotifs, rel_prop) %>%
			filter(nERVs > 20) %>%
			mutate(SEQA=substr(Motif, 2, nbptmp-1),
				SEQB=substr(Motif, nbptmp+3, nchar(Motif)-2),
				MotifP=paste0(SEQA, "(", SEQB, ")")) %>%
			group_by(Type, MotifP) %>%
			arrange(Type, MotifP) %>%
			mutate(exp=sum(nERVs)/sum(nMotifs)*nMotifs,
				p=exp/sum(exp),
				n=n(),
				b1=substr(Motif,1,1),
				b2=substr(Motif,nbptmp,nbptmp)) %>%
			filter(n==16)

		moddat <- parentdat %>%
			do(tidy(glance(lm(ERV_rel_rate ~ b1+b2, data=.))))
		modlist[[i]] <- moddat

		hettests <- parentdat %>%
			summarise(pval=chisq.test(nERVs, p=p)$p.value) %>%
			ungroup() %>%
			mutate(fdr=p.adjust(pval, method="fdr"))
		testlist[[i]] <- hettests
	}

	write.table(gpdat,
		paste0(parentdir, "/output/", nbptmp, "bp_final_rates2.txt"),
		col.names=T, row.names=F, quote=F, sep="\t")
}

tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

##############################################################################
# Compare rates between ERVs and common variants
##############################################################################
if(common){
	commondatadir <- paste0(parentdir,
		"/output/", nbp, "bp_", bink, "k_common_", data)
	common_data <- getData(commondatadir, bin_scheme)

	rates7_common <- common_data$aggseq
	cat("Reading common summary file:", summfile, "...\n")
	common_sites <- read.table(common_summfile, header=T, stringsAsFactors=F)
	common_sites$BIN <- ceiling(sites$POS/binw)

	rates7_common <- common_sites %>%
		dplyr::select(CHR, POS, Motif=Sequence, Type=Category, BIN) %>%
		group_by(Type, Motif) %>%
		# summarise(n=n()) %>%
		summarise(nCommon=n())
	rates7_common <- merge(rates7_common, mct, by="Motif")
	rates7_common$Common_rel_rate <- rates7_common$nCommon/rates7_common$nMotifs

	rates7 <- ratelist[[4]]
	r5m <- merge(rates7, rates7_common, by=c("Type", "Motif"))

	avrates <- read.table(paste0(parentdir, "/posterior_7bp.txt"),
		header=T, stringsAsFactors=F, sep="\t")

	names(avrates)[1] <- "Motif"
	avrates$CAT <- paste0(substr(avrates$Motif,4,4), substr(avrates$alt, 4,4))
	avrates$Category[avrates$CAT=="AC" | avrates$CAT=="TG"] <- "AT_CG"
	avrates$Category[avrates$CAT=="AG" | avrates$CAT=="TC"] <- "AT_GC"
	avrates$Category[avrates$CAT=="AT" | avrates$CAT=="TA"] <- "AT_TA"
	avrates$Category[avrates$CAT=="GA" | avrates$CAT=="CT"] <- "GC_AT"
	avrates$Category[avrates$CAT=="GC" | avrates$CAT=="CG"] <- "GC_CG"
	avrates$Category[avrates$CAT=="GT" | avrates$CAT=="CA"] <- "GC_TA"

	avrates$Motif <- paste0(avrates$Motif, "(", avrates$refrev, ")")

	avrates <- avrates %>%
		dplyr::select(Type=Category, Motif, eur)

	r5m <- merge(rates7, avrates, by=c("Type", "Motif"))

	r5m$Category2 <- ifelse(
		substr(r5m$Motif,4,5)=="CG",
		paste0("cpg_",r5m$Type),
		r5m$Type)

	r5m <- r5m %>%
		mutate(Category2 = plyr::mapvalues(Category2, orderedcats1, orderedcats2))

	r5m$Category2 <- factor(r5m$Category2, levels=orderedcats2)

	r5m$prop_diff <- (r5m$eur/(mean(r5m$eur)/mean(r5m$ERV_rel_rate)))/r5m$ERV_rel_rate
	r5m$prop_diff4 <- r5m$prop_diff
	r5m$prop_diff4[r5m$prop_diff< 0.5] <- 0.5
	r5m$prop_diff4[r5m$prop_diff>2] <- 2

	r5m$v2 <- substr(r5m$Motif,1,3)
	r5m$v2a <- as.character(lapply(as.vector(r5m$v2), reverse_chars))
	r5m$v2a <- factor(r5m$v2a)
	r5m$v3 <- substr(r5m$Motif,3+2,3*2+1)

	##############################################################################
	# Plot heatmap of change in relative rates
	##############################################################################
	nbox<-length(unique(r5m$v2a))
	nint<-nbox/4
	xhi <- rep(1:4,4)*nint+0.5
	xlo <- xhi-nint
	yhi <- rep(1:4,each=4)*nint+0.5
	ylo <- yhi-nint
	f <- data.frame(xlo,xhi,ylo,yhi)

	levs <- as.character(lapply(as.vector(levels(r5m$v2a)), reverse_chars))

	for(j in 1:6){
		categ<-orderedcats[j]
		p1a<-rrheat3(r5m[r5m$Type==categ,])

		png(paste0(parentdir, "/images/rare_common_diff2", categ, "_panel.png"),
			height=5, width=5, units="in", res=300)
		pushViewport(viewport(width=unit(5, "in"), height=unit(5, "in")))
		grid.draw(ggplotGrob(p1a))
		dev.off()
	}

	# trim whitespace on panels with imagemagick mogrify
	trimcmd <- paste0("mogrify -trim ",
		parentdir, "/images/rare_common_diff*_panel.png")
	system(trimcmd)

	plegend<-ggplot()+
		geom_tile(data=r5m[r5m$Category==categ,],
			aes(x=v3, y=v2a, fill=prop_diff4))+
		scale_fill_gradientn("Rp/Rs\n",
			colours=myPaletteBrBG(nbp),
			trans="log",
			breaks=c(0.5, 1, 2),
			labels=c("<0.5", "1", ">2"),
			limits=c(0.5, 2.2))+
		xlab("3' flank")+
		ylab("5' flank")

	legend <- get_legend(plegend)
	png(paste0(parentdir, "/images/rare_common_diff2_legend.png"),
		height=8, width=3, units="in", res=300)
	grid.draw(legend)
	dev.off()
}

##############################################################################
# Scripts below assume 7-mers+features model is already complete
##############################################################################
sitefile <- paste0(parentdir, "/output/rocdat.sort_new.txt")
if(file.exists(sitefile)){

	# Validation model
	# -must have parentdir/output/rocdat.7bp.2.txt containing list of sites output
	# from this model and pre-annotated with 7-mer rates from
	# 	-MAC10+ (MU_C)
	# 	-ERVs (MU_S)
	# 	-AV (MU_A)
	ptm <- proc.time()
	cat("Validating models on de novo mutations...\n")
	source("./R/validation.r")
	tottime <- (proc.time()-ptm)[3]
	cat("Done (", tottime, "s)\n")


	# Analyze effects of genomic features
	ptm <- proc.time()
	cat("Analyzing genomic features...\n")
	source("./R/coef_summary.r")
	tottime <- (proc.time()-ptm)[3]
	cat("Done (", tottime, "s)\n")


	# Run negative binomial regression models and plot
	if(negbin_model){
		cat("Initializing negbin regression model...\n")
		ptm <- proc.time()
		source("./R/negbin_mod.r")
		tottime <- (proc.time()-ptm)[3]
		cat("Done. Finished in", tottime, "seconds \n")
	}
} else {
	cat("Looks like you haven't run the logit models yet.
	Refer to the following README files for instructions:

	-data_mgmt/per-site_dp/README.md
	-data_mgmt/logit_scripts/README.md
	-data_mgmt/process_predicted/README.md \n")
}
