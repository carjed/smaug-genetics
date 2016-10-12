#! /usr/bin/env Rscript

##############################################################################
# Compare correlation between 100kb windowed distribution of singletons
# and common variants.
# -Should be merged with mutation_method_corr.R
##############################################################################

#!/usr/bin/Rscript

##############################################################################
# Script for running genome-wide models
# Built under R version 3.2.2
#
# Usage:
# ./mod_shell.r --summfile="/path/to/summaryfile"
#				--binfile="/path/to/binfile"
#				--binw=binwidth (integer)
# 				--adj=number of bases in either direction to look
# 				--run_agg=TRUE: aggregate data
# 				--pcs=FALSE: use principal components for modeling
# 				--negbin_model==TRUE: run negbin regression
# 				--log_model=FALSE: run logit regression
# 				--run_predict=FALSE: get predicted values
# 				--categ="NN_NN" category to use for logit model
##############################################################################

##############################################################################
# Setup: process args, load function helper script, load packages
##############################################################################
ptm <- proc.time()

parentdir<-dirname(getwd())
cat("Loading functions and packages...\n")

source("R/get_functions.r")

# Get args from command line; defaults defined below
args <- getArgs(
	defaults=list(adj=3,
		binw=1000000,
		summfile=paste0(parentdir, "/output/7bp_1000k_common/full.summary"),
		binfile=paste0(parentdir, "/output/7bp_1000k_common/full_bin.txt"),
		# summfile=paste0(parentdir, "/output/5bp_100k/full.summary"),
		# binfile=paste0(parentdir, "/output/5bp_100k/full_bin.txt"),
		# summfile=paste0(parentdir, "/output/7bp_1000k/chrX.expanded.summary"),
		# binfile=paste0(parentdir, "/output/7bp_1000k/chrX.bin_out.txt"),
		run_agg=TRUE,
		pcs=FALSE,
		categ="GC_CG",
		negbin_model=TRUE,
		log_model=FALSE,
		run_predict=FALSE))

# The usePackage function loads packages if they already exist,
# otherwise installs from default CRAN repository
suppressMessages(usePackage(ggplot2))
suppressMessages(usePackage(dplyr))
suppressMessages(usePackage(tidyr))
suppressMessages(usePackage(reshape2))
suppressMessages(usePackage(RColorBrewer))
suppressMessages(usePackage(MASS))
suppressMessages(usePackage(speedglm))
suppressMessages(usePackage(boot))
suppressMessages(usePackage(devtools))
require(psych)
# suppressMessages(usePackage(ggbio))

# Manual toggle for installing ggbio package
# Uses the install_github() function from devtools to pull latest version,
# due to issue on some clusters where using
# "source("https://bioconductor.org/biocLite.R")"
# does not properly update the installer
#
# If this is not an issue, can simply run:
# source("https://bioconductor.org/biocLite.R")
# biocLite("ggbio", suppressUpdates=TRUE)
install_ggbio <- 0
if(install_ggbio){
	install_github("Bioconductor/BiocInstaller")
	biocLite("ggbio", suppressUpdates=TRUE)
}

# Get ideogram data for plotting track under chromosome
# data(hg19IdeogramCyto, package = "biovizBase")

# Parse args--command line options become objects, so instead of using
# "args$summfile", we can just use "summfile"
# Works as a lightweight version of attach(args), but explicitly defines objects
# in future version, automatically parse numeric variables properly
# Modified from: http://www.r-bloggers.com/extract-objects-from-a-list/

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

# Need to manually coerce binwidth and adj args to numeric
binw <- as.numeric(binw)
adj <- as.numeric(adj)

# Define additional variables for cleaner strings, etc.
bink <- binw/1000
nbp <- adj*2+1

datadir <- dirname(summfile)

# Define color palettes
myPaletteCat <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")

myPaletteRdBu <- colorRampPalette(rev(brewer.pal(9, "RdBu")), space="Lab")
myPaletteOrRd <- colorRampPalette(rev(brewer.pal(9, "OrRd")), space="Lab")
myPalettePurples <- colorRampPalette(rev(brewer.pal(9, "Purples")), space="Lab")
myPaletteBrBG <- colorRampPalette(rev(brewer.pal(9, "BrBG")), space="Lab")

rb <- c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g <- myPaletteG(6)[1:3]
rbg<-c(myPaletteB(6)[1:3],myPaletteR(6)[1:3], myPaletteG(6)[1:3])

darken <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
}

tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

##############################################################################
# Read in data
##############################################################################
ptm <- proc.time()

if(!file.exists(summfile)){
	cat("Merged summary/bin files do not exist---Merging now...\n")

	# Change ^SEQ to ^CHR--needed to fix bug in header of common variant data
	combinecmd <- paste0(
		"awk 'FNR==1 && NR!=1{while(/^SEQ/) getline; } 1 {print} ' ",
		datadir, "/chr*.expanded.summary > ", datadir, "/full.summary")
	combinecmd2 <- paste0(
		"awk 'FNR==1 && NR!=1{while(/^CHR/) getline; } 1 {print} ' ",
		datadir, "/chr*.bin_out.txt > ", datadir, "/full_bin.txt")
	system(combinecmd)
	system(combinecmd2)
}

cat("Reading summary file:", summfile, "...\n")

summ_5bp_100k <- read.table(summfile, header=F, stringsAsFactors=F, skip=1)
names(summ_5bp_100k)<-c(
	"CHR", "POS", "REF", "ALT", "DP", "AN", "SEQ", "ALTSEQ", "GC")

summ_5bp_100k$BIN <- ceiling(summ_5bp_100k$POS/binw)

bins_5bp_100k <- read.table(binfile, header=T, stringsAsFactors=F, check.names=F)

source("R/update_dat.r")
dat_5bp_100k <- updateData(summ_5bp_100k, bins_5bp_100k, adj)
rm(summ_5bp_100k)
rm(bins_5bp_100k)

source("R/agg_dat.r")
aggV <- aggData(dat_5bp_100k, adj) #<-modify the adj value for 3bp data

agg_5bp_100k <- aggV$oe
rates5 <- aggV$agg
summagg2 <- aggV$summagg2

write.table(rates5, "/net/bipolar/jedidiah/mutation/output/5bp_100k_common_rates.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(agg_5bp_100k, "/net/bipolar/jedidiah/mutation/output/5bp_100k_common_bins.txt", col.names=T, row.names=F, quote=F, sep="\t")

agg_5bp_100k_common<-agg_5bp_100k
rates5_common<-rates5

##############################################################################
##############################################################################
##############################################################################
# Repeat with singletons
##############################################################################
##############################################################################
##############################################################################
summfile=paste0(parentdir, "/output/7bp_1000k/full_j.summary")
binfile=paste0(parentdir, "/output/7bp_1000k/full_bin.txt")

summ_5bp_100k <- read.table(summfile, header=F, stringsAsFactors=F, skip=1)
summ_5bp_100k <- summ_5bp_100k[sample(nrow(summ_5bp_100k), 12088037),]
names(summ_5bp_100k)<-c(
	"CHR", "POS", "REF", "ALT", "DP", "AN", "SEQ", "ALTSEQ", "GC")

summ_5bp_100k$BIN <- ceiling(summ_5bp_100k$POS/binw)

bins_5bp_100k <- read.table(binfile, header=T, stringsAsFactors=F, check.names=F)

##############################################################################
# Update data
##############################################################################
dat_5bp_100k <- updateData(summ_5bp_100k, bins_5bp_100k, adj)
rm(summ_5bp_100k)
rm(bins_5bp_100k)

aggV <- aggData(dat_5bp_100k, adj) #<-modify the adj value for 3bp data

agg_5bp_100k <- aggV$oe
rates5 <- aggV$agg
summagg2 <- aggV$summagg2

##############################################################################
# Plot correlation comparing per-bin counts
##############################################################################
names(agg_5bp_100k_common)[4]<-"common_obs"
bincts<-merge(agg_5bp_100k, agg_5bp_100k_common, by=c("CHR", "BIN", "Category2"))

bincts %>% group_by(Category2) %>% summarise(cor=cor(obs, common_obs))

ggplot(data=bincts,
		aes(x=obs, y=common_obs, colour=Category2, group=Category2))+
	geom_point(alpha=0.4)+
	geom_smooth(method=lm, se=FALSE, colour="black")+
	scale_colour_manual(values=rbg)+
	xlab("Singletons")+
	ylab("Polymorphisms (MAC>10)")+
	# xlab("Group 1 Relative Rate")+
	# ylab("Group 2 Relative Rate")+
	theme_bw()+
	theme(legend.position="none")+
	# geom_text(data=dat, aes(x=-Inf, y=Inf, label=val, hjust=0, vjust=1), colour="black")+
	facet_wrap(~Category2, scales="free")
ggsave("/net/bipolar/jedidiah/mutation/images/sing_com_corr.png", width=6.2, height=4.2)

##############################################################################
# Plot heatmap of change in relative rates
##############################################################################
names(rates5_common)[8]<-"common_rel_prop"
r5m<-merge(rates5[,c(2,3,4,7,8)],
	rates5_common[,c(2,3,4,7,8)], by=c("Sequence", "Category2")) %>%
	mutate(Category=gsub("cpg_", "", Category2)) %>%
	#group_by(Category) %>%
	mutate(nsing=sum(num.x),
		ncommon=sum(num.y),
		prop_s=num.x/nsing,
		prop_p=num.y/ncommon,
		prop_diff=prop_p/prop_s)

r5m$prop_diff3<-r5m$prop_diff
r5m$prop_diff3[r5m$prop_diff< 0.5]<- 0.5
r5m$prop_diff3[r5m$prop_diff>2]<- 2

r5m$prop_diff4<-r5m$num.y/r5m$num.x
r5m$prop_diff4<-r5m$prop_diff4/mean(r5m$prop_diff4)

r5m$prop_diff5<-r5m$prop_diff4
r5m$prop_diff4[r5m$prop_diff4< 0.5]<- 0.5
r5m$prop_diff4[r5m$prop_diff4>2]<- 2

r5m$v2 <- substr(r5m$Sequence,1,adj)
r5m$v2a <- as.character(lapply(as.vector(r5m$v2), reverse_chars))
r5m$v2a <- factor(r5m$v2a)
r5m$v3 <- substr(r5m$Sequence,adj+2,adj*2+1)

nbox<-length(unique(r5m$v2a))
nint<-nbox/4
xhi <- rep(1:4,4)*nint+0.5
xlo <- xhi-nint
yhi <- rep(1:4,each=4)*nint+0.5
ylo <- yhi-nint
f <- data.frame(xlo,xhi,ylo,yhi)

levs <- as.character(lapply(as.vector(levels(r5m$v2a)), reverse_chars))

ggplot()+
	geom_tile(data=r5m, aes(x=v2, y=v3, fill=prop_diff4))+
	geom_rect(data=f, size=0.8, colour="grey30",
		aes(xmin=xlo, xmax=xhi, ymin=ylo, ymax=yhi), fill=NA)+
	scale_fill_gradientn("Rp/Rs\n",
		colours=myPaletteBrBG(nbp),
		trans="log",
		breaks=c(0.5, 1, 2),
		labels=c("<0.5", "1", ">2"),
		limits=c(0.5, 2.2))+
	xlab("5' flank")+
	ylab("3' flank")+
	theme(legend.title = element_text(size=18),
	  legend.text = element_text(size=16),
	  strip.text.x = element_text(size=20),
	  axis.title.x = element_text(size=20),
	  axis.title.y = element_text(size=20),
		axis.text.y = element_text(size=6, colour="black"),
	  axis.text.x = element_text(size=6, colour="black", angle=90, hjust=1))+
		# axis.text.y = element_blank(),
		# axis.text.x = element_blank())+
	scale_x_discrete(labels=levs)+
	facet_wrap(~Category, ncol=3, scales="free_x")
ggsave("/net/bipolar/jedidiah/mutation/images/rare_common_diff2.png",  height=12, width=24)

# t-test for significant difference between Rs/Rp
r5ma<-r5m
r5ma$pval<-0
for(i in 1:nrow(r5m)){
  row<-r5ma[i,]
  r5ma[i,]$pval<-prop.test(c(row$num.x, row$num.y),
    c(row$COUNT.x, row$COUNT.y))$p.value
}

##############################################################################
# Test for heteroscedasticity
##############################################################################
orderedcats<-c("AT_CG", "AT_GC", "AT_TA",
	"GC_AT", "GC_CG", "GC_TA",
	"cpg_GC_AT", "cpg_GC_CG", "cpg_GC_TA")

dat<-data.frame()
for(i in orderedcats){
  r5m1 <- r5m %>% filter(Category2==i)
  r5m2 <- r5m1 %>% filter(num.x>=100 & num.y>=100)
  p1 <- lm(rel_prop~common_rel_prop+num.x+num.y, data=r5m1)
  p2 <- lm(rel_prop~common_rel_prop+num.x+num.y, data=r5m2)
  t1 <- bptest(p1)
  t2 <- bptest(p2)
  fullcor <- cor(r5m1$rel_prop, r5m1$common_rel_prop, method="spearman")
  subcor <- cor(r5m2$rel_prop, r5m2$common_rel_prop, method="spearman")

  fullcor <- cor(r5m1$rel_prop, r5m1$common_rel_prop)
  subcor <- cor(r5m2$rel_prop, r5m2$common_rel_prop)

  newrow1 <- data.frame(Category=i, Data="full", Cor=fullcor, Dim=nrow(r5m1), pval=t1$p.value)
  newrow2 <- data.frame(Category=i, Data="sub", Cor=subcor, Dim=nrow(r5m2), pval=t2$p.value)
  dat<-rbind(dat, newrow1)
  dat<-rbind(dat, newrow2)
}


r5m$Category2<-factor(r5m$Category2, levels=orderedcats)
ggplot(r5m, aes(x=Category2, y=prop_diff5, fill=Category2))+
	geom_violin()+
	geom_boxplot(width=0.2,alpha=0.4)+
	geom_hline(yintercept=1, linetype="longdash")+
	ylab("Rs/Rc")+
	theme_bw()+
	theme(legend.position="none",
		axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
		axis.title.x=element_blank())
ggsave("/net/bipolar/jedidiah/mutation/images/rare_common_dist.png")

nsing<-sum(r5m$num.x)
ncommon<-sum(r5m$num.y)
pvals<-rep(0, nrow(r5m))
for(i in 1:nrow(r5m)){
	row<-r5m[i,]

	test<-prop.test(c(ceiling(row$num.x/3.073), row$num.y), c(ceiling(nsing/3.073), ncommon))
	pvals[i]<-test$p.value
}
r5m$pval<-pvals

# Plot scatterplots of correlation
r5m %>% group_by(Category2) %>% summarise(cor=cor(rel_prop, common_rel_prop))

myPaletteCat <- colorRampPalette(brewer.pal(12, "Paired"))
cols <- myPaletteCat(12)[c(8,10,12,2,4,6,1,3,5)]

r5m$log.Rp.Rs <- ifelse(abs(log(r5m$prop_diff5)) > 0.223, ">0.223", "<0.223")
ggplot(data=r5m, aes(x=rel_prop, y=common_rel_prop, group=Category2, colour=log.Rp.Rs))+
  geom_point(alpha=0.3, size=2)+
  # geom_point(data=dat[dat$cluster!=2,], aes(x=eur, y=rel_prop, group=Category2), alpha=0.3)+
  # geom_point(data=dat[dat$cluster==2,], aes(x=eur, y=rel_prop, group=Category2, colour=factor(ats)), alpha=0.3)+
  # geom_point(alpha=0.3)+
  geom_smooth(method="lm", se=F)+
	geom_abline(intercept=0, linetype="dashed")+
  #geom_text_repel(data=dat[dat$cluster==2,], aes(label=SEQUENCE, colour=factor(ats)))+
  facet_wrap(~Category2, scales="free")+
  xlab("Relative rate (ERVs)")+
  ylab("Relative rate (Common)")+
	scale_x_log10()+
	scale_y_log10()+
  # scale_colour_manual(values=cols)+
	scale_colour_brewer("|log(Rp/Rs)|", palette="Dark2")+
  theme_classic()+
  theme(#legend.position="none",
    strip.text.x=element_text(size=14),
    legend.title=element_text(size=16),
    axis.title.x=element_text(size=18),
    axis.title.y=element_text(size=18))
ggsave("/net/bipolar/jedidiah/mutation/images/sing_com_cor_rates.png", width=9, height=9)

r5m<-r5m %>%
	mutate(Category2 = plyr::mapvalues(Category2, orderedcats1, orderedcats2))

r5m$Category2 <- factor(r5m$Category2, levels=orderedcats2)

r5m2 <- r5m %>% dplyr::select(-Category2)
ggplot(data=r5m, aes(x=rel_prop, y=common_rel_prop))+
  geom_point(data=r5m2, alpha=0.3, size=2, colour="grey70")+
	geom_point(aes(group=Category2, colour=Category2), alpha=0.3, size=2)+
  xlab("Relative rate (ERVs)")+
  ylab("Relative rate (Common)")+
	geom_smooth(method="lm", se=F, colour="black")+
	geom_abline(intercept=0, linetype="dashed")+
	scale_x_log10()+
	scale_y_log10()+
	scale_colour_manual("Mutation Type", values=cols)+
	facet_wrap(~Category2, scales="free", dir="v")+
  theme_bw()+
  theme(
		legend.position="none",
    strip.text.x=element_text(size=14),
    legend.title=element_text(size=16),
    axis.title.x=element_text(size=18),
    axis.title.y=element_text(size=18))
ggsave("/net/bipolar/jedidiah/mutation/images/sing_com_cor_rates_facet.png", width=9, height=9)

ggplot(data=r5m, aes(x=rel_prop, y=common_rel_prop))+
  # geom_point(data=r5m2, alpha=0.3, size=2, colour="grey70")+
	geom_point(aes(group=Category2, colour=Category2), alpha=0.3, size=2)+
  xlab("Relative rate (ERVs)")+
  ylab("Relative rate (Common)")+
	scale_x_log10()+
	scale_y_log10()+
	scale_colour_manual("Mutation Type", values=cols)+
	guides(colour = guide_legend(nrow=3, override.aes = list(alpha=1)))+
	# facet_wrap(~Category2, scales="free")+
  theme_bw()+
  theme(
		legend.position="bottom",
    strip.text.x=element_text(size=14),
    legend.title=element_text(size=16),
    axis.title.x=element_text(size=18),
    axis.title.y=element_text(size=18))
ggsave("/net/bipolar/jedidiah/mutation/images/sing_com_cor_rates_no_facet.png", width=6, height=6)


pcadat <- prcomp(r5m[r5m$Category2=="A>G",c(5,8)])

dat <- data.frame(cbind(pcadat$x, r5m[r5m$Category2=="A>G",]))

# k<-kmeans(pca$x, 5)
# dat$kmcl<-k$cluster
dat$cluster <- ifelse(dat$PC2>0.003, 2,1)
dat$Sequence <- substr(dat$Sequence, 0, 7)
dat$cgs <- ifelse(substr(dat$Sequence, 1,2)%in%c("CC", "CG", "GG", "GC"), 1,0)
dat$ts <- ifelse(grepl("T", dat$Sequence), 0,1)
dat$num_GC_flank <- ifelse(str_count(dat$Sequence, "A|T")<=3,">=4","<4")
dat$maxc <- ifelse(abs(log(dat$prop_diff5))>0.223, "N", "Y")
dat$maxc <- factor(dat$maxc, levels=c("Y", "N"))

ggplot(data=dat,
		aes(x=rel_prop, y=common_rel_prop, group=Category2, colour=maxc, shape=num_GC_flank))+
  geom_point(alpha=0.3, size=4)+
	scale_colour_brewer(palette="Set2")+
	# geom_point(data=dat[dat$num_AT_flank!=">2",], alpha=0.3, size=2)+
  # geom_point(data=dat[dat$cluster!=2,], aes(x=eur, y=rel_prop, group=Category2), alpha=0.3)+
  # geom_point(data=dat[dat$cluster==2,], aes(x=eur, y=rel_prop, group=Category2, colour=factor(ats)), alpha=0.3)+
  # geom_point(alpha=0.3)+
  #geom_text_repel(data=dat[dat$cluster==2,], aes(label=SEQUENCE, colour=factor(ats)))+
  #facet_wrap(~Category2, scales="free")+
  xlab("Relative rate (ERVs)")+
  ylab("Relative rate (Common)")+
	guides(shape = guide_legend("#G/C bases in flanking region"),
		colour = guide_legend("0.8<Rp/Rs<1.25 ?", override.aes = list(alpha=1)))+
  # scale_colour_brewer("# G/C bases in flank", palette="Dark2")+
  theme_bw()+
  theme(legend.position = "bottom",
    legend.title=element_text(size=16),
    axis.title.x=element_text(size=18),
    axis.title.y=element_text(size=18))
ggsave("/net/bipolar/jedidiah/mutation/images/AT_GC_scatter.png", width=8, height=6)

##############################################################################
# Simulate correlations
##############################################################################
output<-data.frame()
for(i in 1:50){
	cat("running simulation", i, "of 50...\n")
	simdat1<-dat_5bp_100k$summ[sample(nrow(dat_5bp_100k$summ), 24176074),]
	simdat1$gp<-sample(0:1, nrow(simdat1), replace=T)

	simdat<-simdat1 %>%
		group_by(CHR, BIN, Category2, gp) %>%
		summarise(count=n())

	simdat$gp<-ifelse(simdat$gp==0,"gp1", "gp2")
	simdat2<-simdat %>% spread(gp, count)
	simdat2<-simdat2[complete.cases(simdat2),]

	simcor<-simdat2 %>% group_by(Category2) %>% summarise(cor=cor(gp1, gp2))
	output<-rbind(output, simcor)
}

simcor<-output %>%
	group_by(Category2) %>%
	summarise(corsim=mean(cor)) %>%
	mutate(gp="sim")

obscor<-bincts %>%
	group_by(Category2) %>%
	summarise(corobs=cor(obs, common_obs)) %>%
	mutate(gp="obs")

tmp<-r.test(2897, simcor$corsim, obscor$corobs)$z

corcomb<-cbind(merge(simcor, obscor, by="Category2"), tmp) %>%
	mutate(cor.p=-2*pnorm(-abs(tmp), log.p=T))

names(simcor)[2]<-"cor"
names(obscor)[2]<-"cor"
corplot<-rbind(simcor, obscor)


corplot$Category2<-factor(corplot$Category2, levels=orderedcats)
corplot$gp<-as.factor(corplot$gp)
corplot$gp<-relevel(corplot$gp, "sim")
corplot<-corplot %>%
	group_by(Category2) %>%
	mutate(cormax=max(cor)) %>%
	arrange(Category2)

corplot$xst<-seq(0.75,9.25,0.5)
corplot$xend<-rep(1:9,each=2)

corplot2<-merge(corplot, corcomb)
corplot2$cor.p<-round(corplot2$cor.p, 1)

ggplot(corplot2, aes(x=Category2, y=cor, fill=gp))+
	geom_bar(position="dodge", stat="identity")+
	geom_segment(data=corplot2, aes(x=xst, xend=xst, yend=cormax+0.1))+
	geom_segment(data=corplot2, aes(x=xst, xend=xend, y=cormax+0.1, yend=cormax+0.1))+
	geom_text(data=corplot2, aes(y=cormax+0.2, label=cor.p))+
	scale_fill_manual("Correlation", values=myPaletteCat(4)[1:2], labels=c("Simulated", "Observed"))+
	scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1))+
	# geom_errorbar(limits, position="dodge", width=1)+
	theme_classic()+
	theme(legend.title=element_blank(),
		  legend.text=element_text(size=16),
		  axis.title.x=element_text(size=20),
		  axis.title.y=element_text(size=20),
		axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
		  axis.text.y=element_text(size=16))+
	xlab("Category")+
	ylab("Correlation")
ggsave("/net/bipolar/jedidiah/mutation/images/sing_com_cor_bars.png", width=9.6, height=4.8)
