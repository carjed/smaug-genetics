#! /usr/bin/env Rscript

##############################################################################
# Setup: process args, load function helper script, load packages
##############################################################################
ptm <- proc.time()

parentdir<-dirname(getwd())
cat("Loading functions and packages...\n")

source("R/get_functions.r")

adj <- 3
binw <- 1000000
summfile <- paste0(parentdir, "/output/7bp_1000k_common/full.summary")
binfile <- paste0(parentdir, "/output/7bp_1000k_common/full_bin.txt")

# Define additional variables for cleaner strings, etc.
bink <- binw/1000
nbp <- adj*2+1

datadir <- dirname(summfile)

tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

##############################################################################
# Read in MAC10+ data
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

agg_5bp_100k_common <- agg_5bp_100k
rates5_common <- rates5
names(rates5_common)[8] <- "common_rel_prop"

##############################################################################
# Repeat with singletons
##############################################################################

summfile <- paste0(parentdir, "/output/7bp_1000k/full_j.summary")
binfile <- paste0(parentdir, "/output/7bp_1000k/full_bin.txt")

summ_5bp_100k <- read.table(summfile, header=F, stringsAsFactors=F, skip=1)
summ_5bp_100k <- summ_5bp_100k[sample(nrow(summ_5bp_100k), 12088037),]
names(summ_5bp_100k)<-c(
	"CHR", "POS", "REF", "ALT", "DP", "AN", "SEQ", "ALTSEQ", "GC")

summ_5bp_100k$BIN <- ceiling(summ_5bp_100k$POS/binw)

bins_5bp_100k <- read.table(binfile, header=T, stringsAsFactors=F, check.names=F)

dat_5bp_100k <- updateData(summ_5bp_100k, bins_5bp_100k, adj)
rm(summ_5bp_100k)
rm(bins_5bp_100k)

aggV <- aggData(dat_5bp_100k, adj) #<-modify the adj value for 3bp data

agg_5bp_100k <- aggV$oe
rates5 <- aggV$agg
summagg2 <- aggV$summagg2

##############################################################################
# Merge data
##############################################################################
r5m <- merge(rates5[,c(2,3,4,7,8)],
	rates5_common[,c(2,3,4,7,8)], by=c("Sequence", "Category2")) %>%
	mutate(Category=gsub("cpg_", "", Category2)) %>%
	#group_by(Category) %>%
	mutate(nsing=sum(num.x),
		ncommon=sum(num.y),
		prop_s=num.x/nsing,
		prop_p=num.y/ncommon,
		prop_diff=prop_p/prop_s)

r5m <- r5m %>%
	mutate(Category2 = plyr::mapvalues(Category2, orderedcats1, orderedcats2))

r5m$Category2 <- factor(r5m$Category2, levels=orderedcats2)

r5m$prop_diff3 <- r5m$prop_diff
r5m$prop_diff3[r5m$prop_diff< 0.5] <- 0.5
r5m$prop_diff3[r5m$prop_diff>2] <- 2

r5m$prop_diff4 <- r5m$num.y/r5m$num.x
r5m$prop_diff4 <- r5m$prop_diff4/mean(r5m$prop_diff4)

r5m$prop_diff5 <- r5m$prop_diff4
r5m$prop_diff4[r5m$prop_diff4< 0.5] <- 0.5
r5m$prop_diff4[r5m$prop_diff4>2] <- 2

r5m$v2 <- substr(r5m$Sequence,1,adj)
r5m$v2a <- as.character(lapply(as.vector(r5m$v2), reverse_chars))
r5m$v2a <- factor(r5m$v2a)
r5m$v3 <- substr(r5m$Sequence,adj+2,adj*2+1)

##############################################################################
# Compare subtypes
##############################################################################
# r5m %>% group_by(Category2) %>% summarise(cor=cor(rel_prop, common_rel_prop))

ggplot(data=r5m, aes(x=rel_prop, y=common_rel_prop))+
	geom_point(aes(group=Category2, colour=Category2), alpha=0.3, size=2, shape=1)+
  xlab("Relative mutation rate (ERVs)")+
  ylab("Relative mutation rate (MAC10+)")+
	geom_abline(intercept=0, linetype="dashed")+
	scale_x_log10(labels=c(0.0001, 0.001, 0.01, 0.1),
		breaks=c(0.0001, 0.001, 0.01, 0.1),
		limits=c(0.000025, 0.1))+
  scale_y_log10(labels=c(0.0001, 0.001, 0.01, 0.1),
		breaks=c(0.0001, 0.001, 0.01, 0.1),
		limits=c(0.000025, 0.1))+
	coord_fixed()+
	scale_colour_manual(values=cols)+
	guides(colour = guide_legend(title=NULL, nrow=3, override.aes = list(alpha=1, shape=16)))+
  theme_bw()+
  theme(
		legend.position="bottom",
    strip.text.x=element_text(size=14),
    legend.title=element_text(size=12),
    axis.title.x=element_text(size=14),
    axis.title.y=element_text(size=14),
		axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1),
		axis.text.y=element_text(size=12))
ggsave(paste0(parentdir, "/images/sing_com_cor_rates_no_facet.png"),
	width=4, height=5)

r5m2 <- r5m %>% dplyr::select(-Category2)

r5mcor <- r5m %>%
	group_by(Category2) %>%
	summarise(cor=cor(rel_prop, common_rel_prop, method="spearman")) %>%
	mutate(label=paste0("rho=", round(cor,2)))

ggplot(data=r5m, aes(x=rel_prop,  y=common_rel_prop))+
  geom_point(data=r5m2, alpha=0.1, size=2, colour="grey70")+
  geom_hex(bins=50)+
	geom_text(data=r5mcor, aes(label=paste0("rho==", round(cor,2))),
            x=-Inf, y=Inf, hjust=0, vjust=1, parse=TRUE, size=4)+
	xlab("Relative mutation rate (ERVs)")+
	ylab("Relative mutation rate (MAC10+)")+
  # geom_smooth(se=F, colour="red", linetype="dashed", aes(alpha=0.6))+
  geom_abline(intercept=0, linetype="dashed")+
  scale_x_log10(labels=c(0.0001, 0.001, 0.01, 0.1),
		breaks=c(0.0001, 0.001, 0.01, 0.1),
		limits=c(0.000025, 0.1))+
  scale_y_log10(labels=c(0.0001, 0.001, 0.01, 0.1),
		breaks=c(0.0001, 0.001, 0.01, 0.1),
		limits=c(0.000025, 0.1))+
	coord_fixed()+
  scale_colour_manual("Mutation Type", values=cols)+
  scale_fill_gradientn(colours = myPalette(6))+
	facet_wrap(~Category2, dir="v")+
  theme_bw()+
  theme(
    legend.position="none",
    strip.text.x=element_text(size=14),
    legend.title=element_text(size=12),
    axis.title.x=element_text(size=14),
    axis.title.y=element_text(size=14),
		axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1),
    axis.text.y=element_text(size=12))
ggsave(paste0(parentdir, "/images/sing_com_cor_rates_facet_dens.png"),
	width=6, height=6)

##############################################################################
# Plot 1d densities
##############################################################################
r5ma<-r5m %>%
	dplyr::select(Sequence, Category2, rel_prop, common_rel_prop) %>%
	gather(gp, prop, rel_prop:common_rel_prop) %>%
	mutate(group=recode(gp, rel_prop = "ERVs", common_rel_prop="MAC10+"))

kspvals <- r5m %>%
	group_by(Category2) %>%
	summarise(pval=ks.test(rel_prop, common_rel_prop)$p.value) #%>%
	# mutate()

kspvals$pv2<-ifelse(kspvals$pval==0, "2.2e-16", signif(kspvals$pval,2))

ggplot()+
	geom_density(data=r5ma, aes(x=log10(prop), group=group, colour=group))+
	geom_text(data=kspvals, aes(label=paste0("P<",pv2)),
						x=-Inf, y=Inf, hjust=0, vjust=1, parse=TRUE, size=4)+
	xlab("log10(relative mutation rate)")+
  facet_wrap(~Category2, scales="free", dir="v")+
  theme_bw()+
  theme(
    strip.text.x=element_text(size=14),
    legend.title=element_text(size=16),
    axis.title.x=element_text(size=18),
    axis.title.y=element_text(size=18))
ggsave(paste0(parentdir, "/images/sing_com_cor_rates_facet_densl.png"),
	width=9, height=9)

##############################################################################
# Compare subtypes (restrict to >50 observations)
##############################################################################
r5ma <- r5m %>% filter(num.x>=50 & num.y>=50)

ggplot(data=r5ma, aes(x=rel_prop, y=common_rel_prop))+
	geom_point(aes(group=Category2, colour=Category2),
		alpha=0.3, size=2, shape=1)+
  xlab("Relative mutation rate (ERVs)")+
  ylab("Relative mutation rate (MAC10+)")+
	geom_abline(intercept=0, linetype="dashed")+
	scale_x_log10(labels=c(0.0001, 0.001, 0.01, 0.1),
		breaks=c(0.0001, 0.001, 0.01, 0.1),
		limits=c(0.000025, 0.1))+
  scale_y_log10(labels=c(0.0001, 0.001, 0.01, 0.1),
		breaks=c(0.0001, 0.001, 0.01, 0.1),
		limits=c(0.000025, 0.1))+
	coord_fixed()+
	scale_colour_manual(values=cols)+
	guides(colour = guide_legend(title=NULL, nrow=3,
		override.aes = list(alpha=1, shape=16)))+
  theme_bw()+
  theme(
		legend.position="bottom",
    strip.text.x=element_text(size=14),
    legend.title=element_text(size=12),
    # axis.title.x=element_text(size=14),
    # axis.title.y=element_text(size=14),
		axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1),
		axis.text.y=element_text(size=12))
ggsave(paste0(parentdir, "/images/sing_com_cor_rates_no_facet_sub.png"),
	width=4, height=5)

r5m2a <- r5ma %>% dplyr::select(-Category2)

r5mcora <- r5ma %>%
	group_by(Category2) %>%
	summarise(cor=cor(rel_prop, common_rel_prop, method="spearman")) %>%
	mutate(label=paste0("rho=", round(cor,2)))

ggplot(data=r5ma, aes(x=rel_prop,  y=common_rel_prop))+
  geom_point(data=r5m2a, alpha=0.1, size=2, colour="grey70")+
  geom_hex(bins=50)+
	geom_text(data=r5mcora, aes(label=paste0("rho==", round(cor,2))),
            x=-Inf, y=Inf, hjust=0, vjust=1, parse=TRUE, size=4)+
	xlab("Relative mutation rate (ERVs)")+
	ylab("Relative mutation rate (MAC10+)")+
  # geom_smooth(se=F, colour="red", linetype="dashed", aes(alpha=0.6))+
  geom_abline(intercept=0, linetype="dashed")+
	scale_x_log10(labels=c(0.0001, 0.001, 0.01, 0.1),
		breaks=c(0.0001, 0.001, 0.01, 0.1),
		limits=c(0.000025, 0.1))+
  scale_y_log10(labels=c(0.0001, 0.001, 0.01, 0.1),
		breaks=c(0.0001, 0.001, 0.01, 0.1),
		limits=c(0.000025, 0.1))+
	coord_fixed()+
  scale_colour_manual("Mutation Type", values=cols)+
  scale_fill_gradientn(colours = myPalette(6))+
	facet_wrap(~Category2, dir="v")+
  theme_bw()+
  theme(
    legend.position="none",
    strip.text.x=element_text(size=14),
    legend.title=element_text(size=12),
    axis.title.x=element_text(size=14),
    axis.title.y=element_text(size=14),
		# axis.text.x=element_text(size=12),
    # axis.text.y=element_text(size=12),
		axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1),
		axis.text.y=element_text(size=12))
ggsave(paste0(parentdir, "/images/sing_com_cor_rates_facet_dens_sub.png"),
	width=6, height=6)

for(j in 1:6){
	categ<-orderedcats[j]
	p1a<-rrheat3(r5ma[r5ma$Category==categ,])
	 png(paste0(parentdir, "/images/rare_common_diff2", categ, "_sub.png"),
	 	height=5, width=5, units="in", res=300)
	 pushViewport(viewport(width=unit(5, "in"), height=unit(5, "in")))
	 grid.draw(ggplotGrob(p1a))
	 dev.off()
}

# trim whitespace on panels with imagemagick mogrify
trimcmd <- paste0("mogrify -trim ", parentdir, "/images/rare_common_diff*_sub.png")
system(trimcmd)

##############################################################################
# Compare MAC10+ with A&V rates
##############################################################################

avrates <- read.table(paste0(parentdir, "/posterior_7bp.txt"),
	header=T, stringsAsFactors=F, sep="\t")
head(avrates)
names(avrates)[1] <- "SEQUENCE"
avrates$CAT <- paste0(substr(avrates$SEQUENCE,4,4), substr(avrates$alt, 4,4))
avrates$Category[avrates$CAT=="AC" | avrates$CAT=="TG"] <- "AT_CG"
avrates$Category[avrates$CAT=="AG" | avrates$CAT=="TC"] <- "AT_GC"
avrates$Category[avrates$CAT=="AT" | avrates$CAT=="TA"] <- "AT_TA"
avrates$Category[avrates$CAT=="GA" | avrates$CAT=="CT"] <- "GC_AT"
avrates$Category[avrates$CAT=="GC" | avrates$CAT=="CG"] <- "GC_CG"
avrates$Category[avrates$CAT=="GT" | avrates$CAT=="CA"] <- "GC_TA"
cat("Assigning avrates CpG categories...\n")
# Second category column to include +3 CpG categories
# avrates$Category <- ifelse(substr(avrates$SEQUENCE,4,5)=="CG",
#                              paste0("cpg_",avrates$Category),
#                              avrates$Category)
names(avrates)[1] <- "Sequence"
# avrates$Sequence<-paste0(avrates$Sequence, "(", avrates$alt, ")")
r5m1<-r5m
r5m1$Sequence<-substr(r5m1$Sequence, 1, 7)
rm<-merge(r5m1, avrates, by=c("Sequence", "Category"))

rm<-rm %>%
	mutate(Category2 = plyr::mapvalues(Category2, orderedcats1, orderedcats2))
rm$Category2 <- factor(rm$Category2, levels=orderedcats2)
rm2 <- rm %>% dplyr::select(-Category2)

ggplot(data=rm, aes(x=common_rel_prop, y=eur/(mean(eur)/mean(common_rel_prop))))+
	geom_point(aes(group=Category2, colour=Category2),
		alpha=0.3, size=2, shape=1)+
  xlab("Relative mutation rate \n (BRIDGES MAC10+)")+
  ylab("Relative mutation rate \n (1000G EUR intergenic variants)")+
	geom_abline(intercept=0, linetype="dashed")+
	scale_x_log10(labels=c(0.0001, 0.001, 0.01, 0.1),
		breaks=c(0.0001, 0.001, 0.01, 0.1),
		limits=c(0.000025, 0.1))+
  scale_y_log10(labels=c(0.0001, 0.001, 0.01, 0.1),
		breaks=c(0.0001, 0.001, 0.01, 0.1),
		limits=c(0.000025, 0.1))+
	coord_fixed()+
	scale_colour_manual(values=cols)+
	guides(colour = guide_legend(title=NULL, nrow=3,
		override.aes = list(alpha=1, shape=16)))+
  theme_bw()+
  theme(
		legend.position="bottom",
    strip.text.x=element_text(size=14),
    legend.title=element_text(size=12),
    axis.title.x=element_text(size=14),
    axis.title.y=element_text(size=14),
		axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1),
		axis.text.y=element_text(size=12))
ggsave(paste0(parentdir, "/images/com_av_cor_rates_no_facet.png"),
	width=4, height=5)

rmcora <- rm %>%
	group_by(Category2) %>%
	summarise(cor=cor(common_rel_prop, eur, method="spearman")) %>%
	mutate(label=paste0("rho=", round(cor,2)))

ggplot(data=rm, aes(x=common_rel_prop, y=eur/(mean(eur)/mean(common_rel_prop))))+
  geom_point(data=rm2, alpha=0.3, size=2, colour="grey70")+
	# geom_point(aes(group=Category2, colour=Category2), alpha=0.3, size=2)+
  xlab("Relative rate (BRIDGES MAC10+)")+
  ylab("Relative mutation rate \n (1000G EUR intergenic variants)")+
	geom_hex(bins=50)+
	geom_text(data=r5mcora, aes(label=paste0("rho==", round(cor,2))),
            x=-Inf, y=Inf, hjust=0, vjust=1, parse=TRUE, size=4)+
  # geom_smooth(se=F, colour="red", linetype="dashed", aes(alpha=0.6))+
  geom_abline(intercept=0, linetype="dashed")+
	scale_x_log10(labels=c(0.0001, 0.001, 0.01, 0.1),
		breaks=c(0.0001, 0.001, 0.01, 0.1),
		limits=c(0.000025, 0.1))+
  scale_y_log10(labels=c(0.0001, 0.001, 0.01, 0.1),
		breaks=c(0.0001, 0.001, 0.01, 0.1),
		limits=c(0.000025, 0.1))+
	coord_fixed()+
  scale_colour_manual("Mutation Type", values=cols)+
  scale_fill_gradientn(colours = myPalette(6))+
	facet_wrap(~Category2, dir="v")+
  theme_bw()+
  theme(
    legend.position="none",
    strip.text.x=element_text(size=14),
    legend.title=element_text(size=12),
    axis.title.x=element_text(size=14),
    axis.title.y=element_text(size=14),
		axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1),
		axis.text.y=element_text(size=12))
ggsave(paste0(parentdir, "/images/com_av_cor_rates_facet.png"),
  width=6, height=6)

rm$prop_diff <- (rm$eur/(mean(rm$eur)/mean(rm$common_rel_prop)))/rm$common_rel_prop
rm$prop_diff4 <- rm$prop_diff
rm$prop_diff4[rm$prop_diff< 0.5] <- 0.5
rm$prop_diff4[rm$prop_diff>2] <- 2

# rm$prop_diff4 <- rm$num.y/rm$num.x
# rm$prop_diff4 <- rm$prop_diff4/mean(rm$prop_diff4)
#
# rm$prop_diff5 <- rm$prop_diff4
# rm$prop_diff4[rm$prop_diff4< 0.5] <- 0.5
# rm$prop_diff4[rm$prop_diff4>2] <- 2

# r5m$v2 <- substr(r5m$Sequence,1,adj)
# r5m$v2a <- as.character(lapply(as.vector(r5m$v2), reverse_chars))
# r5m$v2a <- factor(r5m$v2a)
# r5m$v3 <- substr(r5m$Sequence,adj+2,adj*2+1)

for(j in 1:6){
	categ<-orderedcats[j]
	dat <- rm[rm$Category==categ,]
	dat$prop_diff <- (dat$eur/dat$common_rel_prop)*(mean(dat$common_rel_prop)/mean(dat$eur))
	dat$prop_diff4 <- dat$prop_diff
	dat$prop_diff4[dat$prop_diff< 0.5] <- 0.5
	dat$prop_diff4[dat$prop_diff>2] <- 2
	print(nrow(dat[(dat$prop_diff4==2 | dat$prop_diff4==0.5),]))
	p1a<-rrheat3(dat)
	 png(paste0(parentdir, "/images/rare_common_diff2", categ, "_com_av.png"),
	 	height=5, width=5, units="in", res=300)
	 pushViewport(viewport(width=unit(5, "in"), height=unit(5, "in")))
	 grid.draw(ggplotGrob(p1a))
	 dev.off()
}

# trim whitespace on panels with imagemagick mogrify
trimcmd <- paste0("mogrify -trim ", parentdir, "/images/rare_common_diff*_av.png")
system(trimcmd)

##############################################################################
# Plot A>G subtypes
##############################################################################
pcadat <- prcomp(r5m[r5m$Category2=="A>G",c(5,8)])

dat <- data.frame(cbind(pcadat$x, r5m[r5m$Category2=="A>G",]))
dat$cluster <- ifelse(dat$PC2>0.003, 2,1)
dat$Sequence <- substr(dat$Sequence, 0, 7)
dat$cgs <- ifelse(substr(dat$Sequence, 1,2)%in%c("CC", "CG", "GG", "GC"), 1,0)
dat$ts <- ifelse(grepl("T", dat$Sequence), 0,1)
dat$num_GC_flank <- ifelse(str_count(dat$Sequence, "A|T")<=3,">=4","<4")
dat$maxc <- ifelse(log(dat$prop_diff5)>0.69, ">2", "<2")
dat$maxc <- factor(dat$maxc, levels=c("<2", ">2"))

ggplot(data=dat,
		aes(x=rel_prop, y=common_rel_prop,
			group=Category2, colour=maxc, shape=num_GC_flank))+
  geom_point(alpha=0.7, size=2)+
	scale_colour_brewer(palette="Set2")+
	scale_x_log10()+
	scale_y_log10()+
	scale_shape(solid = FALSE)+
  xlab("Relative rate (downsampled BRIDGES ERVs)")+
  ylab("Relative rate (BRIDGES MAC10+ variants)")+
	guides(shape = guide_legend("#G/C bases in flanking region"),
		colour = guide_legend("MAC10+:ERV ratio", override.aes = list(alpha=1)))+
  theme_bw()+
  theme(legend.position = "bottom",
    legend.title=element_text(size=16),
    axis.title.x=element_text(size=18),
    axis.title.y=element_text(size=18))
ggsave(paste0(parentdir, "/images/AT_GC_scatter.png"), width=8, height=6)
