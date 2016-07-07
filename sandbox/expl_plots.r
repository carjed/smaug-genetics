##############################################################################
# Exploratory plots
##############################################################################

##############################################################################
# First plot paneled histograms of each mutation category
##############################################################################
agg_oe <- gather(agg_5bp_100k, var, value, c(obs,exp))

ggplot(agg_oe, aes(x=value, fill=var))+
	geom_histogram(alpha=0.5, position="identity")+
	facet_wrap(~Category2, scales="free")+
	theme_bw()

cathistfile <- paste0(parentdir, "/images/obs_hist.png")
ggsave(cathistfile)

##############################################################################
# plot raw 5bp motif predictions vs observed
##############################################################################
p2 <- ggplot(agg_5bp_100k, aes(x=obs, y=nmotifs))+
	geom_point(alpha=0.2, size=3)+
	scale_colour_manual("Model", values=myPaletteCat(8)[6:7])+
	facet_wrap(~Category2, scales="free", ncol=3)+
	ylab("Mutable motifs")+
	xlab("Observed count")+
	theme_bw()+
	theme(axis.title.y=element_text(size=16),
		axis.text.y=element_text(size=14),
		axis.title.x=element_text(size=16),
		axis.text.x=element_text(size=14),
		legend.title=element_text(size=16),
		legend.text=element_text(size=14))

hierfile8 <- paste0(parentdir, "/images/raw_motif_pred_vs_obs.png")
ggsave(hierfile8, width=18, height=18)

##############################################################################
# For each motif, calculates correlation between observed and mutable sites,
# then plot this correlation against the motif frequency
#
# In AT>NN and CpG GC>NN categories, appears to be a relationship where more
# common motifs are associated with negative or positive correlations,
# respectively, e.g., for AT>CG mutations, more common motifs have a strong
# negative correlation with
##############################################################################
agg_cov2 <- summagg2 %>%
	group_by(Category2, Sequence) %>%
		summarise(corgc=cor(obs, prop_GC, use="complete.obs"),
		cornm=cor(obs, Count, use="complete.obs"),
		nmotifs=sum(Count),
		rel_prop=mean(rel_prop, na.rm=T))

# agg_cov2 %>% summarise(cor=cor(cor, nmotifs, method="spearman"))

agg_cov2$mgc <- nchar(gsub("[AT]", "", agg_cov2$Sequence))-2

ggplot(agg_cov2, aes(x=corgc, y=cornm, colour=mgc, size=nmotifs))+
	geom_point(alpha=0.3)+
	scale_colour_continuous(name="Relative Mutation Rate")+
	facet_wrap(~Category2, scales="free")+
	xlab("# mutable sites~singleton correlation")+
	ylab("%GC~singleton correlation")+
	theme_bw()

gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)

##############################################################################
# Plots %GC against variance of motifs in each bin
##############################################################################
ggplot(agg_cov2, aes(x=nmotifs, y=corgc))+
	geom_point(alpha=0.3)+
	facet_wrap(~Category2, scales="free")+
	xlab("# mutable sites")+
	ylab("Singleton~GC correlation")+
	theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor2.png")
ggsave(gccorfile)

##############################################################################
# Plots %GC against observed mutations, colored by motif
##############################################################################
agg_cov2 <- merge(summagg2, mut_cov, by=c("CHR", "BIN"))

agg_cov2$diff <- agg_cov2$exp-agg_cov2$obs
agg_cov2$err <- agg_cov2$diff/agg_cov2$obs

err_motif <- agg_cov2 %>%
	group_by(Category, Sequence) %>%
	summarise(rel_prop=mean(rel_prop, na.rm=T),
		corgc=cor(err, prop_GC.x, use="complete.obs"),
		corrc=cor(err, RATE, use="complete.obs"),
		coroe=cor(obs, exp, use="complete.obs"),
		nmotifs=sum(Count),
		meanerr=mean(err, na.rm=T)) %>%
	arrange(Category, corgc) %>%
	mutate(rank=rank(corgc))

# Create separate categories of CpG-specific motifs
err_motif$Category2 <- ifelse(
	substr(err_motif$Sequence,adj+1,adj+2) == "CG",
	paste0("cpg_",err_motif$Category),
	err_motif$Category)

# plotdat<-err_motif[err_motif$Category2=="GC_AT",]
plotdat<-err_motif
ggplot(plotdat,
		aes(x=rank, y=corgc, group=Category, colour=Category2, size=log(rel_prop)))+
	geom_point(alpha=0.7)+
	scale_y_continuous(breaks=seq(-0.8, 0.8, 0.1))+
	scale_x_continuous(breaks=c(0,64,128,192,256))+
	# scale_size_continuous(breaks=seq(log(min(plotdat$rel_prop)),
		# log(max(plotdat$rel_prop)),by=1))+
	xlab("Motif")+
	ylab("Error~GC Correlation")+
	theme_bw()

gcerrfile <- paste0(parentdir, "/images/err-GC_cor.png")
ggsave(gcerrfile)

ggplot(plotdat,
		aes(x=Category2, y=meanerr, group=Category2, fill=Category2))+
	geom_boxplot()+
	# scale_y_continuous(breaks=seq(-0.1, 0.8, 0.1))+
	# scale_x_continuous(breaks=c(0,64,128,192,256))+
	# scale_size_continuous(breaks=seq(log(min(plotdat$rel_prop)),
		# log(max(plotdat$rel_prop)),by=1))+
	ylab("Mean % Error")+
	theme_bw()

gcboxfile <- paste0(parentdir, "/images/err-GC_box.png")
ggsave(gcboxfile)

# [err_motif$Category2!="cpg_GC_AT",]

##############################################################################
# Plots %GC against variance of motifs in each bin
##############################################################################
# ggplot(agg_cov, aes(x=prop_GC, y=motifvar))+
# 	geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
# 	facet_wrap(~Category2, scales="free")+
# 	theme_bw()
#
# gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
# ggsave(gcfile)

##############################################################################
# Plot GC content vs GC>TA singleton density to show negative correlation
# -must be done after a3 dataframe has been created
##############################################################################
# plotdat<-filter(a3, Category2=="GC_TA" & EXON<0.1)
# plotdat<-filter(a3, Category2=="GC_TA")
# ggplot(plotdat, aes(x=prop_GC, y=obs))+
# 	geom_point(alpha=0.3)+
# 	facet_wrap(~Category2, scales="free")+
# 	theme_bw()
#
# gcobsfile <- paste0(parentdir, "/images/gc_vs_obs.png")
# ggsave(gcobsfile, width=12, heigh=12)
#
# mc2 <- a3 %>%
# 	group_by(Category2) %>%
# 	summarise(cor1=cor(exp, obs, method="pearson"),
# 		cor2=cor(exp1, obs, method="pearson"))
