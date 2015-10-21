##############################################################################
# Negative binomial regression models and plots
##############################################################################

# First plot paneled histograms of each mutation category
agg_oe <- gather(agg_5bp_100k, var, value, c(obs,exp))

ggplot(agg_oe, aes(x=value, fill=var))+
	geom_histogram(alpha=0.5, position="identity")+
	facet_wrap(~Category2, scales="free")+
	theme_bw()

cathistfile <- paste0(parentdir, "/images/obs_hist.png")
ggsave(cathistfile)

##############################################################################
# For each motif, calculates correlation between observed and mutable sites,
# then plots this correlation against the motif frequency
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
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
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

agg_cov <- merge(agg_5bp_100k, mut_cov, by=c("CHR", "BIN"))

agg_cov$ratio <- agg_cov$exp/agg_cov$obs
agg_cov <- filter(agg_cov, ratio<5) %>% mutate(med=(maxn+minn)/2)


##############################################################################
# Plots %GC against variance of motifs in each bin
##############################################################################
ggplot(agg_cov, aes(x=prop_GC, y=motifvar))+
	geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
	# geom_errorbar(aes(ymax=maxn, ymin=minn))+
	facet_wrap(~Category2, scales="free")+
	theme_bw()

gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
ggsave(gcfile)

##############################################################################
# Run negbin models
##############################################################################

# Get 1bp relative rates
sa2<-summagg2 %>%
	group_by(CHR, BIN, Category2) %>%
	summarise(ct=sum(Count, na.rm=T), obs=sum(obs, na.rm=T))
sa3<-sa2 %>%
	group_by(Category2) %>%
	summarise(ct=sum(ct), obs=sum(obs), rel_prop=obs/ct)

a3<-merge(agg_cov, sa3[,c(1,2,4)], by="Category2")
a3$exp1<-a3$nmotifs*a3$rel_prop

# plotdat<-filter(a3, Category2=="GC_TA" & EXON<0.1)
plotdat<-filter(a3, Category2=="GC_TA")
ggplot(plotdat, aes(x=prop_GC, y=obs))+
	geom_point(alpha=0.3)+
	facet_wrap(~Category2, scales="free")+
	theme_bw()

gcobsfile <- paste0(parentdir, "/images/gc_vs_obs.png")
ggsave(gcobsfile, width=12, heigh=12)

mc2 <- a3 %>%
	group_by(Category2) %>%
	summarise(cor1=cor(exp, obs, method="pearson"),
		cor2=cor(exp1, obs, method="pearson"))

# Initialize data of motif counts to use directly as covariates
# in negbin model
# dat_5bp_100k$bin$CHR <- as.integer(substring(dat_5bp_100k$bin$CHR, 4))

names(dat_5bp_100k$bin) <- gsub('\\(', '_', names(dat_5bp_100k$bin))
names(dat_5bp_100k$bin) <- gsub('\\)', '_', names(dat_5bp_100k$bin))

atcols <- c(names(dat_5bp_100k$bin)[1:5],
  names(dat_5bp_100k$bin)[which(substr(names(dat_5bp_100k$bin), 3, 3)=="A")])

binsAT <- dat_5bp_100k$bin %>%
  select_(.dots = atcols) %>%
	arrange(CHR, BIN)

gcdn <- c("CA", "CC", "CT")

gccols <- c(names(dat_5bp_100k$bin)[1:5],
  names(dat_5bp_100k$bin)[which(substr(names(dat_5bp_100k$bin), 3, 4) %in% gcdn)])

binsGC <- dat_5bp_100k$bin %>%
  select_(.dots = gccols) %>%
	arrange(CHR, BIN)

cpggccols <- c(names(dat_5bp_100k$bin)[1:5],
  names(dat_5bp_100k$bin)[which(substr(names(dat_5bp_100k$bin), 3, 4)=="CG")])

binscpgGC <- dat_5bp_100k$bin %>%
  select_(.dots = cpggccols) %>%
	arrange(CHR, BIN)

# Motifs + genomic features
compare.all <- data.frame()
compare.err <- data.frame()
logliks <- data.frame()
mut_cats <- unique(agg_5bp_100k$Category2)

for(i in 1:length(mut_cats)){
	cat1 <- mut_cats[i]
	aggcat <- a3[a3$Category2==mut_cats[i],]

	if(grepl("^AT", cat1)){
		aggcatm <- merge(aggcat, binsAT, by=c("CHR", "BIN"))
		mcols <- atcols
	} else if (grepl("^GC", cat1)){
		aggcatm <- merge(aggcat, binsGC, by=c("CHR", "BIN"))
		mcols <- gccols
	} else{
		aggcatm <- merge(aggcat, binscpgGC, by=c("CHR", "BIN"))
		mcols <- cpggccols
	}

	aggcat$exp_s <- sample(aggcat$obs, length(aggcat$obs), replace=T)
	aggcatm$exp_s <- sample(aggcatm$obs, length(aggcatm$obs), replace=T)

	# Fit models with all data
	feat_mod_formula <- as.formula(paste("obs~",
		paste(covnames, collapse="+")))
	full_mod_formula <- as.formula(paste("obs~",
		paste(mcols[-(1:5)], collapse="+"), "+",
		paste(covnames, collapse="+")))
	motif_mod_formula <- as.formula(paste("obs~",
		paste(mcols[-(1:5)], collapse="+")))
	mut_lm_gc <- glm.nb(obs~prop_GC, data=aggcat)
	mut_lm_feat <- glm.nb(feat_mod_formula, data=aggcat)
	mut_lm_motif <- glm.nb(obs~exp, data=aggcat)
	mut_lm_motif2 <- glm.nb(motif_mod_formula, data=aggcatm)
	mut_lm_1bp <- glm.nb(obs~exp1, data=aggcat)
	# mut_lm_3bp <- glm.nb(obs~exp3, data=aggcat)
	mut_lm_full <- glm.nb(full_mod_formula, data=aggcatm)

	ll_5bp <- data.frame(
		Category2=cat1,
		params=256,
		ll=summary(mut_lm_motif)$twologlik)
	ll_1bp <- data.frame(
		Category2=cat1,
		params=6,
		ll=summary(mut_lm_1bp)$twologlik)

	logliks <- rbind(logliks, ll_5bp, ll_1bp)

	# 10-fold cross-validation--may need to update so expected counts are
	# re-calculated for each 1/N subset
	gc_cv <- cv.glm(data=aggcat, glmfit=mut_lm_gc, K=10)
	feat_cv <- cv.glm(data=aggcat, glmfit=mut_lm_feat, K=10)
	motif_cv <- cv.glm(data=aggcat, glmfit=mut_lm_motif, K=10)
	full_cv <- cv.glm(data=aggcatm, glmfit=mut_lm_full, K=10)

	mspe <- c(gc_cv$delta[2], feat_cv$delta[2], motif_cv$delta[2], full_cv$delta[2])
	rmse <- sqrt(mspe)
	meanct <- mean(aggcat$obs)
	pcterr <- rmse/meanct
	mspe.res <- c("GC", "features", "motifs", "motifs+features")
	mspe.dat <- data.frame(Category2=cat1, res=mspe.res, mspe, rmse, meanct, pcterr)
	compare.err <-rbind(compare.err, mspe.dat)

	fits_gc <- mut_lm_gc$fitted.values
	names(fits_gc) <- paste0(aggcat$CHR,".",aggcat$BIN)

	fits_feat <- mut_lm_feat$fitted.values
	names(fits_feat) <- paste0(aggcat$CHR,".",aggcat$BIN)

	fits_motif <- mut_lm_motif$fitted.values
	names(fits_motif) <- paste0(aggcat$CHR,".",aggcat$BIN)

	fits_motif2 <- mut_lm_motif2$fitted.values
	names(fits_motif2) <- paste0(aggcatm$CHR,".",aggcatm$BIN)

	fits_full <- mut_lm_full$fitted.values
	names(fits_full) <- paste0(aggcatm$CHR,".",aggcatm$BIN)

	BIN <- as.integer(gsub(".*\\.", "", names(fits_feat)))
	CHR <- as.integer(gsub("\\..*", "", names(fits_feat)))

	BINA <- as.integer(gsub(".*\\.", "", names(fits_motif2)))
	CHRA <- as.integer(gsub("\\..*", "", names(fits_motif2)))

	model_dat_gc <- data.frame(CHR,Category2=cat1, BIN,
		exp=fits_gc,
		obs=aggcat$obs,
		exp_s=aggcat$exp_s,
		stringsAsFactors=F)
	model_dat_gc$res <- "GC"

	model_dat_feat <- data.frame(CHR,Category2=cat1, BIN,
		exp=fits_feat,
		obs=aggcat$obs,
		exp_s=aggcat$exp_s,
		stringsAsFactors=F)
	model_dat_feat$res <- "features"

	model_dat_motif <- data.frame(CHR,Category2=cat1, BIN,
		exp=fits_motif,
		# exp=aggcat$exp,
		obs=aggcat$obs,
		exp_s=aggcat$exp_s,
		stringsAsFactors=F)
	model_dat_motif$res <- "motifs"

	model_dat_motif2 <- data.frame(CHR=CHRA,Category2=cat1, BIN=BINA,
		exp=fits_motif2,
		# exp=aggcat$exp,
		obs=aggcatm$obs,
		exp_s=aggcatm$exp_s,
		stringsAsFactors=F)
	model_dat_motif2$res <- "motifs2"

	model_dat_full <- data.frame(CHR=CHRA,Category2=cat1, BIN=BINA,
		exp=fits_full,
		obs=aggcat$obs,
		exp_s=aggcat$exp_s,
		stringsAsFactors=F)
	model_dat_full$res <- "motifs+features"

	compare.all <- rbind(compare.all,
		model_dat_gc,
		model_dat_feat,
		model_dat_motif,
		model_dat_motif2,
		model_dat_full)

	# z <- summary(mut_lm_full)
	# print(cat1)
	# print(z)
}

ggplot(logliks, aes(x=params, y=ll, group=Category2, colour=Category2))+
	geom_line()+
	theme_bw()+
	xlab("# parameters")+
	ylab("2 log likelihood")+
ggsave("/net/bipolar/jedidiah/mutation/images/ll.png")

# Update comparison data frame
# compare.all<-rbind(agg_5bp_100k[,c(1,2,3,4,5,8)], d, d1)
#^ this version uses direct estimates from motif rates; not fair comparison
compare.all$res <- factor(compare.all$res,
	levels = c("GC", "features", "motifs", "motifs2", "motifs+features"))
compare.all$diff <- compare.all$obs-compare.all$exp
compare.all$diff_s <- compare.all$obs-compare.all$exp_s
compare.all$Category2 <- factor(compare.all$Category2)
# levels(compare.all$Category2) <- c(
	# "AT>CG", "AT>GC", "AT>TA",
	# "(CpG)GC>AT", "(CpG)GC>CG", "(CpG)GC>TA",
	# "GC>AT", "GC>CG", "GC>TA")

ca.gp <- group_by(compare.all, Category2, res)
compare.all <- mutate(ca.gp, dm=mean(diff), dsd=sd(diff), zscore=(diff-dm)/dsd)
compare.all$err <- abs(compare.all$diff)/compare.all$obs
compare.all$err_s <- abs(compare.all$diff_s)/compare.all$obs

# get mean absolute error for each category/model
ca <- compare.all %>%
	group_by(Category2, res) %>%
	summarise(meanerr_s=mean(err_s),
		meanerr=mean(err),
		sderr=std(err),
		meanobs=mean(obs),
		sdobs=std(obs),
		fold=meanerr_s/meanerr)

# Get mean fold-improvement of model error vs. error from random draws from dist
ca.fold <- ca %>% group_by(res) %>% summarise(mean=mean(fold))

caout <- paste0(parentdir, "/output/", nbp, "bp_err.txt")
write.table(ca, caout, col.names=T, row.names=F, quote=F, sep="\t")
cad <- dcast(data.frame(ca), Category2~res, value.var="meanerr")

# ac2 <- merge(compare.all, mut_cov, by=c("CHR", "BIN")) %>%
	# group_by("Category2") %>%
	# summarise(corgc=cor(diff, prop_GC))


# Function to compute standard error for correlations
corSE<-function(corval, ct){
	sqrt((1-corval^2)/(ct-2))
}

# plot raw 5bp motif predictions vs observed
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

a2 <- mutate(agg_5bp_100k, rel=obs/nmotifs)
ggplot(a2, aes(x=rel, y=nmotifs))+
	geom_point()+
	facet_wrap(~Category2, scales="free")
ggsave("/net/bipolar/jedidiah/mutation/images/rel_ct_cor.png")

a3 <- a2 %>% group_by(Category2) %>% summarise(cor=cor(rel, nmotifs))

# plot negbin model 5bp motif predictions vs observed
plotdat <- compare.all[compare.all$res=="motifs",]
p2 <- ggplot(plotdat, aes(x=obs, y=exp, colour=res))+
	geom_point(alpha=0.2, size=3)+
	scale_colour_manual("Model", values=myPaletteCat(8)[7])+
	facet_wrap(~Category2, ncol=3, scales="free")+
	ylab("Predicted count")+
	xlab("Observed count")+
	theme_bw()+
	theme(axis.title.y=element_text(size=16),
		axis.text.y=element_text(size=14),
		axis.title.x=element_text(size=16),
		axis.text.x=element_text(size=14),
		legend.title=element_text(size=16),
		legend.text=element_text(size=14))

hierfile7 <- paste0(parentdir, "/images/negbin_motif_pred_vs_obs.png")
ggsave(hierfile7, width=18, height=18)

# plot negbin model 5bp motif predictions vs observed
p2 <- ggplot(compare.all, aes(x=obs, y=exp, colour=res))+
	geom_point(alpha=0.2, size=3)+
	geom_point(alpha=0.2, size=3, data=filter(plotdat, res=="features"))+
	geom_point(alpha=0.2, size=3, data=filter(plotdat, res=="motifs+features"))+
	scale_colour_manual("Model", values=myPaletteCat(8)[5:8])+
	facet_wrap(~Category2, ncol=3, scales="free")+
	ylab("Predicted count")+
	xlab("Observed count")+
	theme_bw()+
	theme(axis.title.y=element_text(size=16),
		axis.text.y=element_text(size=14),
		axis.title.x=element_text(size=16),
		axis.text.x=element_text(size=14),
		legend.title=element_text(size=16),
		legend.text=element_text(size=14))

hierfile4 <- paste0(parentdir, "/images/negbin_all_pred_vs_obs.png")
ggsave(hierfile4, width=18, height=18)

# Plot scatterplot across chr of model errors (using negbin motif predictions)
plotdat <- compare.all[compare.all$CHR==2 & compare.all$diff>-150,]
p2 <- ggplot(plotdat, aes(x=BIN, y=zscore, colour=res))+
	geom_point(alpha=0.4, size=4)+
	geom_point(alpha=0.4, size=4, data=filter(plotdat, res=="features"))+
	geom_point(alpha=0.4, size=4, data=filter(plotdat, res=="motifs+features"))+
	scale_colour_manual("Model", values=myPaletteCat(8)[4:8])+
	facet_wrap(~Category2, scales="free", ncol=3)+
	ylab("Error")+
	xlab(NULL)+
	theme_bw()+
	theme(axis.title.y=element_text(size=16),
		axis.text.y=element_text(size=14),
		legend.title=element_text(size=16),
		legend.text=element_text(size=14),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank())

hierfile5 <- paste0(parentdir, "/images/hier_diffs_pt5.png")
ggsave(hierfile5, width=18, height=18)

# Plot barcharts comparing obs/exp correlation for different models
mod.corr <- compare.all %>%
	filter(exp>100) %>%
	group_by(Category2, res) %>%
	summarise(num=length(exp),
		cor=cor(exp, obs, method="pearson"),
		cor.p=cor.test(exp, obs, method="pearson")$p.value)
mod.corr$SE <- corSE(mod.corr$cor, mod.corr$num)
limits <- aes(ymax = mod.corr$cor + mod.corr$SE,
	ymin=mod.corr$cor - mod.corr$SE)
dodge <- position_dodge(width=0.9)

ggplot(mod.corr, aes(x=Category2, y=cor, fill=res))+
	geom_bar(stat="identity", position=dodge)+
	scale_colour_manual("Predictor", values=myPaletteCat(8)[4:8])+
	scale_fill_manual("Predictor", values=myPaletteCat(8)[4:8])+
	xlab("Category")+
	ylab("Correlation with observed count")+
	geom_errorbar(limits, position=dodge, width=0.25)+
	theme_bw()+
	theme(legend.title = element_text(size=18),
		legend.text = element_text(size=16),
		axis.title.x = element_text(size=20),
		axis.title.y = element_text(size=20),
		axis.text.y = element_text(size=16),
		axis.text.x = element_text(size=16, angle = 45,  vjust=1, hjust=1.01))

modelbar <- paste0(parentdir, "/images/gw_5bp_vs_mod.png")
ggsave(modelbar, width=7, height=7)

# Plot barcharts comparing 10-fold cross validation MSPE
compare.err$res <- factor(compare.err$res,
	levels = c("GC", "features", "motifs", "motifs+features"))

ggplot(compare.err, aes(x=Category2, y=mspe, fill=res))+
	geom_bar(stat="identity", position=dodge)+
	scale_colour_manual("Predictor", values=myPaletteCat(8)[5:8])+
	scale_fill_manual("Predictor", values=myPaletteCat(8)[5:8])+
	facet_wrap(~Category2, scales="free")+
	ylab("MSPE")+
	theme_bw()+
	theme(legend.title = element_text(size=18),
		legend.text = element_text(size=16),
		axis.title.y = element_text(size=20),
		axis.title.x = element_blank(),
	  axis.text.y = element_text(size=16),
	  axis.text.x = element_blank(),
	  axis.ticks.x = element_blank())

mspebar <- paste0(parentdir, "/images/compare_mspe.png")
ggsave(mspebar, width=12, height=12)
