extra_plots <- 0
if(extra_plots){source("expl_plots.r")}

##############################################################################
# Merge aggregate data with covariate data
##############################################################################
agg_cov <- merge(agg_5bp_100k, mut_cov, by=c("CHR", "BIN"))
agg_cov$ratio <- agg_cov$exp/agg_cov$obs
agg_cov <- filter(agg_cov, ratio<5) %>% mutate(med=(maxn+minn)/2)

##############################################################################
# Get relative rates for 1bp and 3bp motifs and compare likelihoods
##############################################################################
rates5$Seq3 <- substr(rates5$Sequence, adj, nbp-adj+1)

rates1 <- rates5 %>%
  dplyr::select(Sequence, Category2, COUNT, num) %>%
  group_by(Category2) %>%
  summarise(COUNT1=sum(COUNT), num1=sum(num), rel_prop1=num1/COUNT1)

rates3 <- rates5 %>%
  dplyr::select(Seq3, Category2, COUNT, num) %>%
  group_by(Seq3, Category2) %>%
  summarise(COUNT3=sum(COUNT), num3=sum(num), rel_prop3=num3/COUNT3)

if(nbp==7){
	rates5$Seq5 <- substr(rates5$Sequence, adj-1, nbp-1)

	rates5a <- rates5 %>%
	  dplyr::select(Seq5, Category2, COUNT, num) %>%
	  group_by(Seq5, Category2) %>%
	  summarise(COUNT5=sum(COUNT), num5=sum(num), rel_prop5=num5/COUNT5)
}

rates_full <- merge(rates5, rates1, by="Category2")
rates_full <- merge(rates_full, rates3, by=c("Seq3", "Category2"))

if(nbp==7){
	rates_full <- merge(rates_full, rates5a, by=c("Seq5", "Category2"))
}

rates_full$logLik1 <- dbinom(rates_full$num,
  rates_full$COUNT,
  rates_full$rel_prop1, log=T)

rates_full$logLik3 <- dbinom(rates_full$num,
  rates_full$COUNT,
  rates_full$rel_prop3, log=T)

if(nbp==7){
	rates_full$logLik5 <- dbinom(rates_full$num,
	  rates_full$COUNT,
	  rates_full$rel_prop5, log=T)

	rates_full$logLik7 <- dbinom(rates_full$num,
	  rates_full$COUNT,
	  rates_full$rel_prop, log=T)

} else {
	rates_full$logLik5 <- dbinom(rates_full$num,
	  rates_full$COUNT,
	  rates_full$rel_prop, log=T)
}

if(nbp==7){
	ra1<-rates_full %>%
	  group_by(Category2) %>%
	  summarise(L1=-2*sum(logLik1),
	    L3=-2*sum(logLik3),
	    L5=-2*sum(logLik5),
			L7=-2*sum(logLik7))

	ra1a<-gather(ra1, model, log, L1:L7)

	ra1a$k <- c(rep(1,9),
	  c(rep(16,3), rep(12,3), rep(4,3)),
	  c(rep(256,3), rep(192,3), rep(64,3)),
		c(rep(4096,3), rep(3072,3), rep(1024,3)))

	ra1a$AIC <- 2*ra1a$k+2*ra1a$log

	ra1b <- gather(ra1a, stat, L, c(log, AIC))
	levels(ra1b$model) <- c("1", "3", "5", "7")

} else {
	ra1<-rates_full %>%
	  group_by(Category2) %>%
	  summarise(L1=-2*sum(logLik1),
	    L3=-2*sum(logLik3),
	    L5=-2*sum(logLik5))

	ra1a<-gather(ra1, model, log, L1:L5)

	ra1a$k <- c(rep(1,9),
	  c(rep(16,3), rep(12,3), rep(4,3)),
	  c(rep(256,3), rep(192,3), rep(64,3)))

	ra1a$AIC <- 2*ra1a$k+2*ra1a$log

	ra1b <- gather(ra1a, stat, L, c(log, AIC))
	levels(ra1b$model) <- c("1", "3", "5")
}

levels(ra1b$stat) <- c("-2ln(L)", "AIC")
names(ra1b) <- c("Category", "Motif_Length", "k", "Stat", "L")

# Plot AIC and -2log(L) for each category
ggplot(ra1b, aes(x=Motif_Length, y=L, group=Stat, colour=Stat))+
  scale_colour_brewer(palette="Set1")+
  geom_point()+
  geom_line()+
  facet_wrap(~Category, scales="free")+
  theme_bw()+
  xlab("Motif Length")+
  theme(axis.title.y=element_blank(),
    legend.title=element_blank())

ggsave("/net/bipolar/jedidiah/mutation/images/compare_AIC.png")

##############################################################################
# Initialize data of motif counts to use as covariates in negbin model
##############################################################################
dat_5bp_100k$bin$CHR <- as.integer(substring(dat_5bp_100k$bin$CHR, 4))

names(dat_5bp_100k$bin) <- gsub('\\(', '_', names(dat_5bp_100k$bin))
names(dat_5bp_100k$bin) <- gsub('\\)', '_', names(dat_5bp_100k$bin))

atcols <- c(names(dat_5bp_100k$bin)[1:5],
  names(dat_5bp_100k$bin)[which(substr(names(dat_5bp_100k$bin), adj+1, adj+1)=="A")])

binsAT <- dat_5bp_100k$bin %>%
  select_(.dots = atcols) %>%
	arrange(CHR, BIN)

gcdn <- c("CA", "CC", "CT")
gccols <- c(names(dat_5bp_100k$bin)[1:5],
  names(dat_5bp_100k$bin)[which(substr(names(dat_5bp_100k$bin), adj+1, adj+2) %in% gcdn)])

binsGC <- dat_5bp_100k$bin %>%
  select_(.dots = gccols) %>%
	arrange(CHR, BIN)

cpggccols <- c(names(dat_5bp_100k$bin)[1:5],
  names(dat_5bp_100k$bin)[which(substr(names(dat_5bp_100k$bin), adj+1, adj+2)=="CG")])

binscpgGC <- dat_5bp_100k$bin %>%
  select_(.dots = cpggccols) %>%
	arrange(CHR, BIN)

##############################################################################
# Run combined negbin model
##############################################################################

# Testing--use this to calculate weighted means for all 3 motif lengths
rates_full_s <- rates_full %>%
	dplyr::select(Sequence, Category2, rel_prop, rel_prop1, rel_prop3)

# Add expected #singletons per bin using 1bp motifs
# -used in comparing category-specific models with weighted mean method
a3 <- merge(agg_cov, rates1[,c(1,2,4)], by="Category2")
a3$exp1 <- a3$nmotifs*a3$rel_prop1

# Aggregate across categories for total #singletons per bin
a3a <- a3 %>%
  group_by_(.dots=lapply(names(a3)[c(2,3,10:22)], as.symbol)) %>%
  summarise(obs=sum(obs))

# Add motif counts per bin
a3a1 <- merge(a3a, dat_5bp_100k$bin, by=c("CHR", "BIN"))

motif_mod_form <- as.formula(paste("obs~",
	paste(names(a3a1)[-c(1:19)], collapse="+")))

feat_mod_form <- as.formula(paste("obs~",
  paste(names(a3a1)[3:15], collapse="+")))

full_mod_form <- as.formula(paste("obs~",
	paste(names(a3a1)[-c(1,2,16:19)], collapse="+")))

# Uses poisson regression instead of negbin, since negbin fails to converge
# with default parameters of glm.nb()
mut_lm_m_all <- glm(motif_mod_form, data=a3a1, family="poisson")
mut_lm_fe_all <- glm(feat_mod_form, data=a3a1, family="poisson")
mut_lm_f_all <- glm(full_mod_form, data=a3a1, family="poisson")

# Calculate McFadden's pseudo R-squared (unadjusted)
mrsq <- 1-mut_lm_m_all$deviance/mut_lm_m_all$null.deviance
fersq <- 1-mut_lm_fe_all$deviance/mut_lm_fe_all$null.deviance
frsq <- 1-mut_lm_f_all$deviance/mut_lm_f_all$null.deviance
maic <- AIC(mut_lm_m_all)
feaic <- AIC(mut_lm_fe_all)
faic <- AIC(mut_lm_f_all)

cat("Adjusted R-squared of combined model (motifs): ", mrsq, "\n")
cat("AIC of combined model (motifs): ", maic, "\n")

cat("Adjusted R-squared of combined model (features): ", fersq, "\n")
cat("AIC of combined model (features): ", feaic, "\n")

cat("Adjusted R-squared of combined model (full): ", frsq, "\n")
cat("AIC of combined model (full): ", faic, "\n")

##############################################################################
# Model each category independently
##############################################################################
compare.all <- data.frame()
compare.err <- data.frame()
logliks <- data.frame()
mut_cats <- unique(agg_5bp_100k$Category2)

for(i in 1:length(mut_cats)) {
	cat1 <- mut_cats[i]
	aggcat <- a3[a3$Category2==mut_cats[i],]

	if(grepl("^AT", cat1)) {
		aggcatm <- merge(aggcat, binsAT, by=c("CHR", "BIN"))
		mcols <- atcols
	} else if(grepl("^GC", cat1)) {
		aggcatm <- merge(aggcat, binsGC, by=c("CHR", "BIN"))
		mcols <- gccols
	} else {
		aggcatm <- merge(aggcat, binscpgGC, by=c("CHR", "BIN"))
		mcols <- cpggccols
	}

	aggcat$exp_s <- sample(aggcat$obs, length(aggcat$obs), replace=T)
	aggcatm$exp_s <- sample(aggcatm$obs, length(aggcatm$obs), replace=T)

	# Fit models with all data
	feat_mod_formula <- as.formula(paste("obs~",
		paste(covnames, collapse="+")))
	full_mod_formula <- as.formula(paste("obs~",
		paste(mcols[-(1:5)], collapse="+"), "+prop_GC.x+",
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
	model_dat_motif2$res <- "motifs_ext"

	model_dat_full <- data.frame(CHR=CHRA,Category2=cat1, BIN=BINA,
		exp=fits_full,
		obs=aggcatm$obs,
		exp_s=aggcatm$exp_s,
		stringsAsFactors=F)
	model_dat_full$res <- "motifs_ext+features"

	compare.all <- rbind(compare.all,
		model_dat_gc,
		model_dat_feat,
		model_dat_motif,
		model_dat_motif2,
		model_dat_full)
}

##############################################################################
# Plot model log-likelihoods
##############################################################################
ggplot(logliks, aes(x=params, y=ll, group=Category2, colour=Category2))+
	geom_line()+
	theme_bw()+
	xlab("# parameters")+
	ylab("2 log likelihood")+
ggsave("/net/bipolar/jedidiah/mutation/images/ll.png")

##############################################################################
# Process resulting data from models
##############################################################################
# Update comparison data frame
compare.all$res <- factor(compare.all$res,
	levels = c("GC", "features", "motifs", "motifs_ext", "motifs_ext+features"))
compare.all$diff <- compare.all$obs-compare.all$exp
compare.all$diff_s <- compare.all$obs-compare.all$exp_s # based on random sample
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

caout <- paste0(parentdir, "/output/", nbp, "bp_err.txt")
write.table(ca, caout, col.names=T, row.names=F, quote=F, sep="\t")

cad <- dcast(data.frame(ca), Category2~res, value.var="meanerr")

# Calculate correlation between relative rates and #motifs
rate_motif_cor <- agg_5bp_100k %>%
	mutate(rel=obs/nmotifs) %>%
	group_by(Category2) %>%
	summarise(cor=cor(rel, nmotifs, use="complete.obs"))

##############################################################################
# Plot barcharts comparing obs/exp correlation for different models
##############################################################################
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

##############################################################################
# Plot barcharts comparing 10-fold cross validation MSPE
##############################################################################
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
