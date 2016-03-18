extra_plots <- 0
if(extra_plots){source("expl_plots.r")}

##############################################################################
# Merge aggregate data with covariate data
##############################################################################
cat("Merging bin-wise data with covariate data...\n")
agg_cov <- merge(agg_5bp_100k, mut_cov, by=c("CHR", "BIN"))

ratiofilter <- 0
if(ratiofilter){
  agg_cov$ratio <- agg_cov$exp/agg_cov$obs
  agg_cov <- filter(agg_cov, ratio<5) %>% mutate(med=(maxn+minn)/2)
}

##############################################################################
# Get relative rates for 1bp and 3bp motifs and check for homogeneity
# among constituent motifs
##############################################################################
cat("Getting relative rates for submotifs...\n")
# Get submotifs
rates5$Seq3 <- substr(rates5$Sequence, adj, nbp-adj+1)
rates5$Seq5 <- substr(rates5$Sequence, adj-1, nbp-1)

# Recalculate rates for submotifs
rates1 <- rates5 %>%
  dplyr::select(Sequence, Category2, COUNT, num) %>%
  group_by(Category2) %>%
  summarise(COUNT1=sum(COUNT), num1=sum(num), rel_prop1=num1/COUNT1)

rates3 <- rates5 %>%
  dplyr::select(Seq3, Category2, COUNT, num) %>%
  group_by(Seq3, Category2) %>%
  summarise(COUNT3=sum(COUNT), num3=sum(num), rel_prop3=num3/COUNT3)

rates5a <- rates5 %>%
  dplyr::select(Seq5, Category2, COUNT, num) %>%
  group_by(Seq5, Category2) %>%
  summarise(COUNT5=sum(COUNT), num5=sum(num), rel_prop5=num5/COUNT5)

# Merge to single data frame
rates_full <- merge(rates5, rates1, by="Category2")
rates_full <- merge(rates_full, rates3, by=c("Seq3", "Category2"))
rates_full <- merge(rates_full, rates5a, by=c("Seq5", "Category2"))

# Identify motifs whose components have most similar relative rates
cat("Checking for homogeneity of similar motifs...\n")
count<-0
results<-data.frame()

for(i in unique(rates5$Category2)){
  aggcat <- rates5[rates5$Category2==i,]
  # names(aggcat)[1] <- "Sequence"
  nbpt <- 7

  tmpresults <- data.frame()
  for(j in unique(substr(aggcat$Sequence,2,nbpt-1))){

    ref <- unique(substr(aggcat$Sequence,adj+1,adj+1))
    dat <- aggcat[substr(aggcat$Sequence,2,nbpt-1)==j,]
    dat$exp <- dat$COUNT*(sum(dat$n)/sum(dat$COUNT))

    test <- chisq.test(dat$n, p=dat$exp/sum(dat$exp))

    datrow <- data.frame(Category2=i, Sequence=j, pval=test$p.value)
    tmpresults <- rbind(tmpresults, datrow)
  }

  tmpresults$Q <- ntile(tmpresults$pval, 10)
  results <- rbind(results, tmpresults)
}

##############################################################################
# Run likelihood analysis on full data
##############################################################################
# ra1a %>% filter(model=="L5")
logit_curves<-0
if(logit_curves){
  source("logit_curves.r")
}

##############################################################################
# Initialize data of motif counts to use as covariates in negbin model
##############################################################################
cat("Initializing data for negbin models...\n")
dat_5bp_100k$bin$CHR <- as.integer(substring(dat_5bp_100k$bin$CHR, 4))

names(dat_5bp_100k$bin) <- gsub('\\(', '_', names(dat_5bp_100k$bin))
names(dat_5bp_100k$bin) <- gsub('\\)', '_', names(dat_5bp_100k$bin))

atcols <- c(names(dat_5bp_100k$bin)[1:5],
  names(dat_5bp_100k$bin)[which(substr(names(dat_5bp_100k$bin),
    adj+1, adj+1)=="A")])

# For binsAT, etc., must subselect variables

binsAT <- dat_5bp_100k$bin %>%
  select_(.dots = atcols) %>%
	arrange(CHR, BIN)

gcdn <- c("CA", "CC", "CT")
gccols <- c(names(dat_5bp_100k$bin)[1:5],
  names(dat_5bp_100k$bin)[which(substr(names(dat_5bp_100k$bin),
    adj+1, adj+2) %in% gcdn)])

binsGC <- dat_5bp_100k$bin %>%
  select_(.dots = gccols) %>%
	arrange(CHR, BIN)

cpggccols <- c(names(dat_5bp_100k$bin)[1:5],
  names(dat_5bp_100k$bin)[which(substr(names(dat_5bp_100k$bin),
   adj+1, adj+2)=="CG")])

binscpgGC <- dat_5bp_100k$bin %>%
  select_(.dots = cpggccols) %>%
	arrange(CHR, BIN)

##############################################################################
# Run combined negbin model
##############################################################################
# Add expected #singletons per bin using 1bp motifs
# -used in comparing category-specific models with weighted mean method
a3 <- merge(agg_cov, rates1[,c(1,2,4)], by="Category2")
# a3$exp1 <- a3$nmotifs*a3$rel_prop1

# Aggregate across categories for total #singletons per bin
a3a <- a3 %>%
  # group_by_(.dots=lapply(names(a3)[c(2,3,5:17)], as.symbol)) %>%
  group_by_(.dots=lapply(danames, as.symbol)) %>%
  summarise(obs=sum(obs))

# Add motif counts per bin
a3a1 <- merge(a3a, dat_5bp_100k$bin, by=c("CHR", "BIN", "prop_GC"))

bases <- c("A", "C", "G", "T")
nts <- c("A", "C")

# Uses poisson regression instead of negbin, since negbin fails to converge
# with default parameters of glm.nb()
overall<-0
if(overall){
  cat("Running overall models...\n")

  cat("Formatting data...\n")
  a3a1 <- getSubMotifs(a3a1, nts, bases)

  # Create model formulas
  # Must fix to use filtered motifs
  motif_mod_form <- as.formula(paste("obs~",
  	paste(names(a3a1)[(ncol(a3a1)-511):ncol(a3a1)], collapse="+")))

  feat_mod_form <- as.formula(paste("obs~",
    paste(c(covnames, "prop_GC"), collapse="+")))

  full_mod_form <- as.formula(paste("obs~",
  	paste(names(a3a1)[c(3:15, (ncol(a3a1)-511):ncol(a3a1))], collapse="+")))


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
}

##############################################################################
# Model each category independently
##############################################################################
ptm <- proc.time()
cat("Running per-category models...\n")
mut_cats <- unique(agg_5bp_100k$Category2)
compare.all <- data.frame()
compare.err <- data.frame()
compare.aic <- data.frame()

for(i in 1:length(mut_cats)) {
	cat1 <- mut_cats[i]
  cat("Running ", cat1, "models...\n")
	aggcat <- a3[a3$Category2==mut_cats[i],]

	if(grepl("^AT", cat1)) {
		bindat <- binsAT
		mcols <- atcols
	} else if(grepl("^GC", cat1)) {
		bindat <- binsGC
		mcols <- gccols
	} else {
		bindat <- binscpgGC
		mcols <- cpggccols
	}

  aggcatm1 <- merge(aggcat, bindat, by=c("CHR", "BIN", "prop_GC"), all.x=T)

  # Subset 7bp rates for category i and sort
  rcat<-rates5 %>% filter(Category2==cat1) %>% arrange(Sequence)

  # Get expected num per window per bin
  z<-as.vector(rcat$rel_prop)*as.matrix(bindat[,6:ncol(bindat)])

  # Merge row sums with CHR/BIN
  # CHR BIN EXP
  r6<-cbind(bindat[,c(1,5)],marg=rowSums(z))

  bases <- c("A", "C", "G", "T")
  nts <- ifelse(grepl("^AT", cat1), "A", "C")

  b3 <- bases
  if(grepl("^cpg", cat1)){
    b3 <- c("G")
  } else if (grepl("^GC", cat1)){
    b3 <- c("A", "C", "T")
  }

  aggcatm2 <- getSubMotifs(aggcatm1, nts, b3)

  # Merge data to include column of marginals
  aggcatm<-merge(aggcatm2, r6, by=c("CHR", "BIN"))

  # Fix issue where a single bin in Chr5 with 15 AT>GC observations
  # causes glm.nb() to fail to converge
  if(cat1=="AT_GC"){
    aggcatm <- aggcatm[aggcatm$obs>15,]
  }

  # Get all 5bp motifs to use
  pset <- results %>%
    filter(Category2==cat1, Q!=1)

  # Get 5bp motifs to be expanded to constituents
  psetq1 <- results %>%
    filter(Category2==cat1, Q==1)

  # Expand significant set to 7
  hierset <- apply(expand.grid(bases, psetq1$Sequence, bases),
    1, paste, collapse="")
  hiersetrev <- sapply(hierset, revcomp)
  hierset <- paste0(hierset, "_", hiersetrev, "_")

  m5set <- c(as.character(pset$Sequence), as.character(psetq1$Sequence))
  m7set <- c(as.character(pset$Sequence), hierset)

	# Fit models with all data
  # wm_form <- as.formula("obs~exp")
  # m1_form <- as.formula("obs~exp1")
  gc_form <- as.formula("obs~prop_GC")
	feat_form <- as.formula(paste("obs~",
		paste(covnames, collapse="+")))
	full_form <- as.formula(paste("obs~",
		paste(m5set, collapse="+"), "+prop_GC+",
		paste(covnames, collapse="+")))
	motif_form <- as.formula(paste("obs~",
		paste(m5set, collapse="+")))
  motif2_form <- as.formula(paste("obs~",
		paste(m7set, collapse="+")))
  full_form_int <- as.formula(paste("obs~prop_GC*(",
		paste(m5set, collapse="+"), ")+",
		paste(covnames, collapse="+")))
  marg_form <- as.formula("obs~marg")

  # Add formulas to list
  forms <- c(gc_form, feat_form,
    full_form, full_form_int,
    motif_form, motif2_form,
    marg_form)
  names(forms) <- c("gc", "features",
    "full", "full_gc_inter",
    "motifs5", "motifs5_top7",
    "marginal7")

  # Run models for each formula in list
  models <- runMod(forms, aggcatm)

  aics <- sapply(models, function(x) AIC(x))
  aicdf <- data.frame(Category2=cat1, model=names(aics), AIC=aics)
  compare.aic <- rbind(compare.aic, aicdf)
	# 5-fold cross-validation--may need to update so expected counts are
	# re-calculated for each 1/N subset
	# gc_cv <- cv.glm(data=aggcatm, glmfit=models$gc, K=5)
	# feat_cv <- cv.glm(data=aggcatm, glmfit=models$feat, K=5)
	# motif_cv <- cv.glm(data=aggcatm, glmfit=models$motif, K=5)
	# full_cv <- cv.glm(data=aggcatm, glmfit=models$full, K=5)
  #
	# mspe <- c(gc_cv$delta[2], feat_cv$delta[2], motif_cv$delta[2], full_cv$delta[2])
	# rmse <- sqrt(mspe)
	# meanct <- mean(aggcat$obs)
	# pcterr <- rmse/meanct
	# mspe.res <- c("GC", "features", "motifs", "motifs+features")
	# mspe.dat <- data.frame(Category2=cat1, res=mspe.res, mspe, rmse, meanct, pcterr)
	# compare.err <-rbind(compare.err, mspe.dat)

  # Get fitted values from each model and name with CHR/BIN
  fits <- getFits(models, aggcatm)

	BIN <- as.integer(gsub(".*\\.", "", names(fits$feat)))
	CHR <- as.integer(gsub("\\..*", "", names(fits$feat)))

  # Build list of dataframes for each model
  moddat <- buildDF(fits, aggcatm)

  # Add column specifying model
  for(i in 1:length(names(moddat))){
    moddat[[i]]$res <- names(moddat)[i]
  }

  # Append model predictions to full df
	compare.all <- rbind(compare.all, bind_rows(moddat))
}
tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

# Specify category as factor
compare.all$Category2 <- factor(compare.all$Category2)
# levels(compare.all$Category2) <- c(
	# "AT>CG", "AT>GC", "AT>TA",
	# "(CpG)GC>AT", "(CpG)GC>CG", "(CpG)GC>TA",
	# "GC>AT", "GC>CG", "GC>TA")

overallcor<-compare.all %>%
  group_by(res) %>%
  summarise(cor=cor(exp,obs), corsq=cor^2)

ggplot(compare.all, aes(x=obs, y=exp, group=res, colour=res))+
  geom_point(alpha=0.3)+
  scale_colour_brewer("Model",palette="Dark2")+
  facet_wrap(~res, scales="free")+
  theme_bw()
ggsave("/net/bipolar/jedidiah/mutation/images/1mb_scatter.png")


##############################################################################
# Plot barcharts comparing obs/exp correlation for different models
##############################################################################
mod.corr <- compare.all %>%
	# filter(exp > 100) %>%
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
  scale_colour_brewer("Predictor",palette="Dark2")+
  scale_fill_brewer("Predictor", palette="Dark2")+
	xlab("Category")+
	ylab("Correlation with observed count")+
	# geom_errorbar(limits, position=dodge, width=0.25)+
	theme_bw()+
	theme(legend.title = element_text(size=18),
		legend.text = element_text(size=16),
		axis.title.x = element_text(size=20),
		axis.title.y = element_text(size=20),
		axis.text.y = element_text(size=16),
		axis.text.x = element_text(size=16, angle = 45,  vjust=1, hjust=1.01))

modelbar <- paste0(parentdir, "/images/gw_5bp_vs_mod_3.png")
ggsave(modelbar, width=7, height=7)

# Get AIC for each model/category
compare.aic <- compare.aic %>% spread(Category2, AIC)

##############################################################################
# Plot per-chromosome variation
##############################################################################
cat("Plotting per-chromosome variation...\n")
require(ggbio)
data(hg19IdeogramCyto, package = "biovizBase")
ic2<-as.data.frame(hg19IdeogramCyto)

# Read in chromosome lengths for specifying plot range
chrlen<-read.table("/net/bipolar/jedidiah/mutation/reference_data/hg19.genome", header=T, stringsAsFactors=F)

for(i in 1:22){
  # Specify chromosome
  chrname<-paste0("chr",i)

  #Subset and update data
  d2<-compare.all %>%
    filter(CHR==i, res=="full") %>%
    mutate(ratio=exp/obs)
  d2$BIN<-d2$BIN*1000000

  # Get chromosome band
  ic2chr<-ic2[ic2$seqnames==chrname,]
  ic2chr$mean<-(ic2chr$start+ic2chr$end)/2

  # Plot scatterplot & loess curve, with banding overlay
  p2<-ggplot()+
    geom_rect(data=ic2chr,
      aes(xmin=start, xmax=end,
        # ymin=min(d2$ratio)-0.3, ymax=min(d2$ratio)-0.1, fill=gieStain),
        ymin=0.9, ymax=1.1, fill=gieStain),
      alpha=0.6)+
    scale_fill_manual(values=getOption("biovizBase")$cytobandColor)+
    geom_point(data=d2[d2$ratio<2,],
        aes(x=BIN, y=ratio, colour=Category2, group=Category2), alpha=0.3)+
    geom_smooth(data=d2[d2$ratio<2,],
        aes(x=BIN, y=ratio, colour=Category2, group=Category2), span=0.2, se=FALSE)+
    scale_colour_manual("Category", values=myPaletteCat(9))+
    scale_x_continuous(limits=c(0,chrlen[chrlen$chrom==chrname,]$size))+
    # scale_y_continuous(limits=c(0.25,2))+
    geom_hline(yintercept=1)+
    geom_hline(yintercept=0.9, linetype="dashed")+
    geom_hline(yintercept=1.1, linetype="dashed")+
    geom_text(data=ic2chr,
      aes(x=mean, y=min(d2$ratio)-0.3, label=name),
      hjust=0, angle=90, size=2)+
    # facet_wrap(~Category2, scales="free", ncol=1)+
    ylab("Predicted:Observed Ratio")+
    xlab(NULL)+
    theme_bw()+
    theme(legend.position="top")

  p2
  ratiofile<-paste0("/net/bipolar/jedidiah/mutation/images/chr",i,"_ratio.png")
  ggsave(ratiofile, width=10, height=5)
}
