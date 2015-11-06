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
# Get relative rates for 1bp and 3bp motifs and compare likelihoods
##############################################################################
cat("Getting relative rates for submotifs...\n")
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

cat("Calculating likelihoods under different motif lengths...\n")
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
cat("Plotting likelihood curves...\n")
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

# Create model formulas
# Must fix to use filtered motifs
motif_mod_form <- as.formula(paste("obs~",
	paste(names(a3a1)[-c(1:18)], collapse="+")))

feat_mod_form <- as.formula(paste("obs~",
  paste(c(covnames, "prop_GC"), collapse="+")))

full_mod_form <- as.formula(paste("obs~",
	paste(names(a3a1)[-c(1,2,16:18)], collapse="+")))

# Uses poisson regression instead of negbin, since negbin fails to converge
# with default parameters of glm.nb()
overall<-0
if(overall){
  cat("Running overall models...\n")
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
compare.all <- data.frame()
compare.err <- data.frame()
mut_cats <- unique(agg_5bp_100k$Category2)

for(i in 1:length(mut_cats)) {
	cat1 <- mut_cats[i]
  cat("Running ", cat1, "models...\n")
	aggcat <- a3[a3$Category2==mut_cats[i],]

	if(grepl("^AT", cat1)) {
		aggcatm <- merge(aggcat, binsAT, by=c("CHR", "BIN", "prop_GC"))
		mcols <- atcols
	} else if(grepl("^GC", cat1)) {
		aggcatm <- merge(aggcat, binsGC, by=c("CHR", "BIN", "prop_GC"))
		mcols <- gccols
	} else {
		aggcatm <- merge(aggcat, binscpgGC, by=c("CHR", "BIN", "prop_GC"))
		mcols <- cpggccols
	}

  # Append columns to windowed count data for all motif lengths
  bases <- c("A", "C", "G", "T")
  nts <- ifelse(grepl("AT", cat1), "A", "C")
  # Loop currently just runs for 7->5bp motif aggregation;
  # can run over 7->5->3 by setting last index to :1
  for(j in ((nbp-1)/2-1):2){

    # Specify iteration motif length
    mlength <- (j+1)*2+1

    # Define rule for substring evaluation
    griddef <- paste(c(rep("bases", j), "nts", rep("bases", j)), collapse=",")

    # Evaluate substring rule and get vector of submotifs
    tris <- apply(eval(parse(text=paste("expand.grid(",griddef,")"))),
      1, paste, collapse="")

    # Loop through each substring and append column of
    # aggregated counts
    for(k in tris){
      # Generate regex string; j is fixed per iteration
      # (e.g., looking for internal 3-mers or 5-mers)
      # so we search for all 3-mers or 5-mers by allowing
      # any base preceding or following the internal motif
      # regtri <- paste0("^", "[A-Z]{", j, "}", i, "[A-Z]{", j, "}")
      regtri <- paste0("^[A-Z]", k, "[A-Z]")

      # Extract sequences matching this submotif
      z <- names(aggcatm)[grepl(regtri, names(aggcatm))]

      # Ensure motif match vector includes only sequences
      # corresponding to the appropriate motif length
      z <- z[nchar(head(gsub("_[A-Z]*", "", z)))==mlength]

      # Create column and append to df
      tripct <- aggcatm %>%
        mutate_(.dots=setNames(paste(z, collapse="+"), k)) %>%
        select_(.dots=k)
      aggcatm <- cbind(aggcatm, tripct)
    }
  }

  # Get all 5bp motifs to use
  pset <- results %>%
    filter(Category2==cat1, Q!=1)

  # Get 5bp motifs to be expanded to constituents
  psetq1 <- results %>%
    filter(Category2==cat1, Q==1)

  # Expand significant set to 7
  hierset <- apply(expand.grid(bases, r2$Sequence, bases),
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

  # Add formulas to list
  forms <- c(gc_form, feat_form, full_form, motif_form, motif2_form)
  names(forms) <- c("gc", "feat", "full", "motif", "motif2")

  # Run models for each formula in list
  models <- runMod(forms, aggcatm)

	# 5-fold cross-validation--may need to update so expected counts are
	# re-calculated for each 1/N subset
	gc_cv <- cv.glm(data=aggcat, glmfit=mut_lm_gc, K=5)
	feat_cv <- cv.glm(data=aggcat, glmfit=mut_lm_feat, K=5)
	motif_cv <- cv.glm(data=aggcatm, glmfit=mut_lm_motif, K=5)
	full_cv <- cv.glm(data=aggcatm, glmfit=mut_lm_full, K=5)

	mspe <- c(gc_cv$delta[2], feat_cv$delta[2], motif_cv$delta[2], full_cv$delta[2])
	rmse <- sqrt(mspe)
	meanct <- mean(aggcat$obs)
	pcterr <- rmse/meanct
	mspe.res <- c("GC", "features", "motifs", "motifs+features")
	mspe.dat <- data.frame(Category2=cat1, res=mspe.res, mspe, rmse, meanct, pcterr)
	compare.err <-rbind(compare.err, mspe.dat)

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
	compare.all <- rbind(compare.all, rbind_all(moddat))
}
tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

##############################################################################
# Process resulting data from models
##############################################################################
# Update comparison data frame
# compare.all$res <- factor(compare.all$res,
# 	levels = c("GC", "features", "motifs_5bp", "motifs_5bp+features"))
compare.all$Category2 <- factor(compare.all$Category2)
# levels(compare.all$Category2) <- c(
	# "AT>CG", "AT>GC", "AT>TA",
	# "(CpG)GC>AT", "(CpG)GC>CG", "(CpG)GC>TA",
	# "GC>AT", "GC>CG", "GC>TA")

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
	# scale_colour_manual("Predictor", values=myPaletteCat(8)[4:8])+
	# scale_fill_manual("Predictor", values=myPaletteCat(8)[4:8])+
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
  scale_colour_brewer("Predictor",palette="Dark2")+
  scale_fill_brewer("Predictor", palette="Dark2")+
	# scale_colour_manual("Predictor", values=myPaletteCat(8)[5:8])+
	# scale_fill_manual("Predictor", values=myPaletteCat(8)[5:8])+
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
