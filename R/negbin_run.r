mut_cats <- unique(agg_5bp_100k$Category2)
compare.all <- data.frame()
compare.err <- data.frame()
compare.aic <- data.frame()

modsumm<-list()

for(j in 1:length(mut_cats)) {
	cat1 <- mut_cats[j]
  cat("Running ", cat1, "models...\n")
	aggcat <- a3[a3$Category2==cat1,]

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
  logit_form <- as.formula("obs~LSUM")

  # Add formulas to list
  forms <- c(gc_form, feat_form,
    full_form, full_form_int,
    motif_form, motif2_form,
    marg_form, logit_form)
  names(forms) <- c("gc", "features",
    "full", "full_gc_inter",
    "motifs5", "motifs5_top7",
    "marginal7", "logit")

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
	marg_nm <- aggcatm %>%
		dplyr::select(CHR, Category2, BIN, exp=marg, obs) %>%
		mutate(res="marg_nm")
	logit_nm <- aggcatm %>%
		dplyr::select(CHR, Category2, BIN, exp=LSUM, obs) %>%
		mutate(res="logit_nm")
	compare.all <- rbind(compare.all, marg_nm, logit_nm)

  modsumm[[j]] <- summary(models$full)
}
