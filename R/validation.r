##############################################################################
# function for running models
##############################################################################
runDNMLogit <- function(data, group){
	models <- list()
	if(group=="COMBINED"){

    # Nested models (1-5)
		models$mod_1mers_n <- glm(OBS~MU_1,
      data=data, family=binomial())
		models$mod_3mers_n <- glm(OBS~MU_1+resid3,
      data=data, family=binomial())
		models$mod_5mers_n <- glm(OBS~MU_1+resid3+resid5,
      data=data, family=binomial())
		models$mod_7mers_n <- glm(OBS~MU_1+resid3+resid5+resid7,
			data=data, family=binomial())
		models$mod_7mers_features_n <- glm(OBS~MU_1+resid3+resid5+resid7+residL,
			data=data, family=binomial())

    # Non-nested models (6-8)
		models$mod_1mers <- glm(OBS~MU_1, data=data, family=binomial())
    models$mod_3mers <- glm(OBS~MU_3, data=data, family=binomial())
		models$mod_5mers <- glm(OBS~MU_5, data=data, family=binomial())

    # Comparison of 7-mer models (9-13)
		models$mod_7mers <- glm(OBS~MU_7, data=data, family=binomial())
    models$mod_7mers_ERVs_down <- glm(OBS~MU_7D, data=data, family=binomial())
    models$mod_7mers_MAC10 <- glm(OBS~MU_7C, data=data, family=binomial())
    models$mod_7mers_AV <- glm(OBS~MU_7A, data=data, family=binomial())
    models$mod_7mers_masked <- glm(OBS~MU_7M, data=data, family=binomial())

    # Non-nested 7-mer+features model (14)
    models$mod_7mers_features <- glm(OBS~MU, data=data, family=binomial())
		models$mod_7mers_adj <- glm(OBS~MU_7adj, data=data, family=binomial())
    # models$mod_7mers_anc <- glm(OBS~MU_7AN, data=data, family=binomial())
    # models$mod_7mers_sig <- glm(OBS~MU_7+X1+X2+X3+X4+X5, data=data, family=binomial())
    # models$mod_9mers <- glm(OBS~MU_9, data=data, family=binomial())

	} else {
    # Nested models (1-4)
		models$mod_3mers_n <- glm(OBS~MU_3,
      data=data, family=binomial())
		models$mod_5mers_n <- glm(OBS~MU_3+resid5,
      data=data, family=binomial())
		models$mod_7mers_n <- glm(OBS~MU_3+resid5+resid7,
			data=data, family=binomial())
		models$mod_7mers_features_n <- glm(OBS~MU_3+resid5+resid7+residL,
			data=data, family=binomial())

    # Non-nested models (5-6)
		models$mod_3mers <- glm(OBS~MU_3, data=data, family=binomial())
		models$mod_5mers <- glm(OBS~MU_5, data=data, family=binomial())

    # Comparison of 7-mer models (7-11)
    models$mod_7mers <- glm(OBS~MU_7, data=data, family=binomial())
    models$mod_7mers_ERVs_down <- glm(OBS~MU_7D, data=data, family=binomial())
    models$mod_7mers_MAC10 <- glm(OBS~MU_7C, data=data, family=binomial())
    models$mod_7mers_AV <- glm(OBS~MU_7A, data=data, family=binomial())
    models$mod_7mers_masked <- glm(OBS~MU_7M, data=data, family=binomial())

    # Non-nested 7-mer+features model (12)
    models$mod_7mers_features <- glm(OBS~MU, data=data, family=binomial())
    # models$mod_7mers_anc <- glm(OBS~MU_7AN, data=data, family=binomial())
    # models$mod_7mers_sig <- glm(OBS~MU_7+X1+X2+X3+X4+X5, data=data, family=binomial())
    # models$mod_9mers <- glm(OBS~MU_9, data=data, family=binomial())
	}

  return(models)
}

##############################################################################
# get Nagelkerke's R^2 from model list
##############################################################################
getNR2 <- function(list){
  unlist(lapply(list, function(x)
	  NagelkerkeR2(x)))[seq(2, 2*length(list), by=2)]
}

##############################################################################
# get AIC from model list
##############################################################################
getAIC <- function(list){
  unlist(lapply(list, function(x)
  	AIC(x)))[1:length(list)]
}

##############################################################################
# function for simulating mutation data
##############################################################################
simMu <- function(data, nobs, chunksize=50000, rseed){
  success <- FALSE
  mutated <- data.frame()
  seedit <- 1
	while(!success){
    # set.seed(rseed+seedit)
		rowind <- sample(nrow(data), chunksize)
		batch <- data[rowind,]
		mu <- batch$MU

		batch$SIMOBS <- sapply(mu, function(x) rbinom(1,1,x))
    mutated <- rbind(mutated, batch[batch$SIMOBS==1,])

    success <- nrow(mutated) > nobs
    seedit <- seedit+1
	}

  # last batch will go over; sample to desired # of simulated sites
  # set.seed(rseed)
  mutated <- mutated[sample(nrow(mutated), nobs),] %>%
    mutate(OBS=0, SIM="b")
  return(mutated)
}

##############################################################################
# function samples user-defined number of non-mutated sites
##############################################################################
buildValidationData <- function(data, nsites){
  set.seed(nsites)
  outdat <- data[sample(nrow(data), nsites),] %>%
    mutate(SIM="ab", SIMOBS=0) # include in simulation analysis

  return(outdat)
}

##############################################################################
# Function merges rate estimates with validation data
##############################################################################
mergeRates <- function(chrp_c){
  # rates9 <- ratelist[[5]] %>%
  #   dplyr::select(Category=Type, SEQ=Motif, MU_9=ERV_rel_rate)
	adj2 <- 3
  rates7 <- ratelist[[4]] %>%
    mutate(SEQ7=substr(Motif, 1, 7)) %>%
    dplyr::select(Category=Type, SEQ7, MU_7=ERV_rel_rate)

  rates5 <- ratelist[[3]] %>%
    mutate(SEQ5=substr(Motif, 1, 5)) %>%
    dplyr::select(Category=Type, SEQ5, MU_5=ERV_rel_rate)

  rates3 <- ratelist[[2]] %>%
    mutate(SEQ3=substr(Motif, 1, 3)) %>%
    dplyr::select(Category=Type, SEQ3, MU_3=ERV_rel_rate)

  rates1 <- ratelist[[1]] %>%
    dplyr::select(Category=Type, MU_1=ERV_rel_rate)

  rates7A <- avrates %>%
    mutate(SEQ7=substr(Motif, 1, 7)) %>%
    dplyr::select(Category=Type, SEQ7, MU_7A=eur)

  # chrp_c <- merge(chrp_c, rates9, by=c("Category", "SEQ"), all.x=T)
  chrp_c <- merge(chrp_c, rates7A, by=c("Category", "SEQ7"), all.x=T)
  chrp_c <- merge(chrp_c, rates7, by=c("Category", "SEQ7"), all.x=T)
  chrp_c <- merge(chrp_c, rates5, by=c("Category", "SEQ5"), all.x=T)
  chrp_c <- merge(chrp_c, rates3, by=c("Category", "SEQ3"))
  chrp_c <- merge(chrp_c, rates1, by=c("Category"))

  if(exists("maskgpdat")){
    rates_mask <- maskgpdat %>%
      mutate(SEQ7=substr(Motif, 1, 7)) %>%
      dplyr::select(Category=Type, SEQ7, MU_7M=ERV_rel_rate_mask)

      chrp_c <- merge(chrp_c, rates_mask, by=c("Category", "SEQ7"), all.x=T)
  }

  rates_7C <- r5m %>%
    mutate(SEQ7=substr(Motif, 1, 7)) %>%
    dplyr::select(Category=Type, SEQ7, MU_7C=MAC10_rel_rate)

  rates_7D <- r5m %>%
    mutate(SEQ7=substr(Motif, 1, 7)) %>%
    dplyr::select(Category=Type, SEQ7, MU_7D=ERV_down_rel_rate)

  rates_anc <- ancgpdat %>%
    mutate(SEQ7=substr(Motif, 1, 7)) %>%
    dplyr::select(Category=Type, SEQ7, MU_7AN=ERV_rel_rate_anc)

  chrp_c <- merge(chrp_c, rates_7C, by=c("Category", "SEQ7"), all.x=T)
  chrp_c <- merge(chrp_c, rates_7D, by=c("Category", "SEQ7"), all.x=T)
  chrp_c <- merge(chrp_c, rates_anc, by=c("Category", "SEQ7"), all.x=T)

	chrp_c <- chrp_c %>%
		mutate(MU_7adj=-lambertW(-MU_7))

  chrp_c <- chrp_c %>%
    mutate(Category=ifelse(substr(SEQ,adj2+1,adj2+2)=="CG",
                  paste0("cpg_",Category), Category)) %>%
    mutate(BIN=ceiling(POS/binw)) %>%
    # mutate(Category =
    #     plyr::mapvalues(Category, orderedcats1, orderedcats2)) %>%
    mutate(resid3=MU_3-MU_1,
      resid5=MU_5-MU_3,
      resid7=MU_7-MU_5,
      residL=MU-MU_7,
      resid5a=MU_5-MU_1,
      resid7a=MU_7-MU_1,
      residLa=MU-MU_1)

  return(chrp_c)
}

##############################################################################
# Pipeline:
# 1. sample non-mutated input sites
# 2. merge with DNMs
# 3. append columns for additional rate estimates
# 4. Simulate DNMs with same data
##############################################################################
validationPipe <- function(nsites){
  cat("Sampling non-mutated sites...\n")
  sampled_sites <- buildValidationData(input_sites, nsites)

  cat("Merging DNMs with sub-sampled background...\n")
  eval_sites <- bind_rows(list(sampled_sites, input_dnms)) %>%
    group_by(Category) %>%
    mutate(prop=cumsum(OBS)/sum(OBS)) %>%
    arrange(MU, prop) %>%
    mutate(SEQ7 = substr(SEQ, 1, 7),
      SEQ5 = substr(SEQ, 2, 6),
      SEQ3 = substr(SEQ, 3, 5))

  cat("Appending additional rate estimates...\n")
  eval_sites <- mergeRates(eval_sites)

  cat("Generating simulated dataset...\n")
  simulated_dnms <- simMu(data=eval_sites, nobs=nrow(input_dnms), rseed=rseed)

  eval_sites_sim <- bind_rows(
      list(eval_sites[eval_sites$OBS==0,], simulated_dnms)) %>%
    group_by(Category) %>%
    mutate(prop=cumsum(OBS)/sum(OBS)) %>%
    arrange(MU, prop) %>%
    mutate(OBS=SIMOBS)

  ##############################################################################
  # run overall models
  ##############################################################################
  cat("Running combined models...\n")
  combined_models_sim <- runDNMLogit(eval_sites_sim, "COMBINED")
  combined_models <- runDNMLogit(eval_sites, "COMBINED")

  # Likelihood ratio test between nested models
  test13 <- lrtest(combined_models[[1]], combined_models[[2]])
  test35 <- lrtest(combined_models[[2]], combined_models[[3]])
  test57 <- lrtest(combined_models[[3]], combined_models[[4]])
  test7L <- lrtest(combined_models[[4]], combined_models[[5]])
  fulllist<-list(test13, test35, test57, test7L)

  pvals <- lapply(fulllist, function(x) x[[5]][2])

  combined_models_summary <- data.frame(group="COMBINED",
    category="ALL",
    nsites=nsites,
    rsq=getNR2(combined_models),
    rsq_sim=getNR2(combined_models_sim),
    aic=getAIC(combined_models), stringsAsFactors=F) %>%
    tibble::rownames_to_column("mod")

  fitfile <- paste0(parentdir, "/output/model_fit_", nsites, ".txt")
  write.table(combined_models_summary, fitfile,
    col.names=T, row.names=F, quote=F, sep="\t")

  ##############################################################################
  # run type-specific models
  ##############################################################################
  cat("Running type-specific models...\n")
  type_models_summary <- data.frame()
  lrtestdat <- data.frame()
  for(i in 1:length(orderedcats)){
    categ <- orderedcats2[i]

    type_dat_s <- eval_sites_sim %>%
      ungroup() %>%
      mutate(Category =
          plyr::mapvalues(Category, orderedcats1, orderedcats2)) %>%
      filter(Category==categ) %>%
      mutate(resid5=MU_5-MU_3, resid7=MU_7-MU_5, residL=MU-MU_7)
    type_models_sim <- runDNMLogit(type_dat_s, "TYPE")

    type_dat <- eval_sites %>%
  		mutate(Category =
  				plyr::mapvalues(Category, orderedcats1, orderedcats2)) %>%
      filter(Category==categ) %>%
      mutate(resid5=MU_5-MU_3, resid7=MU_7-MU_5, residL=MU-MU_7)
    type_models <- runDNMLogit(type_dat, "TYPE")

    test53 <- lrtest(type_models[[1]], type_models[[2]])
    test75 <- lrtest(type_models[[2]], type_models[[3]])
    testL7 <- lrtest(type_models[[3]], type_models[[4]])
    fulllist <- list(test53, test75, testL7)
    pvals <- lapply(fulllist, function(x) x[[5]][2])

    mods <- c("5-mers", "7-mers", "Logit")
    lrtests <- data.frame(category=categ, mod=mods, pvals)
    lrtestdat <- bind_rows(lrtestdat, lrtests)

    summ_tmp <- data.frame(group="TYPE",
      category=categ,
      rsq=getNR2(type_models),
      rsq_sim=getNR2(type_models_sim),
      aic=getAIC(type_models), stringsAsFactors=F) %>%
      tibble::rownames_to_column("mod")

    type_models_summary <- bind_rows(type_models_summary, summ_tmp)
  }

  typefitfile <- paste0(parentdir, "/output/model_fit_by_type", nsites, ".txt")
  write.table(type_models_summary, typefitfile,
    col.names=T, row.names=F, quote=F, sep="\t")
  return(type_models_summary)
}

##############################################################################
# Read and process data
##############################################################################
cat("Reading data...\n")
validation_file <- paste0(parentdir, "/output/predicted/validation_sites.txt")
input_sites <- read.table(validation_file, header=F, stringsAsFactors=F)
names(input_sites) <- c("CHR", "POS", "MU", "OBS", "Category", "SEQ", "ID")

# Remove sites with mu=0
input_sites <- input_sites[input_sites$MU>0,]

# Subset to non-DNMs and DNMs
cat("Splitting by DNM status...\n")

input_dnms <- input_sites[input_sites$ID!="all",] %>%
  mutate(SIM="a", SIMOBS=0)

input_sites <- input_sites[input_sites$ID=="all",]

m500k <- validationPipe(500000)
m1m <- validationPipe(1000000)
m2m <- validationPipe(2000000)
m3m <- validationPipe(3000000)
