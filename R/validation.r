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
    models$mod_7mers_anc <- glm(OBS~MU_7AN, data=data, family=binomial())
    # models$mod_7mers_sig <- glm(OBS~MU_7+X1+X2+X3+X4+X5, data=data, family=binomial())
    # models$mod_9mers <- glm(OBS~MU_9, data=data, family=binomial())

		models$mod_int <- glm(OBS~1,
			data=data, family=binomial())

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
    models$mod_7mers_anc <- glm(OBS~MU_7AN, data=data, family=binomial())
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

  # rates_anc <- ancgpdat %>%
  #   mutate(SEQ7=substr(Motif, 1, 7)) %>%
  #   dplyr::select(Category=Type, SEQ7, MU_7AN=ERV_rel_rate_anc)

	rates_anc <- test_anc %>%
    mutate(SEQ7=substr(Motif, 1, 7)) %>%
    dplyr::select(Category=Type, SEQ7, MU_7AN=estimate2)

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
    aic=getAIC(combined_models),
		aic_sim=getAIC(combined_models_sim), stringsAsFactors=F) %>%
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
      aic=getAIC(type_models),
			aic_sim=getAIC(type_models_sim),
			stringsAsFactors=F) %>%
      tibble::rownames_to_column("mod")

    type_models_summary <- bind_rows(type_models_summary, summ_tmp)
  }

  typefitfile <- paste0(parentdir, "/output/model_fit_by_type", nsites, ".txt")
  write.table(type_models_summary, typefitfile,
    col.names=T, row.names=F, quote=F, sep="\t")
  return(list(combined_models_summary, type_models_summary))
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


newmodnames <- c("1-mers", "3-mers", "5-mers", "7-mers", "7-mers+features", "1000G 7-mers")

combined_barplots <- m1m[[1]] %>%
	filter(grepl("_n|AV", mod)) %>%
	mutate(daic=aic-min(aic)) %>%
	# mutate(max=max(aic), aic_col=normalize(aic)+0.1) %>% #data.frame
	# filter(!grepl("down", mod)) %>%
	# filter(!grepl("masked", mod)) %>%
	# filter(!grepl("adj", mod)) %>%
	# mutate(mod=factor(mod, levels=unique(m1m[[1]]$mod)[c(6,7,8,9,14,11,12)])) %>%
	mutate(mod=plyr::mapvalues(mod, unique(m1m[[1]]$mod)[c(1:5,12)], newmodnames)) %>%
	# mutate(mod=factor(mod, levels=unique(m1m[[1]]$mod)[c(1:5)])) %>%
	mutate(mod=factor(mod, levels=newmodnames)) %>%
	mutate(category=paste0(category, "\n")) %>%
	filter(!grepl("features",mod))
	# mutate(data=ifelse(mod %in% unique(m1m[[1]]$mod)[c(11,12)], "poly", "ERV")) %>%

# cbp1 <- ggplot(combined_barplots, aes(x=mod, y=aic_col, group=mod, fill=mod))+
cbp1 <- ggplot(combined_barplots, aes(x=mod, y=daic, group=mod, fill=mod))+
	# geom_point(aes(shape=data), size=4)+
	# geom_point(size=4)+
	geom_col(position="dodge", width=0.8)+
	annotate("rect", xmin=0, xmax=4.5, ymin=0, ymax=Inf, alpha=0.2, fill="grey20")+
	# geom_vline(xintercept=1.5, linetype="dashed")+
	geom_hline(yintercept=10, linetype="dashed")+
	# geom_vline(xintercept=1.5, linetype="dashed")+
	# scale_y_continuous(expand = c(0,0), limits=c(0,1.2))+
	scale_y_continuous(trans = "asinh",
		breaks=c(1,10,100,1000,10000),
		expand = c(0,0))+
	facet_wrap(~category, scales="free")+
	# scale_fill_manual(values=iwhPalette[c(8:9,2:5)])+
	scale_fill_manual(values=c(viridis(4), viridis(25)[24]))+
	# scale_fill_viridis(discrete=TRUE)+
	# scale_colour_viridis(discrete=TRUE)+
	ylab(expression(paste(Delta,AIC)))+
	theme_bw()+
	theme(
		# axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
		axis.text.x=element_blank(),
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.text.y=element_text(size=14),
		# axis.text.y=element_blank(),
		axis.title.x=element_blank(),
		legend.title=element_blank(),
		axis.title.y=element_text(size=16),
		strip.text.x=element_text(size=16),
		legend.text=element_text(size=14))

cbp1 + theme(legend.position="none")
ggsave("/net/bipolar/jedidiah/mutation/images/testbar.png", width=2, height=6)

cbp1a <- cbp1 + theme(legend.position="bottom")
legend <- get_legend(cbp1a)
png(paste0(parentdir, "/images/barplot_legend.png"),
	height=2, width=6, units="in", res=300)
grid.draw(legend)
dev.off()

newmodnames_facet <- c("1-mers", "3-mers", "5-mers", "7-mers", "7-mers+features", "1000G 7-mers")
#
# dummy <- m1m[[2]] %>%
# 	# dplyr::select(mod, category, aic) %>%
# 	filter(grepl("_n|AV", mod)) %>%
# 	group_by(category) %>%
# 	mutate(mod="dummy", aic=min(aic)-10) %>%
# 	data.frame
oc3 <- c("ALL", gsub(" ", "\n", orderedcats2))
# rbind(dummy, m1m[[2]]) %>%
facet_barplots <- rbind(m1m[[1]][,-4], m1m[[2]]) %>%
	dplyr::select(mod, category, aic, rsq) %>%
	filter(grepl("_n|AV|dummy", mod)) %>%
	group_by(category) %>%
	mutate(daic=aic-min(aic)) %>%
	mutate(drsq=-(rsq-max(rsq))) %>%
	mutate(mod=plyr::mapvalues(mod, unique(m1m[[1]]$mod)[c(1:5,12)], newmodnames_facet)) %>%
	mutate(mod=factor(mod, levels=newmodnames_facet)) %>%
	ungroup() %>%
	mutate(category=gsub(" ", "\n", category)) %>%
	mutate(category=factor(category, levels=oc3)) %>%
	mutate(Data=factor(ifelse(grepl("1000G", mod), "1000G EUR Intergenic SNVs", "BRIDGES ERVs"), levels=c("BRIDGES ERVs", "1000G EUR Intergenic SNVs")))

m500c <- rbind(m500k[[1]][,-4], m500k[[2]])
m500c$nsites <- 500000

m1c <- rbind(m1m[[1]][,-4], m1m[[2]])
m1c$nsites <- 1000000

m2c <- rbind(m2m[[1]][,-4], m2m[[2]])
m2c$nsites <- 2000000

m3c <- rbind(m3m[[1]][,-4], m3m[[2]])
m3c$nsites <- 3000000

facet_barplots2 <- rbind(m500c, m1c, m2c, m3c) %>%
	dplyr::select(mod, category, aic, rsq, nsites) %>%
	filter(grepl("_n|AV|dummy", mod)) %>%
	group_by(category) %>%
	mutate(daic=aic-min(aic)) %>%
	mutate(drsq=-(rsq-max(rsq))) %>%
	mutate(mod=plyr::mapvalues(mod, unique(m1m[[1]]$mod)[c(1:5,12)], newmodnames_facet)) %>%
	mutate(mod=factor(mod, levels=newmodnames_facet)) %>%
	ungroup() %>%
	mutate(category=gsub(" ", "\n", category)) %>%
	mutate(category=factor(category, levels=oc3)) %>%
	mutate(Data=factor(ifelse(grepl("1000G", mod), "1000G EUR Intergenic SNVs", "BRIDGES ERVs"), levels=c("BRIDGES ERVs", "1000G EUR Intergenic SNVs")))


library(scales)
asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x),
            inverse = function(x) sinh(x))
}

ratelist[[2]] %>%
	mutate(category=ifelse(substr(Motif,2,3)=="CG",
								paste0("cpg_",Type), Type)) %>%
	group_by(category) %>%
	summarise(nERVs=sum(nERVs), nMotifs=sum(nMotifs)) %>%
	ungroup() %>%
	mutate(prop=nERVs/sum(nERVs))

dnm_type_props <- input_dnms %>%
	mutate(category=ifelse(substr(SEQ,adj2+1,adj2+2)=="CG",
								paste0("cpg_",Category), Category)) %>%
	group_by(category) %>%
	summarise(n=n()) %>%
	ungroup() %>%
	mutate(category =
			plyr::mapvalues(category, orderedcats1, orderedcats2)) %>%
	mutate(category=gsub(" ", "\n", category)) %>%
	mutate(prop=paste0(round(n/sum(n)*100, 1), "%")) %>%
	mutate(propn=paste0(n, "\n(", round(n/sum(n)*100, 1), "%)")) %>%
	rbind(c("ALL", 46813, "100%", "46813\n(100%)")) %>%
	mutate(category=factor(category, levels=oc3)) %>%
	arrange(category)

breakvals <- unique(c(0, rep(10^(c(0,1,2,3)), each=10)*1:10))
breaklabs <- ifelse(grepl("1", breakvals), breakvals, "")
breaklabs[1] <- "0\n(optimal model)"
breaklabs[2] <- ""

pointfig <- ggplot()+
	geom_vline(xintercept=10, linetype="dashed")+
	geom_point(data=facet_barplots, aes(x=daic, y=category), colour="white")+
	annotate("rect", xmin=-0.2, xmax=0.2, ymin=-Inf, ymax=Inf, alpha=0.6, fill="grey30")+
	annotate("rect", xmin=0.2, xmax=10, ymin=-Inf, ymax=Inf, alpha=0.6, fill="grey80")+
	annotate("rect", xmin=30000, xmax=Inf, ymin=-Inf, ymax=Inf, fill="white", colour="black")+
	# geom_point(aes(colour=mod, shape=mod),
	# 	alpha=0.9, size=5)+
	# geom_segment(data=facet_barplots,
	# 	aes(x=0,xend=daic, y=mod, yend=mod, colour=mod))+
	geom_point(data=facet_barplots,
		aes(x=daic, y=mod, group=mod, colour=mod, fill=mod, shape=mod),
		alpha=0.9, size=4, stroke=1)+
	geom_text(data=dnm_type_props, aes(x=50000, y=3, label=prop))+
		# shape=21, alpha=0.9, size=2, stroke=3)+
	facet_wrap(~category, ncol=1, scales="free_y", strip.position="left")+
	# scale_y_discrete(position = "right", breaks=dnm_type_props$category, labels=dnm_type_props$prop)+
	scale_x_continuous(trans = "biexp",
		expand=c(0,0),
		breaks=breakvals,
		labels=breaklabs)+
	scale_fill_manual("Model", values=brewer.pal(10, "Spectral")[c(10:8,3,1,3)])+
	# scale_colour_manual("Model", values=viridis(option="magma", 5, begin=0.1, end=0.8)[c(1:5,4)])+
	# scale_colour_manual("Model", values=c("#af8c3a","#57af6c","#6780d8","#8750a6","#c34f41","#8750a6"))+
	# scale_colour_manual("Model", values=c("#005cbf","#6597d3","#6db840","#ee9800","#ff657c","#ee9800"))+
	# scale_colour_manual("Model", values=c("darkblue", "blue", "green", "orange", "red", "orange"))+
	scale_colour_manual("Model", values=c(brewer.pal(10, "Spectral")[c(10:8,3,1)], "black"))+
	scale_shape_manual("Model", values=c(rep(21,5),24))+
	xlab(expression(paste(Delta,AIC)))+
	guides(colour = guide_legend(nrow = 2),
		fill = guide_legend(nrow = 2))+
	coord_cartesian(xlim=c(-0.2, 80000))+
	labs(subtitle="% of total\n de novo mutations")+
	# labs(subtitle="de novo mutations predicted\n(% of total)")+
	# labs(subtitle=expression(paste("% of total\n", italic(de), " ", italic(novo), " mutations")))+
	theme_bw()+
	theme(
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		plot.subtitle = element_text(hjust=1),
		axis.ticks.y=element_blank(),
		axis.text.x=element_text(size=14),
		axis.text.y=element_blank(),
		axis.title.x=element_text(size=16),
		legend.title=element_blank(),
		axis.title.y=element_blank(),
		strip.background = element_rect(colour = "black", fill = "white"),
		strip.text.y=element_text(size=14, angle=180),
		legend.box.background = element_rect(),
		legend.text=element_text(size=14),
		legend.position="bottom")

pointfig
# cvd_grid(pointfig)
ggsave("/net/bipolar/jedidiah/mutation/images/testpoint_facet.png", width=8, height=8)

pointfig2 <- ggplot()+
	geom_vline(xintercept=1e-04, linetype="dashed", colour="white")+
	geom_point(data=facet_barplots, aes(x=drsq, y=category), colour="white")+
	geom_point(data=facet_barplots,
		aes(x=drsq, y=category, group=mod, colour=mod, fill=mod, shape=mod),
		alpha=0.9, size=4, stroke=1)+
	facet_wrap(~category, ncol=1, scales="free_y", strip.position="left")+
	scale_x_log10()+
	scale_fill_manual("Model", values=brewer.pal(10, "Spectral")[c(10:8,3,1,3)])+
	scale_colour_manual("Model", values=c(brewer.pal(10, "Spectral")[c(10:8,3,1)], "black"))+
	scale_shape_manual("Model", values=c(rep(21,5),24))+
	xlab("decrease in R^2")+
	guides(colour = guide_legend(nrow = 2),
		fill = guide_legend(nrow = 2))+
	theme_bw()+
	theme(
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		plot.subtitle = element_text(hjust=1),
		axis.ticks.y=element_blank(),
		axis.text.x=element_text(size=14),
		axis.text.y=element_blank(),
		axis.title.x=element_text(size=16),
		legend.title=element_blank(),
		axis.title.y=element_blank(),
		strip.background = element_rect(colour = "black", fill = "white"),
		strip.text.y=element_text(size=14, angle=180),
		legend.box.background = element_rect(),
		legend.text=element_text(size=14),
		legend.position="bottom")

pointfig2
ggsave("/net/bipolar/jedidiah/mutation/images/testpoint_facet_rsq.png", width=8, height=8)


# ggplot(facet_barplots2, aes(x=rsq, y=mod, group=nsites, colour=mod, fill=mod))+
ggplot()+
	geom_point(data=facet_barplots2,
		aes(x=rsq, y=factor(nsites), group=mod, colour=mod, fill=mod, shape=mod),
		alpha=0.9, size=4, stroke=1)+
	# geom_bar(stat="identity", position="dodge")+
	scale_shape_manual("Model", values=c(rep(21,5),24))+
	scale_fill_manual("Model", values=brewer.pal(10, "Spectral")[c(10:8,3,1,3)])+
	scale_colour_manual("Model", values=c(brewer.pal(10, "Spectral")[c(10:8,3,1)], "black"))+
	# facet_grid(category~nsites)+
	facet_wrap(~category, ncol=1, strip.position="left")+
	theme_bw()+
	theme(
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		plot.subtitle = element_text(hjust=1),
		axis.ticks.y=element_blank(),
		axis.text.x=element_text(size=14),
		# axis.text.y=element_blank(),
		axis.title.x=element_text(size=16),
		legend.title=element_blank(),
		axis.title.y=element_blank(),
		strip.background = element_rect(colour = "black", fill = "white"),
		strip.text.y=element_text(size=14, angle=180),
		legend.box.background = element_rect(),
		legend.text=element_text(size=14),
		legend.position="bottom")

ggsave("/net/bipolar/jedidiah/mutation/images/model_rsq_nsites_pt.png", width=8, height=12)

	# mutate(data=ifelse(mod %in% unique(m1m[[1]]$mod)[c(11,12)], "poly", "ERV")) %>%
	# ggplot(aes(x=mod, y=aic_col, group=category, fill=mod))+
	ggplot(facet_barplots, aes(x=mod, y=daic, group=category))+
		# geom_rect(aes(xmin=0,xmax=4,ymin=0,ymax=1000,fill="grey20"), colour="black")+
		geom_bar(aes(fill=mod), stat="identity", position="dodge", width=0.8)+
		# annotate("rect", xmin=0, xmax=3.5, ymin=0, ymax=Inf, alpha=0.2, fill="grey20")+
		# geom_vline(xintercept=1.5, linetype="dashed")+
		geom_hline(yintercept=0)+
		geom_hline(yintercept=10, linetype="dashed")+
		geom_hline(yintercept=-10, linetype="dashed")+
		# geom_point(size=4, aes(colour=mod))+
		scale_y_continuous(trans = "asinh",
			breaks=c(-1000,-100,-10,-1,0,1,10,100,1000),
			expand = c(0,0))+
		# scale_y_continuous(expand = c(0,0), limits=c(0,1.2))+
		# scale_y_continuous(expand = c(0,0))+
		# scale_y_continuous(breaks=seq(331000,355000,by=2000))+
		facet_wrap(~category, ncol=9)+
		# facet_wrap(~category, scales="free", dir="v")+
		# scale_colour_manual(values=iwhPalette[1:7])+
		# scale_fill_viridis(discrete=TRUE)+
		# scale_fill_manual(values=c(viridis(4)[2:4], viridis(25)[24]))+
		scale_fill_manual(values=viridis(5)[2:5])+
		# scale_colour_viridis(discrete=TRUE)+
		# scale_colour_manual(values=iwhPalette[c(8,2:5)])+
		# scale_fill_manual(values=iwhPalette[c(8,2:5)])+
		ylab(expression(paste(Delta,AIC)))+
		theme_bw()+
		theme(
			# axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
			panel.grid.minor = element_blank(),
			panel.grid.major.x = element_blank(),
			axis.text.x=element_blank(),
			axis.text.y=element_text(size=14),
			# axis.text.y=element_blank(),
			axis.title.x=element_blank(),
			legend.title=element_blank(),
			axis.title.y=element_text(size=16),
			strip.text.x=element_text(size=16),
			legend.text=element_text(size=14),
			legend.position="bottom")

ggsave("/net/bipolar/jedidiah/mutation/images/testbar_facet.png", width=14, height=6)

facet_barplots2 <- m1m[[2]] %>%
	dplyr::select(mod, category, aic) %>%
	filter(grepl("7", mod)) %>%
	filter(grepl("_n", mod)) %>%
	group_by(category) %>%
	spread(mod, aic) %>%
	mutate(aicref=mod_7mers_n.R2) %>%
	gather(mod, aic, mod_7mers_features_n.R2:mod_7mers_n.R2) %>%
	mutate(daic=aicref-aic) %>%
	# mutate(daic=ifelse(daic<0,0.01,daic)) %>%
	# mutate(max=max(aic), aic_col=normalize(aic)+0.1) %>% #data.frame
	# filter(!grepl("down", mod)) %>%
	# filter(!grepl("masked", mod)) %>%
	# filter(!grepl("adj", mod)) %>%
	# mutate(mod=factor(mod, levels=unique(m1m[[1]]$mod)[c(6,7,8,9,14,11,12)])) %>%
	# mutate(mod=plyr::mapvalues(mod, unique(m1m[[2]]$mod)[c(10,1:4)], newmodnames_facet)) %>%
	mutate(mod=plyr::mapvalues(mod, unique(m1m[[2]]$mod)[c(3:4)], newmodnames_facet[3:4])) %>%
	# mutate(mod=factor(mod, levels=unique(m1m[[1]]$mod)[c(1:5)])) %>%
	mutate(mod=factor(mod, levels=newmodnames_facet)) %>%
	ungroup() %>%
	mutate(category=gsub(" ", "\n", category)) %>%
	mutate(category=factor(category, levels=oc3)) %>%
	filter(grepl("features",mod))

ggplot(facet_barplots2, aes(x=mod, y=daic, group=category))+
	# geom_rect(aes(xmin=0,xmax=4,ymin=0,ymax=1000,fill="grey20"), colour="black")+
	geom_bar(aes(fill=mod), stat="identity", position="dodge", width=0.8)+
	# annotate("rect", xmin=0, xmax=3.5, ymin=0, ymax=Inf, alpha=0.2, fill="grey20")+
	# geom_vline(xintercept=1.5, linetype="dashed")+
	geom_hline(yintercept=0)+
	geom_hline(yintercept=10, linetype="dashed")+
	# geom_hline(yintercept=-10, linetype="dashed")+
	# geom_point(size=4, aes(colour=mod))+
	# scale_y_continuous(trans = "asinh",
	scale_y_continuous(#trans = "biexp",
		breaks=c(0,1,10,100,1000),
		expand = c(0,0))+
	# scale_y_continuous(expand = c(0,0), limits=c(0,1.2))+
	# scale_y_continuous(expand = c(0,0))+
	# scale_y_continuous(breaks=seq(331000,355000,by=2000))+
	facet_wrap(~category, ncol=9)+
	# facet_wrap(~category, scales="free", dir="v")+
	# scale_colour_manual(values=iwhPalette[1:7])+
	# scale_fill_viridis(discrete=TRUE)+
	# scale_fill_manual(values=c(viridis(4)[2:4], viridis(25)[24]))+
	scale_fill_manual(values=viridis(5)[5])+
	# scale_colour_viridis(discrete=TRUE)+
	# scale_colour_manual(values=iwhPalette[c(8,2:5)])+
	# scale_fill_manual(values=iwhPalette[c(8,2:5)])+
	ylab(expression(paste(Delta,AIC,"(ERV 7-mer model as reference)")))+
	theme_bw()+
	theme(
		# axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.text.x=element_blank(),
		axis.text.y=element_text(size=14),
		# axis.text.y=element_blank(),
		axis.title.x=element_blank(),
		legend.title=element_blank(),
		axis.title.y=element_text(size=16),
		strip.text.x=element_text(size=16),
		legend.text=element_text(size=14),
		legend.position="bottom")

ggsave("/net/bipolar/jedidiah/mutation/images/testbar_facet_7vF.png", width=14, height=6)


m1m[[2]] %>%
	filter(grepl("_n|AV|dummy", mod)) %>%
	group_by(category) %>%
	do(reshape2::melt(outer(.$aic, .$aic, "-"))) %>%
	filter(Var2>=Var1) %>%
	ggplot(aes(x=Var2, y=Var1, fill=value))+
		geom_tile()+
		facet_wrap(~category)+
		scale_fill_gradientn(colours=brewer.pal(10,"PiYG"),
			trans = "asinh",
			breaks=c(-1000,-100,-10,-1,0,1,10,100,1000))

 ggsave("/net/bipolar/jedidiah/mutation/images/aicheat.png")
