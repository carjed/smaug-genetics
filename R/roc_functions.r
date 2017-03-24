##############################################################################
# Function checks if elements in a exist in b
# Output is binary vector of length same as b
##############################################################################
toBin <- function(a,b){
	as.numeric(is.element(b,a))
}

##############################################################################
# Function for plotting ROC curves
##############################################################################
plotROC <- function(data, outfile){
  colindex <- length(unique(data$group))+1

  ggplot()+
    geom_line(data=data,
  		# aes(x=(5000-ntile)/5000, y=1-prop, group=group, colour=group), size=1.2)+
			aes(x=(ntile)/n, y=prop, group=group, colour=group), size=1.2)+
    geom_abline(intercept=0, slope=1)+
    facet_wrap(~Category.x)+
    # geom_ribbon(data=full_auc_bounds,
  	# 	aes(x=(1000-ntile)/1000, y=1-val, ymin=1-ub, ymax=1-lb), alpha=0.2)+
  	scale_colour_manual(values=cbbPalette[2:colindex])+
    coord_cartesian(xlim=c(0,1))+
  	xlab("False Positive Rate")+
  	ylab("True Positive Rate")+
    theme_bw()+
  	theme(
  		# legend.position="none",
  		axis.title.x=element_text(size=20),
  		axis.title.y=element_text(size=20))
    # scale_y_continuous(limits=c(0,1))+
    # scale_x_continuous(limits=c(0,1000), labels=c())
  ggsave(outfile, height=9, width=12)
}

##############################################################################
# Function subsets data by specified model and re-sorts by corresponding mu
##############################################################################
subByModel <- function(data, mucol, groupname){

  outdat <- data %>%
		mutate_(mutmp=mucol) %>%
		group_by(Category.x) %>%
    arrange(desc(mutmp)) %>%
    mutate(prop=cumsum(OBS)/sum(OBS)) %>%
		# dplyr::select(-mutmp) %>%
    # arrange_(mucol, "prop") %>%
    # mutate(mutmp=mucol, ntile=dense_rank(mutmp), group=groupname) %>%
		mutate(ntile=min_rank(desc(mutmp))) %>%
		# mutate(ntile=rescale(1-mutmp, c(0, 1)), n=n()) %>%
		mutate(group=groupname) %>%
    dplyr::select(-mutmp)

    return(outdat)
}

##############################################################################
# Wrapper function for subByModel
##############################################################################
subWrapper <- function(data, sim=T){
  chrp_gfasc <- subByModel(data, "MU", "GFASC")
  chrp_erv <- subByModel(data, "MU_S", "ERV 7-mers")
  chrp_common <- subByModel(data, "MU_C", "BRIDGES polymorphism 7-mers")
  chrp_av <- subByModel(data, "MU_A", "Aggarwala 7-mers")
	chrp_erv5 <- subByModel(data, "MU_5", "ERV 5-mers")
	chrp_erv3 <- subByModel(data, "MU_3", "ERV 3-mers")
  outdat <- bind_rows(list(chrp_gfasc, chrp_erv, chrp_common, chrp_av, chrp_erv5, chrp_erv3))

  if(sim){
    chrp_sim_max <- simMu(data, nsim=1, nobs=50000, chunksize=50000)
    outdat <- rbind(data.frame(outdat), data.frame(chrp_sim_max))
  }

  return(outdat)
}

##############################################################################
# Function simulates mutations by randomly selecting sites and checking if
# Bernoulli(mu)==1; loop continues until we have simulated the same number
# of sites as in observed data
##############################################################################
simMu <- function(data, nsim, nobs, chunksize){
  nsim <- nsim
  # nobs <- sum(data$OBS)
	nsample <- chunksize

  full_auc_sim <- data.frame()
  full_auc_sim_dat <- data.frame()

  chrptmp <- data

  for(i in 1:nsim){
  	cat("Running simulation", i, "of", nsim, "...\n")
  	nsample <- chunksize
  	mutated <- c()
  	# nsites<-round(rnorm(1, params$mean, params$sd), 0)
    nsites <- nobs

  	while(length(mutated) < nsites){
  		rowind <- sample(nrow(chrptmp), nsample)
  		row <- chrptmp[rowind,]
  		mu <- row$MU_S

  		batch <- sapply(mu, function(x) rbinom(1,1,x))
  		mutated <- c(mutated, row$POS[which(as.logical(batch))])
  	}

  	chrptmp$OBS <- toBin(mutated, chrptmp$POS)

    grouptmp <- paste0("sim", i)
    chrptmpsub <- subByModel(chrptmp, "MU", grouptmp)

    full_auc_tmp <- chrptmpsub %>%
      group_by(group, Category.x) %>%
      summarise(AUC=1-sum(prop)/n())

    full_auc_sim <- rbind(full_auc_sim, as.data.frame(full_auc_tmp))
    full_auc_sim_dat <- rbind(full_auc_sim_dat, as.data.frame(chrptmpsub))
  }

  # Subset data to simulations with maximum AUC
  full_auc_sim_max <- full_auc_sim %>%
    group_by(Category.x) %>%
    filter(AUC==max(AUC))

  full_auc_sim_dat2 <- merge(full_auc_sim_dat, full_auc_sim_max, by=c("group", "Category.x"))
  full_auc_sim_dat2$group <- "max"
  full_auc_sim_dat2 <- full_auc_sim_dat2 %>%
    dplyr::select_(.dots=c(names(data), "ntile", "group"))

  return(full_auc_sim_dat2)
}
