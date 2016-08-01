require(ggplot2)
require(dplyr)
require(tidyr)

writecats <- 0

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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
  		aes(x=(5000-ntile)/5000, y=1-prop, group=group, colour=group), size=1.2)+
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
    group_by(Category.x) %>%
    arrange_(mucol) %>%
    mutate(prop=cumsum(OBS)/sum(OBS)) %>%
    arrange_(mucol, "prop") %>%
    mutate(mutmp=mucol, ntile=ntile(mutmp, 5000), group=groupname) %>%
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
  outdat <- bind_rows(list(chrp_gfasc, chrp_erv, chrp_common, chrp_av))

  if(sim){
    chrp_sim_max <- simMu(data, nsim=10, nobs=100000, chunksize=50000)
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
  		mu <- row$MU

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

##############################################################################
# Read and process data
##############################################################################

# Read data
cat("Reading data...\n")
maxc <- read.table("/net/bipolar/jedidiah/mutation/maxc_7bp.txt", header=T, stringsAsFactors=F)
maxc$Category <- gsub("cpg_", "", maxc$Category2)
chrpf <- read.table("/net/bipolar/jedidiah/mutation/output/rocdat.7bp.2.txt", header=F)
names(chrpf) <- c("CHR", "POS", "BIN", "MU", "OBS", "Category", "SEQ", "MU_C", "MU_S", "MU_A")

# Remove sites with mu=0
chrpf <- chrpf[chrpf$MU>0,]

# Read DNMs
cat("Reading DNMs...\n")
dnms_full <-  read.table("/net/bipolar/jedidiah/mutation/reference_data/DNMs/GoNL_DNMs.txt", header=T, stringsAsFactors=F)
dnms_full <- dnms_full[,1:5]
names(dnms_full) <- c("ID", "CHR", "POS", "REF", "ALT")

gdnms <- read.table("/net/bipolar/jedidiah/mutation/reference_data/DNMs/goldmann_2016_dnms.txt", header=T, stringsAsFactors=F)
gdnms$ID <- "goldmann"
gdnms <- gdnms[,c(7,1,2,4,5)]
names(gdnms) <- c("ID", "CHR", "POS", "REF", "ALT")
gdnms$CHR <- gsub("chr", "", gdnms$CHR)

dnms_full <- rbind(dnms_full, gdnms)
dnms_full$CAT <- paste(dnms_full$REF, dnms_full$ALT, sep="")

# Manually remove bins near chr20 centromere
# chr22 <- chr22[ which(chr22$BIN<260 | chr22$BIN>300),]
dnms_full$Category[dnms_full$CAT=="AC" | dnms_full$CAT=="TG"] <- "AT_CG"
dnms_full$Category[dnms_full$CAT=="AG" | dnms_full$CAT=="TC"] <- "AT_GC"
dnms_full$Category[dnms_full$CAT=="AT" | dnms_full$CAT=="TA"] <- "AT_TA"
dnms_full$Category[dnms_full$CAT=="GA" | dnms_full$CAT=="CT"] <- "GC_AT"
dnms_full$Category[dnms_full$CAT=="GC" | dnms_full$CAT=="CG"] <- "GC_CG"
dnms_full$Category[dnms_full$CAT=="GT" | dnms_full$CAT=="CA"] <- "GC_TA"

dnms_full <- dnms_full %>%
  dplyr::select(ID, CHR, POS, Category)

# Write DNM data per category with columns ID, CHR, POS, Category
# Necessary if rocdat.7bp.txt not yet created
if(writecats){
  orderedcats <- sort(unique(dnms_full$Category))
  for(i in 1:6){
    catind <- orderedcats[i]
    dnmsub <- dnms_full %>% filter(Category==catind)
    outfile <- paste0(parentdir, "/reference_data/DNMs/GoNL_", catind, ".txt")
    write.table(dnmsub, outfile, col.names=F, row.names=F, quote=F, sep="\t")
  }
}

# Duplicate data, merge with DNMs to get ID
cat("Annotating with ID...\n")
chrpf <- merge(chrpf, dnms_full, by=c("CHR", "POS"), all.x=T)
chrpf$ID[is.na(chrpf$ID)] <- "all"

# Subset to non-DNMs and DNMs
cat("Splitting by DNM status...\n")
chrpfa <- chrpf[chrpf$ID=="all",]
chrpfa <- chrpfa[sample(nrow(chrpfa), 1000000),]
chrpfdnm <- chrpf[chrpf$ID!="all",]

dnm_agg <- chrpfdnm %>%
	group_by(Category.x=as.character(Category.x), SEQ=substr(SEQ, 3, 5)) %>%
	summarise(ndnm=n())

 rates3 <- read.table("/net/bipolar/jedidiah/mutation/output/3bp_1000k_rates.txt", header=T, stringsAsFactors=F)

rates3$Category.x <- gsub("cpg_", "", rates3$Category2)
rates3$SEQ <- substr(rates3$Sequence, 1, 3)
r3m <- merge(dnm_agg, rates3, by=c("Category.x", "SEQ"))

rates_3 <- r3m %>%
	mutate(rel_prop3=ndnm/COUNT) %>%
  dplyr::select(Sequence, Category.x, rel_prop3) %>%
  spread(Category.x, rel_prop3)

rates_3[is.na(rates_3)] <- 0

write.table(rates_3, "/net/bipolar/jedidiah/mutation/dnm_3bp_rates.txt", col.names=T, row.names=F, quote=F, sep="\t")


# Recombine data
cat("Creating combined data...\n")
chrp <- rbind(chrpfdnm, chrpfa) %>%
  group_by(Category.x) %>%
  mutate(prop=cumsum(OBS)/sum(OBS)) %>%
  arrange(MU, prop)

# Process data for each model
full_auc_dat <- subWrapper(chrp, sim=T)

# Get table of AUC for each mutation type/model combo
full_auc <- full_auc_dat %>%
  group_by(group, Category.x) %>%
  summarise(AUC=1-sum(prop)/n()) %>%
  spread(Category.x, AUC)

roc_plotdat <- full_auc_dat[sample(nrow(full_auc_dat), .1*nrow(full_auc_dat)),]
# Plot common vs A&V rates

# Select only AV/BRIDGES common categories
fig_a_dat <- roc_plotdat %>%
  filter(grepl("^A|B", group))

plotROC(roc_plotdat, "/net/bipolar/jedidiah/mutation/images/pseudo_roc_curves.png")

# Plot common vs ERVs
plotROC(full_auc_dat, "/net/bipolar/jedidiah/mutation/images/pseudo_roc_curves.png")

# Plot common vs ERVs
plotROC(full_auc_dat, "/net/bipolar/jedidiah/mutation/images/pseudo_roc_curves.png")

# Plot everything together
plotROC(full_auc_dat, "/net/bipolar/jedidiah/mutation/images/pseudo_roc_curves.png")

# Subset to max contrast
chrp_maxc <- chrp %>%
  filter((SEQ %in% maxc$Sequence) & (Category.x %in% maxc$Category))

maxc_auc_dat <- subWrapper(chrp_maxc, sim=F)

maxc_auc <- maxc_auc_dat %>%
  group_by(group, Category.x) %>%
  summarise(AUC=1-sum(prop)/n()) %>%
  spread(Category.x, AUC)

# Plot all models under max contrast
plotROC(maxc_auc_dat, "/net/bipolar/jedidiah/mutation/images/pseudo_roc_curves_maxc.png")

# Plot common vs ERVs
plotROC(maxc_auc_dat, "/net/bipolar/jedidiah/mutation/images/pseudo_roc_curves.png")


coefs <- read.table("/net/bipolar/jedidiah/mutation/output/logmod_data/coefs/coefs_full.txt", header=F, stringsAsFactors=F)

names(coefs) <- c("Cov", "Est", "SE", "Z", "pval", "Sequence", "Category")

# coefs$qval <- p.adjust(coefs$pval)
splitCoefs <- function(x){
  coefs %>%
    # filter(Cov==x, pval<1e-08) %>%
		filter(Cov==x) %>%
		# mutate(qval=p.adjust(pval)) %>%
		# filter(qval<0.1) %>%
		mutate(ntile=ntile(pval, 100)) %>%
		filter(ntile<=10) %>%
		# filter(ntile<=25 | ntile>=75) %>%
		# filter(exp(Est)>1.5 | exp(Est)<0.67) %>%
    dplyr::select(Sequence, Category)
}

coef_split <- lapply(unique(coefs$Cov), function(x) splitCoefs(x))
names(coef_split) <- unique(coefs$Cov)

matchMotifs <- function(x) {
  chrp %>%
    dplyr::filter((substr(SEQ,0,7) %in% x$Sequence) & (Category.x %in% x$Category))
}

chrp_maxc_features <- lapply(coef_split, function(x) matchMotifs(x))

maxc_auc_feat <- lapply(chrp_maxc_features, function(x) subWrapper(x, sim=F))
# test <- bind_rows(maxc_auc_feat, .id="id")

aucCalc <- function(x){
  x %>%
    dplyr::filter(grepl("^E|^G", group)) %>%
    group_by(group, Category.x) %>%
    summarise(AUC=1-sum(prop)/n())
    # spread(Category.x, AUC)
}

cov_auc <- lapply(maxc_auc_feat, function(x) aucCalc(x))
cov_auc_dat <- bind_rows(cov_auc, .id="id")
cov_auc_dat <- cov_auc_dat[complete.cases(cov_auc_dat),]

cov_auc_dat <- cov_auc_dat %>% filter(!grepl("Intercept|DP", id))

ggplot(cov_auc_dat, aes(x=id, y=AUC, colour=group, fill=group, group=group))+
  geom_bar(stat="identity", position="dodge")+
	# geom_point()+
  facet_wrap(~Category.x)+
  coord_cartesian(ylim=c(0.5, 0.85))+
	coord_flip()+
	# theme_bw()
ggsave("/net/bipolar/jedidiah/mutation/images/test_bar.png")

cov_erv <- cov_auc_dat %>% filter(grepl("^E", group))
cov_logit <- cov_auc_dat %>% filter(!grepl("^E", group))
cov_logit$diff <- cov_logit$AUC-cov_erv$AUC
t.test(cov_logit$diff)

cov_logit$dir <- ifelse(cov_logit$diff<=0, "lower", "higher")
ggplot(cov_logit, aes(x=id, y=diff, colour=dir, fill=dir))+
	geom_bar(stat="identity", position="dodge")+
	scale_fill_brewer(palette="Set1")+
	scale_colour_brewer(palette="Set1")+
	# geom_point()+
  facet_wrap(~Category.x)+
  coord_cartesian(ylim=c(0.5, 0.85))+
	coord_flip()+
	geom_hline(yintercept=0, linetype="dashed")+
	xlab("AUC_[]")
	theme_bw()
ggsave("/net/bipolar/jedidiah/mutation/images/diff_bar.png", height=6, width=12)


# mapply(plotROC, maxc_auc_feat, paste0("/net/bipolar/jedidiah/mutation/images/", seq_along(max_auc_feat), "_roc.png"))
