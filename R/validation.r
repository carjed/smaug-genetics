##############################################################################
# Read and process data
##############################################################################

##############################################################################
# Add columns for 5-mer, 3-mer, and 1-mer rates
##############################################################################

# r7s <- r5m %>%
# 	dplyr::select(Category.x=Category,
# 		SEQ=Sequence,
# 		MU_7S=rel_prop,
# 		MU_7P=common_rel_prop)
# chrp <- merge(chrp, r7s, by=c("Category.x", "SEQ"))

cat("Reading data...\n")

validation_file <- paste0(parentdir, "/output/predicted/validation_sites.txt")
chrp <- read.table(validation_file, header=F, stringsAsFactors=F)
names(chrp) <- c("CHR", "POS", "BIN", "MU", "OBS", "Category", "SEQ", "ID")

# Remove sites with mu=0
chrp <- chrp[chrp$MU>0,]

# Subset to non-DNMs and DNMs
cat("Splitting by DNM status...\n")

# Filter DNMs; set flag to exclude these sites from simulation analysis
chrpdnm <- chrp[chrp$ID!="all",] %>%
  mutate(SIM="a", SIMOBS=0)
chrp <- chrp[chrp$ID=="all",]

set.seed(rseed)
chrp <- chrp[sample(nrow(chrp), 2000000),] %>%
  mutate(SIM="ab", SIMOBS=0) # include in simulation analysis

# Recombine data
cat("Creating combined data...\n")
chrp_c <- bind_rows(list(chrp, chrpdnm)) %>%
  group_by(Category) %>%
  mutate(prop=cumsum(OBS)/sum(OBS)) %>%
  arrange(MU, prop) %>%
  mutate(SEQ7 = substr(SEQ, 2, 8),
    SEQ5 = substr(SEQ, 3, 7),
    SEQ3 = substr(SEQ, 4, 6))

rates9 <- ratelist[[5]] %>%
  dplyr::select(Category=Type, SEQ=Motif, MU_9=ERV_rel_rate)

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

chrp_c <- merge(chrp_c, rates9, by=c("Category", "SEQ"), all.x=T)
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

rates_anc <- ancgpdat %>%
  mutate(SEQ7=substr(Motif, 1, 7)) %>%
  dplyr::select(Category=Type, SEQ7, MU_7AN=ERV_rel_rate_anc)

chrp_c <- merge(chrp_c, rates_anc, by=c("Category", "SEQ7"), all.x=T)

chrp_c <- chrp_c %>%
  mutate(Category=ifelse(substr(SEQ,adj+1,adj+2)=="CG",
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

simMu <- function(data, nobs, chunksize=50000, rseed){
  success <- FALSE
  mutated <- data.frame()
  seedit <- 1
	while(!success){
    set.seed(rseed+seedit)
		rowind <- sample(nrow(data), chunksize)
		batch <- data[rowind,]
		mu <- batch$MU

		batch$SIMOBS <- sapply(mu, function(x) rbinom(1,1,x))
    mutated <- rbind(mutated, batch[batch$SIMOBS==1,])

    success <- nrow(mutated) > nobs
    seedit <- seedit+1
	}

  # last batch will go over; sample to desired # of simulated sites
  set.seed(rseed)
  mutated <- mutated[sample(nrow(mutated), nobs),] %>%
    mutate(OBS=0, SIM="b")
  return(mutated)
}

chrpdnmsim <- simMu(data=chrp_c, nobs=nrow(chrpdnm), rseed=rseed)

chrp_s <- bind_rows(list(chrp_c[chrp_c$OBS==0,], chrpdnmsim)) %>%
  group_by(Category) %>%
  mutate(prop=cumsum(OBS)/sum(OBS)) %>%
  arrange(MU, prop) %>%
  mutate(OBS=SIMOBS)

gc()

##############################################################################
# Add columns for downsampled ERV 7-mers and Common 7-mers
# Must have run sing_vs_com.r
##############################################################################

runDNMLogit<-function(data, group){
	outdat <- list()
	if(nrow(data)>1e6){
		# outdat$logmod1<-glm(OBS~MU_1, data=data, family=binomial())
		# outdat$logmod3<-glm(OBS~MU_1+resid3, data=data, family=binomial())
		# outdat$logmod5<-glm(OBS~MU_1+resid3+resid5, data=data, family=binomial())
		# outdat$logmod7<-glm(OBS~MU_1+resid3+resid5+resid7,
			# data=data, family=binomial())
		# outdat$logmodL<-glm(OBS~MU_1+resid3+resid5+resid7+residL,
			# data=data, family=binomial())
		# outdat$logmodSa<-glm(OBS~MU_7S, data=data, family=binomial())
		# outdat$logmodPa<-glm(OBS~MU_7P, data=data, family=binomial())
		outdat$mod_1mers <- glm(OBS~MU_1, data=data, family=binomial())
		outdat$mod_3mers <- glm(OBS~MU_3, data=data, family=binomial())
		outdat$mod_5mers <- glm(OBS~MU_5, data=data, family=binomial())
		outdat$mod_7mers <- glm(OBS~MU_7, data=data, family=binomial())
    # outdat$mod_7mers_sig <- glm(OBS~MU_7+X1+X2+X3+X4+X5, data=data, family=binomial())
		outdat$mod_7mers_features <- glm(OBS~MU, data=data, family=binomial())
    outdat$mod_7mers_masked <- glm(OBS~MU_7M, data=data, family=binomial())
    outdat$mod_7mers_anc <- glm(OBS~MU_7AN, data=data, family=binomial())
		outdat$mod_7mers_AV <- glm(OBS~MU_7A, data=data, family=binomial())
    outdat$mod_9mers <- glm(OBS~MU_9, data=data, family=binomial())
	} else {
		outdat$logmod3<-glm(OBS~MU_3, data=data, family=binomial())
		outdat$logmod5<-glm(OBS~MU_3+resid5, data=data, family=binomial())
		outdat$logmod7<-glm(OBS~MU_3+resid5+resid7,
			data=data, family=binomial())
		outdat$logmodL<-glm(OBS~MU_3+resid5+resid7+residL,
			data=data, family=binomial())

		# outdat$logmodS<-glm(OBS~MU_7S, data=data, family=binomial())
		# outdat$logmodP<-glm(OBS~MU_7P, data=data, family=binomial())
		# outdat$logmod1a<-glm(OBS~MU_1, data=data, family=binomial())
		outdat$logmod3a<-glm(OBS~MU_3, data=data, family=binomial())
		outdat$logmod5a<-glm(OBS~MU_5, data=data, family=binomial())
		outdat$logmod7a<-glm(OBS~MU_S, data=data, family=binomial())
		outdat$logmodLa<-glm(OBS~MU, data=data, family=binomial())
		outdat$logmodAa<-glm(OBS~MU_A, data=data, family=binomial())
	}

  return(outdat)
}

overall_models <- runDNMLogit(chrp_c, "FULL")
# overall_models_sim <- runDNMLogit(chrp_s, "FULL")

rsq <- unname(unlist(lapply(overall_models, function(x)
	NagelkerkeR2(x)))[seq(2, 2*length(overall_models), by=2)])

rsqsim <- unname(unlist(lapply(overall_models_sim, function(x)
	NagelkerkeR2(x)))[seq(2, 2*length(overall_models_sim), by=2)])

test13 <- lrtest(overall_models[[1]], overall_models[[2]])
test35 <- lrtest(overall_models[[2]], overall_models[[3]])
test57 <- lrtest(overall_models[[3]], overall_models[[4]])
test7L <- lrtest(overall_models[[4]], overall_models[[5]])
fulllist<-list(test13, test35, test57, test7L)

pvals <- lapply(fulllist, function(x) x[[5]][2])


aic <- unlist(lapply(overall_models, function(x)
	AIC(x)))[1:length(overall_models)]
# mod <- c("1-mers", "3-mers", "5-mers", "7-mers", "7-mers+features",
#  "ERVs", "Common", "AV", "1-mers*", "3-mers*", "5-mers*", "7-mers*",
#  "7-mers+features*", "ERVs*", "Common*", "AV*")
combineddat<-data.frame(group="FULL",
	category=categ,
	mod=as.character(names(overall_models)),
	rsq,
	aic, stringsAsFactors=F)

nrsq <- expression(
	paste("Fraction of variance explained (Nagelkerke ", R^2, ")", sep=""))

# Define common theme for Nagelkerke R^2 plots
theme_nrsq <- function(base_size = 12, base_family = ""){
  theme_bw(base_size = base_size) %+replace%
  theme(legend.text=element_text(size=14),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=16, angle=90),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=14))
}

rsq_full_dat <- combineddat %>%
  # filter(mod %in% mod[1:5])
	filter(mod %in% mod[1:6])

rsq_full_dat <- combineddat[1:6,]

Lv7v5v3_full <- ggplot(rsq_full_dat, aes(x=category, y=rsq, fill=mod))+
  geom_bar(stat="identity", position="dodge")+
	# scale_fill_manual("Model", values=c(iwhPalette[c(1,3,4,5,9)]))+
	scale_fill_manual("Model", values=c(iwhPalette[c(1,3,4,5,9,10)]))+
	theme_nrsq()+
	ylab(nrsq)
# ggsave(paste0(parentdir, "/images/Lv7v5v3_rsq_full.png"), width=4, height=6)

# fulldat <- list()
fulldat <- data.frame()
lrtestdat <- data.frame()
for(i in 1:length(orderedcats)){
  categ <- orderedcats2[i]

  overall_dat <- chrp %>%
		mutate(Category=ifelse(substr(SEQ,adj+1,adj+2)=="CG",
										paste0("cpg_",Category), Category)) %>%
		mutate(Category =
				plyr::mapvalues(Category, orderedcats1, orderedcats2)) %>%
    filter(Category==categ) %>%
    mutate(resid5=MU_5-MU_3, resid7=MU_S-MU_5, residL=MU-MU_S)
  overall_models <- runDNMLogit(overall_dat, "FULL")

  test53 <- lrtest(overall_models[[1]], overall_models[[2]])
  test75 <- lrtest(overall_models[[2]], overall_models[[3]])
  testL7 <- lrtest(overall_models[[3]], overall_models[[4]])
  pval53 <- test53[[5]][2]
  pval75 <- test75[[5]][2]
  pvalL7 <- testL7[[5]][2]
  pvals <- c(pval53, pval75, pvalL7)
  mods <- c("5-mers", "7-mers", "Logit")
  lrtests <- data.frame(category=categ, mod=mods, pvals)
  lrtestdat <- bind_rows(lrtestdat, lrtests)

  rsq <- unlist(lapply(overall_models, function(x)
		NagelkerkeR2(x)))[seq(2,14,2)]
  aic <- unlist(lapply(overall_models, function(x) AIC(x)))
  mod <- c("3-mers", "5-mers", "7-mers", "7-mers+features",
		"ERVs", "Common", "AV")
  df <- data.frame(group="FULL", category=categ, mod, rsq, aic)

  fulldat <- bind_rows(fulldat, df)
}

fd2 <- fulldat %>%
	filter(mod %in% c("7-mers", "5-mers", "3-mers", "7-mers+features"))

fd3 <- fd2 %>%
	dplyr::select(-aic) %>%
	spread(mod, rsq)

fd3$fold<-fd3[,5]/fd3[,4]
fd3$fold2<-fd3[,5]/fd3[,3]

buildSegmentData <- function(){
	ngroups <- 9
	nmodels <- 4
	nbars <- ngroups*nmodels

	st_once <- c(seq(1, nbars, nmodels), seq(nmodels, nbars, nmodels))
	st_twice <- setdiff(1:nbars, st_once)
	st_selection <- sort(c(st_once, rep(st_twice, each=2)))

	end_none <- seq(1, nbars, nmodels)
	end_twice <- setdiff(1:nbars, end_none)
	end_selection <- sort(c(rep(end_twice, each=2)))

	vertst <- sort(c(rep(c(1:ngroups-1/9), each=2),
		rep(c(1:ngroups+1/9), each=2),
		c(1:ngroups+1/3),
		c(1:ngroups-1/3)))
	vertend <- vertst

	horizst <- sort(c(
		c(1:ngroups-1/3),
		c(1:ngroups-1/9),
		c(1:ngroups+1/9)))
	horizend <- sort(c(
		c(1:ngroups-1/9),
		c(1:ngroups+1/9),
		c(1:ngroups+1/3)))

	# Build data for geom_segment call in plot
	corplot<-data.frame(
		xst=vertst,
	  xend=vertend,
	  yst=fd2$rsq[st_selection]+.001,
	  yend=fd2$rsq[end_selection]+.01)

	corplot2 <- data.frame(
		xst=horizst,
	  xend=horizend,
	  yst=fd2$rsq[unique(end_selection)]+.01,
	  yend=fd2$rsq[unique(end_selection)]+.01,
	  pval=signif(lrtestdat$pvals, 2))

	corplot2$pval <- ifelse(corplot2$pval<0.001, "***",
		ifelse(corplot2$pval<0.01, "**", ifelse(corplot2$pval<0.05, "*", " ")))

	out <- list()
	out$corplot <- corplot
	out$corplot2 <- corplot2
}

segment_data <- buildSegmentData()
corplot <- segment_data$corplot
corplot2 <- segment_data$corplot2

# Plot pseudo-r^2 for 7-mers vs 5-mers vs 3-mers
rsqdat<-fulldat %>%
  filter(mod %in% c("7-mers", "5-mers", "3-mers", "7-mers+features")) %>%
  filter(group=="FULL") %>%
  mutate(Category =
      factor(plyr::mapvalues(category, orderedcats2, orderedcats2),
				levels=orderedcats2))
Lv7v5v3<-ggplot(rsqdat)+
  geom_bar(aes(x=Category, y=rsq, fill=mod), stat="identity", position="dodge")+
  # scale_fill_manual("Model", values=cbbPalette[c(4,6,7,8)])+
	scale_fill_manual("Model", values=c(iwhPalette[c(3,4,5,9)]))+
  geom_segment(data=corplot, aes(x=xst, xend=xend, y=yst, yend=yend))+
  geom_segment(data=corplot2, aes(x=xst, xend=xend, y=yst, yend=yend))+
  geom_text(data=corplot2,
		aes(x=xst+1/9, y=yst+.005, label=pval, angle=90), size=6)+
  scale_y_continuous(limits=c(0,0.1))+
  theme_bw()+
  theme(legend.text=element_text(size=14),
    axis.title.x=element_text(size=16),
    axis.title.y=element_text(size=16),
    axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
    axis.text.y=element_text(size=14))+
  xlab("Mutation Type")+
  ylab(nrsq)
# ggsave(paste0(parentdir, "/images/Lv7v5v3_rsq.png"), width=12, height=6)

Lv7v5v3_full <- ggplotGrob(Lv7v5v3_full)
Lv7v5v3 <- ggplotGrob(Lv7v5v3)

g <- arrangeGrob(Lv7v5v3_full, Lv7v5v3, nrow=1, widths=c(1,2))
ggsave(paste0(parentdir, "/images/Lv7v5v3_rsq_combined.svg"),
	width=12, height=6, g)

# Plot pseudo-r^2 for 7-mers vs logit
rsqdatL <- fulldat %>%
  filter(mod=="7-mers" | mod=="7-mers+features") %>%
  filter(group=="FULL") %>%
  mutate(Category =
      factor(plyr::mapvalues(category, orderedcats2, orderedcats2),
				levels=orderedcats2))

ggplot(rsqdatL)+
  geom_bar(aes(x=Category, y=rsq, fill=mod),
		stat="identity", position="dodge")+
  scale_fill_manual("Model", values=cbbPalette[c(7,8)])+
	theme_nrsq()+
	ylab(nrsq)+
  xlab("Mutation Type")
ggsave(paste0(parentdir, "/images/7v5v3_rsqL.png"), width=12, height=6)

# Plot pseudo-r^2 for ERVs vs Common
rsqdatEC <- fulldat %>%
	filter(mod=="AV" | mod=="Common" | mod=="ERVs") %>%
  filter(group=="FULL") %>%
  mutate(Category =
      factor(plyr::mapvalues(category, orderedcats2, orderedcats2),
				levels=orderedcats2))
ggplot(rsqdatEC)+
	geom_bar(aes(x=Category, y=rsq, fill=mod),
		stat="identity", position="dodge")+
	scale_fill_manual("Model", values=c("grey30", cbbPalette[c(3,7)]))+
	theme_bw()+
	theme(legend.text=element_text(size=14),
		axis.title.x=element_text(size=16),
		axis.title.y=element_text(size=16),
		axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
		axis.text.y=element_text(size=14))+
	xlab("Mutation Type")+
  ylab(nrsq)
ggsave(paste0(parentdir, "/images/AvCvS_rsq.png"), width=12, height=6)

oldmodnames <- unique(fulldat$mod)

newmodnames <- c("3-mers", "5-mers", "7-mers", "7-mers+features",
	"7-mers (downsampled BRIDGES ERVs)",
	"7-mers (BRIDGES MAC10+ variants)",
	"7-mers (1KG EUR intergenic variants)")

newmodnamesord <- c("3-mers", "5-mers", "7-mers",
	"7-mers (downsampled BRIDGES ERVs)",
	"7-mers (BRIDGES MAC10+ variants)",
	"7-mers (1KG EUR intergenic variants)",
	"7-mers+features")

rsqdatFULL <- fulldat %>%
  filter(group=="FULL") %>%
  mutate(Category =
      factor(plyr::mapvalues(category, orderedcats2, orderedcats2),
				levels=orderedcats2)) %>%
	mutate(mod=
		factor(plyr::mapvalues(mod, oldmodnames, newmodnames),
			levels=newmodnamesord))


alldat <- combineddat[2:8,] %>%
	mutate(mod=factor(plyr::mapvalues(mod, oldmodnames, newmodnames),
		levels=newmodnamesord))

all_full <- ggplot(alldat)+
  geom_bar(aes(x=category, y=rsq, fill=mod),
		stat="identity", position="dodge")+
	scale_fill_manual("Model", values=c(iwhPalette[c(3:9)]))+
  theme_nrsq(legend.position="none")+
  ylab(nrsq)

all <- ggplot(rsqdatFULL)+
  geom_bar(aes(x=Category, y=rsq, fill=mod),
		stat="identity", position="dodge")+
	scale_fill_manual("Model", values=c(iwhPalette[c(3:9)]))+
  theme_bw()+
  theme(
      legend.text=element_text(size=14),
			legend.position=c(0.4,0.7),
			legend.background = element_rect(colour = "black"),
      axis.title.x=element_text(size=16),
      axis.title.y=element_text(size=16),
    axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
      axis.text.y=element_text(size=14))+
  xlab("Mutation Type")+
  ylab(nrsq)

all_full <- ggplotGrob(all_full)
all <- ggplotGrob(all)
g <- arrangeGrob(all_full, all, nrow=1, widths=c(1,2))
ggsave(paste0(parentdir, "/images/all_rsq_combined.png"), width=12, height=6, g)

evcdat <- rsqdatFULL %>%
	filter(grepl("BRIDGES", mod))
EvC <- ggplot(evcdat)+
  geom_bar(aes(x=Category, y=rsq, fill=mod),
		stat="identity", position="dodge")+
  # scale_fill_manual("Model", values=brewer.pal(8, "Set3")[5:6])+
	scale_fill_manual("Model",
		values=c(iwhPalette[c(6,7)]),
		guide = guide_legend(nrow=2))+
  theme_bw()+
  theme(
      legend.text=element_text(size=14),
			legend.position="bottom",
      axis.title.x=element_text(size=16),
      axis.title.y=element_text(size=16),
    axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
      axis.text.y=element_text(size=14))+
  xlab("Mutation Type")+
  ylab(nrsq)
# ggsave(paste0(parentdir, "/images/EvC_rsq.png"), width=12, height=6)

EvC_full_dat <- combineddat %>%
  # filter(mod=="AV" | mod=="Common" | mod=="ERVs") %>%
	filter(mod=="Common" | mod=="ERVs") %>%
	mutate(mod=plyr::mapvalues(mod, c("Common", "ERVs"),
		c("7-mers (BRIDGES MAC10+ variants)", "7-mers (downsampled BRIDGES ERVs)"))) %>%
# combineddat %>%
  # filter(mod=="AV" | mod=="Common" | mod=="ERVs") %>%
	# filter(mod=="Common" | mod=="ERVs") %>%
	mutate(mod=factor(plyr::mapvalues(mod, oldmodnames, newmodnames),
		levels=newmodnamesord)) %>%
	filter(grepl("BRIDGES", mod))

EvC_full <- ggplot(EvC_full_dat)+
  geom_bar(aes(x=category, y=rsq, fill=mod), stat="identity", position="dodge")+
  # scale_fill_manual("Model", values=brewer.pal(8, "Set3")[5:6])+
	scale_fill_manual("Model", values=c(iwhPalette[c(6,7)]))+
	theme_nrsq()+
  ylab(nrsq)
# ggsave(paste0(parentdir, "/images/EvC_rsq_full.png"), width=8, height=6)

EvC_full <- ggplotGrob(EvC_full)
EvC <- ggplotGrob(EvC)
g <- arrangeGrob(EvC_full, EvC, nrow=1, widths=c(1,2))
ggsave(paste0(parentdir, "/images/EvC_rsq_combined.png"), width=12, height=8, g)
