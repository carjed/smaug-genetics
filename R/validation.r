##############################################################################
# Read and process data
##############################################################################
cat("Reading data...\n")
chrpf <- read.table(paste0(parentdir,
		"/output/rocdat.7bp.2.txt"),
	header=F, stringsAsFactors=F)
names(chrpf) <- c("CHR", "POS", "BIN", "MU", "OBS", "Category", "SEQ",
	"MU_C", "MU_S", "MU_A")

# Remove sites with mu=0
chrpf <- chrpf[chrpf$MU>0,]

# Read DNMs
cat("Reading DNMs...\n")
source("read_dnms.r")

# Duplicate data, merge with DNMs to get ID
cat("Annotating with ID...\n")
chrpf <- merge(chrpf, dnms_full, by=c("CHR", "POS"), all.x=T)
chrpf$ID[is.na(chrpf$ID)] <- "all"

# Subset to non-DNMs and DNMs
cat("Splitting by DNM status...\n")
chrpfa <- chrpf[chrpf$ID=="all",]
chrpfa <- chrpfa[sample(nrow(chrpfa), 1000000),]
chrpfdnm <- chrpf[chrpf$ID!="all",]

# Recombine data
cat("Creating combined data...\n")
chrp <- rbind(chrpfdnm[,1:12], chrpfa) %>%
  group_by(Category.x) %>%
  mutate(prop=cumsum(OBS)/sum(OBS)) %>%
  arrange(MU, prop)

# Output 3-mer rates calculated on DNMs
dnm_rates<-0
if(dnm_rates){
	dnm_agg <- chrpfdnm %>%
		group_by(Category.x=as.character(Category.x), SEQ=substr(SEQ, 3, 5)) %>%
		summarise(ndnm=n())

	rates3 <- read.table(paste0(parentdir, "/output/3bp_1000k_rates.txt"),
		header=T, stringsAsFactors=F)
	rates3$Category.x <- gsub("cpg_", "", rates3$Category2)
	rates3$SEQ <- substr(rates3$Sequence, 1, 3)
	r3m <- merge(dnm_agg, rates3, by=c("Category.x", "SEQ"))

	rates_3 <- r3m %>%
		mutate(rel_prop3=ndnm/COUNT) %>%
	  dplyr::select(Sequence, Category.x, rel_prop3) %>%
	  spread(Category.x, rel_prop3)

	rates_3[is.na(rates_3)] <- 0

	write.table(rates_3, paste0(parentdir, "/dnm_3bp_rates.txt"),
		col.names=T, row.names=F, quote=F, sep="\t")
}

##############################################################################
# Add columns for 5-mer, 3-mer, and 1-mer rates
##############################################################################
rates5 <- ratelist[[3]] %>%
	mutate(SEQ5=substr(Motif, 1, 5)) %>%
	dplyr::select(Category.x=Type, SEQ5, MU_5=ERV_rel_rate)

chrp$SEQ5 <- substr(chrp$SEQ, 2, 6)
chrp <- merge(chrp, rates5, by=c("Category.x", "SEQ5"), all.x=T)

rates3 <- ratelist[[2]] %>%
	mutate(SEQ3=substr(Motif, 1, 3)) %>%
	dplyr::select(Category.x=Type, SEQ3, MU_3=ERV_rel_rate)

chrp$SEQ3 <- substr(chrp$SEQ, 3, 5)
chrp <- merge(chrp, rates3, by=c("Category.x", "SEQ3"))

rates1 <- ratelist[[1]] %>%
	dplyr::select(Category.x=Type, MU_1=ERV_rel_rate)

chrp <- merge(chrp, rates1, by=c("Category.x"))
#!!!
##############################################################################
# Add columns for downsampled ERV 7-mers and Common 7-mers
# Must have run sing_vs_com.r
##############################################################################
# r7s <- r5m %>%
# 	dplyr::select(Category.x=Category,
# 		SEQ=Sequence,
# 		MU_7S=rel_prop,
# 		MU_7P=common_rel_prop)
# chrp <- merge(chrp, r7s, by=c("Category.x", "SEQ"))

runDNMLogit<-function(data, group){
	if(nrow(data)>1e6){
		logmod1<-glm(OBS~MU_1, data=data, family=binomial())
		logmod3<-glm(OBS~MU_1+resid3, data=data, family=binomial())
		logmod5<-glm(OBS~MU_1+resid3+resid5, data=data, family=binomial())
		logmod7<-glm(OBS~MU_1+resid3+resid5+resid7, data=data, family=binomial())
		logmodL<-glm(OBS~MU_1+resid3+resid5+resid7+residL, data=data, family=binomial())

		# logmodSa<-glm(OBS~MU_7S, data=data, family=binomial())
		# logmodPa<-glm(OBS~MU_7P, data=data, family=binomial())
		logmod1a<-glm(OBS~MU_1, data=data, family=binomial())
		logmod3a<-glm(OBS~MU_3, data=data, family=binomial())
		logmod5a<-glm(OBS~MU_5, data=data, family=binomial())
		logmod7a<-glm(OBS~MU_S, data=data, family=binomial())
		logmodLa<-glm(OBS~MU, data=data, family=binomial())
		logmodAa<-glm(OBS~MU_A, data=data, family=binomial())

		# outdat<-list(logmod1, logmod3, logmod5, logmod7, logmodL,
		# 	logmodS, logmodP, logmodA,
		outdat <- list(logmod1a, logmod3a, logmod5a, logmod7a, logmodLa, logmodAa)
			# logmodSa, logmodPa, )
	} else {
		logmod3<-glm(OBS~MU_3, data=data, family=binomial())
		logmod5<-glm(OBS~MU_3+resid5, data=data, family=binomial())
		logmod7<-glm(OBS~MU_3+resid5+resid7, data=data, family=binomial())
		logmodL<-glm(OBS~MU_3+resid5+resid7+residL, data=data, family=binomial())

		# logmodS<-glm(OBS~MU_7S, data=data, family=binomial())
		# logmodP<-glm(OBS~MU_7P, data=data, family=binomial())
		# logmod1a<-glm(OBS~MU_1, data=data, family=binomial())
		logmod3a<-glm(OBS~MU_3, data=data, family=binomial())
		logmod5a<-glm(OBS~MU_5, data=data, family=binomial())
		logmod7a<-glm(OBS~MU_S, data=data, family=binomial())
		logmodLa<-glm(OBS~MU, data=data, family=binomial())
		logmodAa<-glm(OBS~MU_A, data=data, family=binomial())

		outdat <- list(logmod3a, logmod5a, logmod7a, logmodLa, logmodAa)

		# outdat<-list(logmod3, logmod5, logmod7, logmodL,
		# 	logmodS, logmodP, logmodA)
			# logmod3a, logmod5a, logmod7a, logmodLa)
	}

  return(outdat)
}

overall_dat <- chrp %>%
	mutate(Category=ifelse(substr(SEQ,adj+1,adj+2)=="CG",
									paste0("cpg_",Category.x), Category.x)) %>%
  # mutate(Category =
  #     plyr::mapvalues(Category, orderedcats1, orderedcats2)) %>%
  mutate(resid3=MU_3-MU_1,
		resid5=MU_5-MU_3,
		resid7=MU_S-MU_5,
		residL=MU-MU_S,
		resid5a=MU_5-MU_1,
		resid7a=MU_S-MU_1,
		residLa=MU-MU_1)
overall_models <- runDNMLogit(overall_dat, "FULL")

test13 <- lrtest(overall_models[[1]], overall_models[[2]])
test35 <- lrtest(overall_models[[2]], overall_models[[3]])
test57 <- lrtest(overall_models[[3]], overall_models[[4]])
test7L <- lrtest(overall_models[[4]], overall_models[[5]])
fulllist<-list(test13, test35, test57, test7L)

pvals <- lapply(fulllist, function(x) x[[5]][2])

rsq<-unlist(lapply(overall_models, function(x)
	NagelkerkeR2(x)))[seq(2, 2*length(overall_models), by=2)]

# rsq<-unlist(lapply(overall_models, function(x)
	# NagelkerkeR2(x)))[c(2,4,6,8,10,12,14,16,18,20,22,24,26)]

aic<-unlist(lapply(overall_models, function(x)
	AIC(x)))[1:length(overall_models)]
mod<-c("1-mers", "3-mers", "5-mers", "7-mers", "7-mers+features",
 "ERVs", "Common", "AV", "1-mers*", "3-mers*", "5-mers*", "7-mers*",
 "7-mers+features*", "ERVs*", "Common*", "AV*")
combineddat<-data.frame(group="FULL", category=categ, mod, rsq, aic)

rsq_full_dat<- combineddat %>%
  filter(mod %in% mod[1:5])
Lv7v5v3_full<-ggplot(rsq_full_dat)+
  geom_bar(aes(x=category, y=rsq, fill=mod), stat="identity", position="dodge")+
  # scale_fill_manual("Model", values=cbbPalette[c(1,5,4,6,7,8)])+
	scale_fill_manual("Model", values=c(iwhPalette[c(1,3,4,5,9)]))+
  theme_bw()+
  theme(
      legend.text=element_text(size=14),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=16),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=14))+
  ylab(expression(paste("Fraction of variance explained (Nagelkerke ", R^2, ")",
		sep="")))
# ggsave(paste0(parentdir, "/images/Lv7v5v3_rsq_full.png"), width=4, height=6)

# fulldat <- list()
fulldat <- data.frame()
lrtestdat <- data.frame()
for(i in 1:length(orderedcats)){
  categ <- orderedcats2[i]

  overall_dat <- chrp %>%
		mutate(Category=ifelse(substr(SEQ,adj+1,adj+2)=="CG",
										paste0("cpg_",Category.x), Category.x)) %>%
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
  theme(
      legend.text=element_text(size=14),
      axis.title.x=element_text(size=16),
      axis.title.y=element_text(size=16),
    axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
      axis.text.y=element_text(size=14))+
  xlab("Mutation Type")+
  ylab(expression(paste("Fraction of variance explained (Nagelkerke ", R^2, ")",
	 	sep="")))
# ggsave(paste0(parentdir, "/images/Lv7v5v3_rsq.png"), width=12, height=6)

Lv7v5v3_full<-ggplotGrob(Lv7v5v3_full)
Lv7v5v3<-ggplotGrob(Lv7v5v3)

g<-arrangeGrob(Lv7v5v3_full, Lv7v5v3, nrow=1, widths=c(1,2))
ggsave(paste0(parentdir, "/images/Lv7v5v3_rsq_combined.svg"), width=12, height=6, g)

# Plot pseudo-r^2 for 7-mers vs logit
rsqdatL<-fulldat %>%
  filter(mod=="7-mers" | mod=="7-mers+features") %>%
  filter(group=="FULL") %>%
  mutate(Category =
      factor(plyr::mapvalues(category, orderedcats2, orderedcats2),
				levels=orderedcats2))
ggplot(rsqdatL)+
  geom_bar(aes(x=Category, y=rsq, fill=mod), stat="identity", position="dodge")+
  scale_fill_manual("Model", values=cbbPalette[c(7,8)])+
  theme_bw()+
  theme(
      legend.text=element_text(size=14),
      axis.title.x=element_text(size=16),
      axis.title.y=element_text(size=16),
    axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
      axis.text.y=element_text(size=14))+
  xlab("Mutation Type")+
  ylab(expression(paste("Fraction of variance explained (Nagelkerke ", R^2, ")",
		sep="")))
ggsave(paste0(parentdir, "/images/7v5v3_rsqL.png"), width=12, height=6)

# Plot pseudo-r^2 for ERVs vs Common
rsqdatEC<-fulldat %>%
	filter(mod=="AV" | mod=="Common" | mod=="ERVs") %>%
  filter(group=="FULL") %>%
  mutate(Category =
      factor(plyr::mapvalues(category, orderedcats2, orderedcats2),
				levels=orderedcats2))
ggplot(rsqdatEC)+
	geom_bar(aes(x=Category, y=rsq, fill=mod), stat="identity", position="dodge")+
	scale_fill_manual("Model", values=c("grey30", cbbPalette[c(3,7)]))+
	theme_bw()+
	theme(
			legend.text=element_text(size=14),
			axis.title.x=element_text(size=16),
			axis.title.y=element_text(size=16),
		axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
			axis.text.y=element_text(size=14))+
	xlab("Mutation Type")+
	ylab(expression(paste("Fraction of variance explained (Nagelkerke ", R^2, ")",
		sep="")))
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

rsqdatFULL<-fulldat %>%
  filter(group=="FULL") %>%
  mutate(Category =
      factor(plyr::mapvalues(category, orderedcats2, orderedcats2),
				levels=orderedcats2)) %>%
	mutate(mod=
		factor(plyr::mapvalues(mod, oldmodnames, newmodnames),
			levels=newmodnamesord))


alldat<-combineddat[2:8,] %>%
	mutate(mod=factor(plyr::mapvalues(mod, oldmodnames, newmodnames),
		levels=newmodnamesord))

all_full<-ggplot(alldat)+
  geom_bar(aes(x=category, y=rsq, fill=mod), stat="identity", position="dodge")+
	scale_fill_manual("Model", values=c(iwhPalette[c(3:9)]))+
  theme_bw()+
  theme(
      legend.text=element_text(size=14),
			legend.position="none",
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=16),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=14))+
  ylab(expression(paste("Fraction of variance explained (Nagelkerke ", R^2, ")",
		sep="")))

all<-ggplot(rsqdatFULL)+
  geom_bar(aes(x=Category, y=rsq, fill=mod), stat="identity", position="dodge")+
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
  ylab(expression(paste("Fraction of variance explained (Nagelkerke ", R^2, ")",
		sep="")))

all_full<-ggplotGrob(all_full)
all<-ggplotGrob(all)
g<-arrangeGrob(all_full, all, nrow=1, widths=c(1,2))
ggsave(paste0(parentdir, "/images/all_rsq_combined.png"), width=12, height=6, g)

evcdat<- rsqdatFULL %>%
	filter(grepl("BRIDGES", mod))
EvC<-ggplot(evcdat)+
  geom_bar(aes(x=Category, y=rsq, fill=mod), stat="identity", position="dodge")+
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
  ylab(expression(paste("Fraction of variance explained (Nagelkerke ", R^2, ")",
		sep="")))
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

EvC_full<-ggplot(EvC_full_dat)+
  geom_bar(aes(x=category, y=rsq, fill=mod), stat="identity", position="dodge")+
  # scale_fill_manual("Model", values=brewer.pal(8, "Set3")[5:6])+
	scale_fill_manual("Model", values=c(iwhPalette[c(6,7)]))+
  theme_bw()+
  theme(
      legend.text=element_text(size=14),
			legend.position="none",
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=16),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=14))+
  ylab(expression(paste("Fraction of variance explained (Nagelkerke ", R^2, ")",
		sep="")))
# ggsave(paste0(parentdir, "/images/EvC_rsq_full.png"), width=8, height=6)

EvC_full<-ggplotGrob(EvC_full)
EvC<-ggplotGrob(EvC)
g<-arrangeGrob(EvC_full, EvC, nrow=1, widths=c(1,2))
ggsave(paste0(parentdir, "/images/EvC_rsq_combined.png"), width=12, height=8, g)

##############################################################################
# Get prediction curves under each model
##############################################################################
tmp <- 0
if(tmp){
	chrp_gfasc <- subByModel(chrp, "MU", "GFASC")
	head(chrp_gfasc) %>% data.frame()

	full_auc_dat <- subWrapper(chrp, sim=F)
	full_auc_dat %>%
		arrange(group, Category.x, ntile) %>%
		filter(Category.x=="GC_AT") %>%
		tail() %>%
		data.frame()

	outdat <- chrp %>%
		mutate_(mutmp="MU") %>%
		group_by(Category.x) %>%
		arrange(desc(mutmp)) %>%
		mutate(prop=cumsum(OBS)/sum(OBS)) %>%
		# dplyr::select(-mutmp) %>%
		# arrange_(mucol, "prop") %>%
		# mutate(mutmp=mucol, ntile=dense_rank(mutmp), group=groupname) %>%
		mutate(ntile=min_rank(desc(mutmp)), group="GFASC") %>%
		dplyr::select(-mutmp)
}

pred_curves <- 0
if(pred_curves){
	full_auc_dat <- subWrapper(chrp, sim=T)

	# Get table of AUC for each mutation type/model combo
	full_auc2 <- full_auc_dat %>%
	  group_by(group, Category.x) %>%
	  summarise(AUC=sum(prop)/n()) %>%
	  spread(Category.x, AUC)

	newrow <- c("Relative change in AUC", (full_auc2[5,-1]-0.5)/(full_auc2[3,-1]-0.5))
	full_auc2 <- data.frame(rbind(as.matrix(full_auc2), newrow))

	roc_plotdat <- full_auc_dat[sample(nrow(full_auc_dat), .1*nrow(full_auc_dat)),] %>%
		group_by(group, Category.x, ntile) %>%
		summarise(prop=min(prop), n=mean(n))

	# Plot everything together
	plotROC(roc_plotdat, paste0(parentdir, "/images/pseudo_roc_curves_all.png"))

	# Plot 3-mers only
	roc_plotdat_seq <- roc_plotdat %>%
		filter(grepl("^ERV 3", group))
	plotROC(roc_plotdat_seq, paste0(parentdir, "/images/pseudo_roc_curves_seq_3.png"))

	# Plot 3-mers vs 5-mers
	roc_plotdat_seq <- roc_plotdat %>%
		filter(grepl("^ERV 5|^ERV 3", group))
	plotROC(roc_plotdat_seq, paste0(parentdir, "/images/pseudo_roc_curves_seq_5_3.png"))

	# Plot 7-mers vs 5-mers vs 3-mers
	roc_plotdat_seq <- roc_plotdat %>%
		filter(grepl("^E", group))
	plotROC(roc_plotdat_seq, paste0(parentdir, "/images/pseudo_roc_curves_seq.png"))

	# Plot common vs A&V
	fig_a_dat <- roc_plotdat %>%
	  filter(grepl("^A|B", group))
	plotROC(fig_a_dat, paste0(parentdir, "/images/av_vs_common_roc.png"))

	# Plot common vs ERVs
	fig_b_dat <- roc_plotdat %>%
	  filter(grepl("^E|B", group))
	plotROC(fig_b_dat, paste0(parentdir, "/images/common_vs_erv_roc.png"))

	# Plot ERVs vs logit
	fig_k_dat <- roc_plotdat %>%
		filter(grepl("^G|^E", group))
	plotROC(fig_k_dat, paste0(parentdir, "/images/erv_vs_logit_roc.png"))

	# Plot ERVs vs logit vs simulated max
	fig_k2_dat <- roc_plotdat %>%
		filter(grepl("^G|^E|^m", group))
	plotROC(fig_k2_dat, paste0(parentdir, "/images/erv_vs_logit_max_roc.png"))
}

##############################################################################
# ROC on 7-mers with max contrast between ERVs and polymorphisms
##############################################################################
maxc_roc<-0
if(maxc_roc){
	maxc <- read.table(paste0(parentdir, "/maxc_7bp.txt"), header=T, stringsAsFactors=F)
	maxc$Category <- gsub("cpg_", "", maxc$Category2)
	chrp_maxc <- chrp %>%
	  dplyr::filter(paste0(SEQ, ".", Category.x) %in%
									paste0(maxc$Sequence, ".", maxc$Category))

	maxc_auc_dat <- subWrapper(chrp_maxc, sim=F)

	maxc_auc <- maxc_auc_dat %>%
	  group_by(group, Category.x) %>%
		# filter(grepl("^B|^E", group)) %>%
	  summarise(AUC=sum(prop)/n()) %>%
	  spread(Category.x, AUC)

	newrow<-c("Relative change in AUC", (maxc_auc[2,-1]-0.5)/(maxc_auc[1,-1]-0.5))
	maxc_auc <- data.frame(rbind(as.matrix(maxc_auc), newrow))

	maxc_plotdat <- maxc_auc_dat[sample(nrow(maxc_auc_dat), .1*nrow(maxc_auc_dat)),]

	# fig_max_dat <- maxc_plotdat %>%
	# 	filter(grepl("^E|B", group))
	# Plot all models under max contrast
	fig_max_dat <- maxc_auc_dat %>%
		filter(grepl("^E|^B", group))
	plotROC(fig_max_dat, paste0(parentdir, "/images/pseudo_roc_curves_maxc.png"))
}

# Additional max contrast analyses, with features
maxc_extra1 <- 0
if(maxc_extra1){
	coefs <- read.table(paste0(parentdir, "/output/logmod_data/coefs/coefs_full.txt"),
		header=F, stringsAsFactors=F)

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
			# ungroup() %>%
			# mutate(Category.x=as.character(Category.x)) %>%
	    dplyr::filter((substr(SEQ,1,7) %in% x$Sequence) & (Category.x %in% x$Category))
	}

	maxc_feat <- data.frame(Sequence=plotdat$Sequence,
		Category=as.character(gsub("cpg_", "", plotdat$Category)))
	chrp_maxc_features<-matchMotifs(maxc_feat)
	maxc_auc_feat <- subWrapper(chrp_maxc_features, sim=F)
	cov_auc<-aucCalc(maxc_auc_feat)

	chrp_maxc_features <- lapply(maxc_feat, function(x) matchMotifs(x))

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
	ggsave(paste0(parentdir, "/images/test_bar.png"))

	cov_erv <- cov_auc_dat %>% filter(grepl("^E", group))
	cov_logit <- cov_auc_dat %>% filter(!grepl("^E", group))
	cov_logit$diff <- (cov_logit$AUC-0.5)/(cov_erv$AUC-0.5)-1
	# t.test(cov_logit$diff)

	cov_logit$dir <- ifelse(cov_logit$diff<=0, "lower", "higher")
	ggplot(cov_logit, aes(x=id, y=diff, colour=dir, fill=dir))+
		geom_bar(stat="identity", position="dodge")+
		scale_fill_brewer(palette="Set1")+
		scale_colour_brewer(palette="Set1")+
		# geom_point()+
	  facet_wrap(~Category.x)+
	  coord_cartesian(ylim=c(-0.15, 0.15))+
		coord_flip()+
		geom_hline(yintercept=0, linetype="dashed")+
		# xlab("AUC_[]")+
		theme_bw()
	ggsave(paste0(parentdir, "/images/diff_bar.png"), height=6, width=12)
}

maxc_extra2 <- 0
if(maxc_extra2){
	maxc_sum <- maxc_auc_dat %>%
		# mutate(SEQ3=substr(SEQ,3,5))%>%
		# group_by(Category.x) %>%
		group_by(Category.x, SEQ) %>%
		summarise(n=sum(OBS)) %>%
		dplyr::filter(n>0)

	# chrp_maxc2 <- merge(maxc_sum[,1:2], chrp_maxc, by=c("Category.x"))
	chrp_maxc$SEQ3 <- substr(chrp_maxc$SEQ, 3, 5)

	chrp_maxc2 <- merge(maxc_sum[,1:2], chrp_maxc, by=c("Category.x", "SEQ"))
	# chrp_maxc2 <- merge(maxc_sum[,1:2],
		# chrp_maxc[!grepl("^gonl", chrp_maxc$ID),], by=c("Category.x", "SEQ"))
	maxc_auc_dat2 <- subWrapper(chrp_maxc2, sim=F)

	maxc_auc2 <- maxc_auc_dat2 %>%
	  # group_by(group, Category.x) %>%
		group_by(group, Category.x, SEQ) %>%
	  summarise(AUC=1-sum(prop)/n()) #%>%
	  #spread(Category.x, AUC)

	# maxc_auc2 %>% arrange(Category.x, group) %>% data.frame()
	maxc_auc2 %>% arrange(Category.x, SEQ, group) %>% head(20)
	maxc_auc2 %>% group_by(group) %>% summarise(AUC=mean(AUC))
	maxc_test<-maxc_auc2 %>%
		group_by(group, Category.x) %>%
		# summarise(AUC=mean(AUC)) %>%
		arrange(Category.x) %>%
		spread(group, AUC) %>%
		data.frame()

	names(maxc_test) <- c("Category.x", "SEQ", "AV", "BP", "ERV", "GFASC")
	maxc_test$diff <- maxc_test$ERV-maxc_test$BP
	maxc_test %>% group_by(Category.x)%>% summarise(diff=mean(diff), n=n())
	# not by seq
	# maxc_auc <- maxc_auc_dat %>%
	#   group_by(group, Category.x) %>%
	#   summarise(AUC=1-sum(prop)/n()) %>%
	#   spread(Category.x, AUC)
}


##############################################################################
# ROC by GoNL individual
##############################################################################
gonl_ind <- 0
if(gonl_ind){
	nperm<-258
	# aucperm<-rep(0,nperm)
	# aucperm3<-aucperm

	ids <- unique(chrpfdnm$ID)[grepl("gonl", unique(chrpfdnm$ID))]
	numind <- length(ids)
	# aucind<-rep(0, numind)
	aucind3 <- aucind

	auc_ind<-data.frame()
	for(i in 1:nperm){
		### Run permutations
		# cat("Permuting AUC", "(", i, "of", nperm, ")...\n")
		ndnms<-round(rnorm(1, 42.7, 10.3), 0)
		nsamp<-1000000

		curid<-ids[i]
		### Get empirical AUC
		cat("Calculating AUC for", curid, "(", i, "of", numind, ")...\n")
		datid<-chrpfdnm %>% filter(ID==curid)

		tmpdat<-rbind(datid, chrpfa) %>% arrange(MU)
		tmp_auc_dat <- subWrapper(tmpdat, sim=F)

		tmp_auc <- tmp_auc_dat %>%
		  group_by(group, Category.x) %>%
		  summarise(AUC=1-sum(prop)/n()) %>%
		  spread(Category.x, AUC)
		tmp_auc$ID <- curid

		auc_ind<-rbind(auc_ind, data.frame(tmp_auc))
	}


	# ggplot(auc_ind, aes(x=group, y=))
	# mapply(plotROC, maxc_auc_feat, paste0(parentdir, "/images/", seq_along(max_auc_feat), "_roc.png"))
}
