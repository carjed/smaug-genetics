##############################################################################
# Dump for old code, comparing normalized mean rates of DNMs
##############################################################################


require(ggplot2)
require(dplyr)
require(tidyr)

writecats <- 0

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


##############################################################################
# Temp code for pred curves
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

##############################################################################
# Read and process data
##############################################################################
{
	# Read data
	cat("Reading data...\n")
	chrpf <- read.table("/net/bipolar/jedidiah/mutation/output/rocdat.7bp.2.txt", header=F, stringsAsFactors=F)
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
	}
	# Write DNM data per category with columns ID, CHR, POS, Category
	# Necessary if rocdat.7bp.txt not yet created
	writecats <- 0
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

# Recombine data
cat("Creating combined data...\n")
chrp <- rbind(chrpfdnm, chrpfa) %>%
  group_by(Category.x) %>%
  mutate(prop=cumsum(OBS)/sum(OBS)) %>%
  arrange(MU, prop)

# Add columns for 5-mer, 3-mer, and 1-mer rates
rates5 <- read.table("/net/bipolar/jedidiah/mutation/output/5bp_1000k_rates.txt", header=T, stringsAsFactors=F)
r5r <- data.frame(Category.x=gsub("cpg_", "", rates5$Category2),
	SEQ5=substr(rates5$Sequence, 1, 5),
	MU_5=rates5$rel_prop, stringsAsFactors=F)
chrp$SEQ5 <- substr(chrp$SEQ, 2, 6)
chrp <- merge(chrp, r5r, by=c("Category.x", "SEQ5"), all.x=T)

rates3 <- read.table("/net/bipolar/jedidiah/mutation/output/3bp_1000k_rates.txt", header=T, stringsAsFactors=F)
r3r <- data.frame(Category.x=gsub("cpg_", "", rates3$Category2),
	SEQ3=substr(rates3$Sequence, 1, 3),
	MU_3=rates3$rel_prop, stringsAsFactors=F)
chrp$SEQ3 <- substr(chrp$SEQ, 3, 5)
chrp <- merge(chrp, r3r, by=c("Category.x", "SEQ3"))

rates1 <- rates3 %>%
	mutate(Category.x=gsub("cpg_", "", Category2)) %>%
	group_by(Category.x) %>%
	summarise(n=sum(num), COUNT=sum(COUNT), MU_1=n/COUNT) %>%
	dplyr::select(Category.x, MU_1)
chrp <- merge(chrp, rates1, by=c("Category.x"))

r7s <- r5m %>%
	dplyr::select(Category.x=Category,
		SEQ=Sequence,
		MU_7S=rel_prop,
		MU_7P=common_rel_prop)
chrp <- merge(chrp, r7s, by=c("Category.x", "SEQ"))


# Data for comparing ERV models
chrp1 <- chrp %>%
	dplyr::select(Category.x, CHR, POS, OBS, SEQ, SEQ3, SEQ5,
		MU_L=MU, MU_7=MU_S, MU_5, MU_3, MU_1) %>%
	gather(Model, value, MU_L:MU_1) %>%
	group_by(Category.x) %>%
	mutate(zscore=scale(value, center=T, scale=T))
	# mutate_each(funs(scale(., center=T, scale=T)),
	# 	c(MU, MU_S, MU_5, MU_3, MU_1, MU_7S, MU_7P))
	# mutate(zscore=scale(MU, center=T, scale=T))

# Data for comparing ERVs/polymorphisms
chrp2 <- chrp %>%
	dplyr::select(Category.x, CHR, POS, OBS, SEQ, MU_7S, MU_7P) %>%
	# gather(Model, value, MU_7:MU_1) %>%
	group_by(Category.x) %>%
	# mutate(zscore=scale(value, center=T, scale=T))
	# mutate_each(funs(scale(., center=T, scale=T)),
		# c(MU_7S, MU_7P)) %>%
	gather(Model, MU, MU_7S:MU_7P)
	# mutate(zscore=scale(MU, center=T, scale=T))

chrp %>%
	group_by(Category.x) %>%
	summarise(MS=mean(MU_7S),
		MP=mean(MU_7P),
		pval=t.test(MU_7S, MU_7P)$p.value)

chrpfdnm2 <- chrpfdnm %>% dplyr::select(-c(MU_7S, MU_7P))
r7s <- r5m %>%
	filter(prop_diff5>1.25 | prop_diff5<0.8) %>%
	dplyr::select(Category.x=Category,
		SEQ=Sequence,
		MU_7S=rel_prop,
		MU_7P=common_rel_prop)
chrpfdnm2 <- merge(chrpfdnm2, r7s, by=c("Category.x", "SEQ"))
chrpfdnm2 %>%
	group_by(Category.x) %>%
	summarise(MS=mean(MU_7S),
		MP=mean(MU_7P),
		pval=t.test(MU_7S, MU_7P)$p.value)


chrpobs <- chrp %>% filter(OBS==1)

spectra3 <- r5m %>%
	# filter(prop_diff5>1.25 | prop_diff5<0.8) %>%
	mutate(SEQ3=substr(Sequence,3,5)) %>%
	group_by(Category, SEQ3) %>%
	summarise(s1=sum(num.x)/12088022, s2=sum(num.y)/12088017,
		r1=sum(num.x)/sum(COUNT.x), r2=sum(num.y)/sum(COUNT.y),
		COUNT=sum(COUNT.x)) %>% ungroup()

dnm_spectra3 <- chrpfdnm %>%
	group_by(Category.x, SEQ3) %>%
	summarise(n=n(), s=n()/nrow(chrpfdnm)) %>%
	ungroup() %>%
	dplyr::select(n,s)

spectra3 <- cbind(spectra3, dnm_spectra3)
spectra3$r <- spectra3$n/spectra3$COUNT

s3a <- spectra3 %>%
	group_by(Category) %>%
	summarise(cor1=cor(r1, r, method="spearman"),
		cor2=cor(r2, r, method="spearman"))

# Data for comparing ERV models
chrp1 <- chrp %>%
	dplyr::select(Category.x, CHR, POS, OBS, SEQ, SEQ3, SEQ5,
		MU_L=MU, MU_7=MU_S, MU_5, MU_3, MU_1) %>%
	gather(Model, value, MU_L:MU_1) %>%
	group_by(Category.x) %>%
	mutate(zscore=scale(value, center=T, scale=T))

# Data for comparing ERVs/polymorphisms
chrp2 <- chrp %>%
	dplyr::select(Category.x, CHR, POS, OBS, SEQ, MU_7S, MU_7P) %>%
	group_by(Category.x) %>%
	mutate_each(funs(scale(., center=T, scale=T)),
		c(MU_7S, MU_7P)) %>%
	gather(Model, zscore, MU_7S:MU_7P)

require(broom)
chrp1sum <- chrp1 %>%
  ungroup() %>%
  filter(OBS==1) %>%
  mutate(zscore=abs(jitter(zscore, 0.0001))) %>%
  group_by(Category.x, Model) %>%
  dplyr::do(broom::tidy(t.test(.$zscore)))
  # summarise(zscore_m=mean(zscore),
  #   low=t.test(zscore)$conf.int[1],
  #   high=t.test(zscore)$conf.int[2])

chrp1sum %>%
  filter(Model!="MU_L") %>%
  ggplot(aes(x=Model, y=estimate, group=Model, colour=Category.x))+
  geom_point()+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, width=0.5), position="dodge")+
  facet_wrap(~Category.x, scales="free")+
  ylab("Z-score")+
  scale_colour_manual("Mutation Type", values=cols[1:6])+
  theme_classic()+
  theme(legend.position="none",
    axis.title.x=element_blank())
ggsave("/net/bipolar/jedidiah/mutation/images/qcomp1.png", width=12, height=6)

chrp1sum %>%
  filter(Model %in% c("MU_L", "MU_7")) %>%
  ggplot(aes(x=Model, y=estimate, group=Model, colour=Category.x))+
  geom_point()+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, width=0.5), position="dodge")+
  facet_wrap(~Category.x, scales="free")+
  ylab("Z-score")+
  scale_colour_manual("Mutation Type", values=cols[1:6])+
  theme_classic()+
  theme(legend.position="none",
    axis.title.x=element_blank())
ggsave("/net/bipolar/jedidiah/mutation/images/qcomp2.png", width=12, height=6)

chrp1 %>%
  # group_by(Model) %>%
    filter(OBS==1) %>%
    filter(Model %in% c("MU_L", "MU_7")) %>%
    # dplyr::select(-c(value, SEQ3, SEQ5)) %>%
    # mutate(row=row_number()) %>%
    # spread(Model, zscore) %>%
    # mutate(diff=MU_L-MU_7) %>%
ggplot(aes(x=zscore, fill=Model, group=Model))+
  # geom_point()+
  geom_histogram(alpha=0.3, position="identity", binwidth=0.2)+
  # geom_violin()+
  # geom_hline(yintercept=0, linetype="dashed")+
  # geom_boxplot()+
  # geom_line()+
  # geom_errorbar(aes(ymin=conf.low, ymax=conf.high, width=0.5), position="dodge")+
  facet_wrap(~Category.x, scales="free")+
  # xlab("Z-score")+
  # scale_colour_manual("Mutation Type", values=cols[1:6])+
  theme_classic()+
  theme(#legend.position="none",
    axis.title.x=element_blank())
ggsave("/net/bipolar/jedidiah/mutation/images/qcomp2a.png", width=12, height=6)

chrp2sum <- chrp2 %>%
  filter(OBS==1) %>%
  mutate(zscore=abs(jitter(zscore, 0.0001))) %>%
  group_by(Category.x, Model) %>%
  dplyr::do(broom::tidy(t.test(.$zscore)))

chrp2test <- chrp2 %>%
  group_by(Model) %>%
    filter(OBS==1) %>%
    mutate(row=row_number()) %>%
    spread(Model, zscore) %>%
    group_by(Category.x) %>%
    dplyr::do(broom::tidy(t.test(.$MU_7P, .$MU_7S)))

chrp2sum %>%
ggplot(aes(x=Model, y=estimate, group=Model, colour=Category.x))+
  geom_point()+
  # geom_boxplot()+
  # geom_line()+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, width=0.5), position="dodge")+
  facet_wrap(~Category.x, scales="free")+
  ylab("Z-score")+
  scale_colour_manual("Mutation Type", values=cols[1:6])+
  theme_classic()+
  theme(legend.position="none",
    axis.title.x=element_blank())
ggsave("/net/bipolar/jedidiah/mutation/images/qcomp3.png", width=12, height=6)

chrp2 %>%
filter(OBS==1) %>%
ggplot(aes( x=zscore, fill=Model, group=Model))+
  # geom_point()+
  geom_histogram(alpha=0.3, position="identity", binwidth=0.2)+
  # geom_violin()+
  # geom_boxplot()+
  # geom_line()+
  # geom_errorbar(aes(ymin=conf.low, ymax=conf.high, width=0.5), position="dodge")+
  facet_wrap(~Category.x, scales="free")+
  xlab("Z-score")+
  # scale_colour_manual("Mutation Type", values=cols[1:6])+
  theme_classic()+
  theme(#legend.position="none",
    axis.title.x=element_blank())
ggsave("/net/bipolar/jedidiah/mutation/images/qcomp3a.png", width=12, height=6)


maxc <- r5m %>%
  filter(prop_diff5<0.8 | prop_diff5>1.25) %>%
  dplyr::select(SEQ=Sequence, Category.x=Category)

chrp2maxc <- merge(maxc, chrp2, by=c("SEQ", "Category.x"))

chrp2a <- chrp2 %>%
  filter(OBS==1) %>%
  mutate(zscore=abs(jitter(zscore, 0.0001)))

chrp2maxc <- merge(maxc, chrp2a, by=c("SEQ", "Category.x"))

chrp2maxcsum <- chrp2maxc %>%
  # filter(OBS==1) %>%
  # mutate(zscore=abs(jitter(zscore, 0.0001))) %>%
  group_by(Category.x, Model) %>%
  dplyr::do(broom::tidy(t.test(.$zscore)))

chrp2maxcsum %>%
ggplot(aes(x=Model, y=estimate, group=Model, colour=Category.x))+
  geom_point()+
  # geom_boxplot()+
  # geom_line()+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, width=0.5), position="dodge")+
  facet_wrap(~Category.x, scales="free")+
  ylab("Z-score")+
  scale_colour_manual("Mutation Type", values=cols[1:6])+
  theme_classic()+
  theme(legend.position="none",
    axis.title.x=element_blank())
ggsave("/net/bipolar/jedidiah/mutation/images/qcomp4.png", width=12, height=6)


#######################
# from snippets.R
#######################


chrpobs <- chrp1 %>% filter(OBS==1)
commonscale <- chrpobs %>%
  group_by(Category.x) %>%
  summarise(scale=mean(MU_S)/mean(MU_C))
chrp_maxcobs <- chrp_maxc %>% filter(OBS==1)

normdat <- data.frame()
for(i in 1:1000){
  subdat <- chrpobs[sample(nrow(chrpobs), 1000),] %>%
  # subdat <- chrp_maxcobs[sample(nrow(chrp_maxcobs), 1000),] %>%
    group_by(Category.x) %>%
    # mutate(MU_C=rescale(MU_C, c(min(MU_S), max(MU_S)))) %>%
    summarise(n=n(),
      MU_1=sum(MU_1),
      MU_3=sum(MU_3),
      MU_5=sum(MU_5),
      MU_7=sum(MU_S),
      MU_L=sum(MU))
      # MU_C=sum(MU_C))
  # subdat <- merge(subdat, commonscale, by=c("Category.x"))
  # subdat$MU_C <- subdat$MU_C/subdat$scale
  normdat <- rbind(subdat, normdat)
}

normdat %>%
  group_by(Category.x) %>%   summarise(n=n(),
      MU_1m=mean(MU_1),
      MU_3m=mean(MU_3),
      MU_5m=mean(MU_5),
      MU_7m=mean(MU_7),
      MU_Lm=mean(MU_L),
      # MU_Cm=mean(MU_C),
      pval13=t.test(MU_1, MU_3)$p.value,
      ratio13=1-MU_1m/MU_3m,
      pval35=t.test(MU_3, MU_5)$p.value,
      ratio35=1-MU_3m/MU_5m,
      pval57=t.test(MU_5, MU_7)$p.value,
      ratio57=1-MU_5m/MU_7m,
      pval7L=t.test(MU_7, MU_L)$p.value,
      ratio7L=1-MU_7m/MU_Lm) %>% data.frame

normdat2 <- normdat %>%
  gather(Model, S, MU_1:MU_7) %>%
  mutate(Model=
    recode(Model, MU_1 = "1-mers", MU_3 = "3-mers", MU_5 = "5-mers", MU_7 = "7-mers"))

ggplot(normdat2, aes(x=Model, y=S))+
  geom_boxplot()+
  facet_wrap(~Category.x, scales="free")+
  theme_classic()
ggsave("/net/bipolar/jedidiah/mutation/images/S_box_K-mers.png")

normdat2L <- normdat %>%
  gather(Model, S, MU_7:MU_L) %>%
  mutate(Model=
    recode(Model, MU_3 = "3-mers", MU_5 = "5-mers", MU_7 = "7-mers"))

ggplot(normdat2L, aes(x=Model, y=S))+
  geom_boxplot()+
  facet_wrap(~Category.x, scales="free")+
  theme_classic()
ggsave("/net/bipolar/jedidiah/mutation/images/S_box_logit_vs_7mers.png")

qdat <- chrp %>%
  filter(OBS==1) %>%
  group_by(Category.x) %>%
  summarise(n=n(),
    MU_1m=mean(MU_1), low1=MU_1m, high1=MU_1m,
    MU_3m=mean(MU_3), low3=t.test(MU_3)$conf.int[1], high3=t.test(MU_3)$conf.int[2],
    MU_5m=mean(MU_5), low5=t.test(MU_5)$conf.int[1], high5=t.test(MU_5)$conf.int[2],
    MU_7m=mean(MU_S), low7=t.test(MU_S)$conf.int[1], high7=t.test(MU_S)$conf.int[2],
    MU_Lm=mean(MU), lowL=t.test(MU)$conf.int[1], highL=t.test(MU)$conf.int[2]) %>%
    gather(Model, S, MU_1m:high1) %>%
    mutate(Model=
      recode(Model,
        MU_1m = "1-mers",
        MU_3m = "3-mers",
        MU_5m = "5-mers",
        MU_7m = "7-mers",
        MU_Lm = "7-mers+features"))

chrpobs %>% group_by(Category.x) %>%   summarise(n=n(),
    # MU_1m=mean(MU_1),
    # MU_3m=mean(MU_3),
    # MU_5m=mean(MU_5),
    # MU_7m=mean(MU_7),
    # MU_Lm=mean(MU_L),
    # MU_Cm=mean(MU_C),
    pval13=t.test(MU_1, MU_3)$p.value,
    # ratio13=1-MU_1m/MU_3m,
    pval35=t.test(MU_3, MU_5)$p.value,
    # ratio35=1-MU_3m/MU_5m,
    pval57=t.test(MU_5, MU_S)$p.value,
    # ratio57=1-MU_5m/MU_7m,
    pval7L=t.test(MU_S, MU)$p.value)
    # ratio7L=1-MU_7m/MU_Lm)

qdat <- chrpobs %>%
  ungroup() %>%
  mutate(MU_1=-abs(jitter(MU_1, 0.001))) %>%
  dplyr::select(Category.x, MU_1, MU_3, MU_5, MU_7=MU_S, MU_L=MU, MU_7S, MU_7P) %>%
  group_by(Category.x) %>%
  gather(Group, MU, c(MU_1, MU_3, MU_5, MU_7, MU_L, MU_7S, MU_7P)) %>%
  group_by(Category.x, Group) %>%
  summarise(n=n(),
    mean=mean(MU),
    low=t.test(MU)$conf.int[1],
    high=t.test(MU)$conf.int[2]) %>%
  mutate(Model= recode(Group,
    MU_1 = "1-mers",
    MU_3 = "3-mers",
    MU_5 = "5-mers",
    MU_7 = "7-mers",
    MU_L = "7-mers+features",
    MU_7S = "Rs",
    MU_7P = "Rp"))
#
# qdat %>%
#   # filter(!grepl("features", Model)) %>%
#   ungroup() %>%
require(broom)
chrp1 %>%
  filter(OBS==1) %>%
  # data.frame() %>%
  ungroup() %>%
  # filter(Model!="MU_L") %>%
  # mutate(zscore=-abs(jitter(zscore, 0.0001))) %>%
  # group_by(Category.x, Model) %>%
  # dplyr::do(broom::tidy(t.test(.$zscore)))
  # summarise(zscore_m=mean(zscore),
  #   low=t.test(zscore)$conf.int[1],
  #   high=t.test(zscore)$conf.int[2]) %>%
  filter(Model!="MU_L") %>%
  ggplot(aes(x=Model, y=value, group=Model, colour=Category.x))+
  # geom_point()+
  # scale_y_log10()+
  geom_boxplot()+
  # geom_line()+
  # geom_errorbar(aes(ymin=low, ymax=high, width=0.5), position="dodge")+
  facet_wrap(~Category.x, scales="free")+
  ylab("Relative rate")+
  scale_colour_manual("Mutation Type", values=cols[c(4,1,7,2,5,8)])+
  theme_classic()+
  theme(legend.position="none",
    axis.title.x=element_blank())
ggsave("/net/bipolar/jedidiah/mutation/images/qcomp1.png", width=12, height=6)

chrp1 %>%
  filter(OBS==1) %>%
  filter(Model %in% c("MU_L", "MU_7")) %>%
  # group_by(Category.x, Model) %>%
  # summarise(zscore_m=mean(zscore),
  #   low=t.test(zscore)$conf.int[1],
  #   high=t.test(zscore)$conf.int[2]) %>%
  ggplot(aes(x=Model, y=value, group=Model, colour=Category.x))+
  geom_point()+
  geom_boxplot()+
  scale_y_log10()+
  # geom_line()+
  # geom_errorbar(aes(ymin=low, ymax=high, width=0.5), position="dodge")+
  facet_wrap(~Category.x)+
  ylab("Z-score")+
  scale_colour_manual("Mutation Type", values=cols[c(4,1,7,2,5,8)])+
  theme_classic()+
  theme(legend.position="none",
    axis.title.x=element_blank())
ggsave("/net/bipolar/jedidiah/mutation/images/qcomp2.png", width=12, height=6)

chrp2 %>%
  filter(OBS==1) %>%
  group_by(Category.x, Model) %>%
  summarise(zscore_m=mean(zscore),
    low=t.test(zscore)$conf.int[1],
    high=t.test(zscore)$conf.int[2]) %>%
ggplot(aes(x=Model, y=zscore_m, group=Model, colour=Category.x))+
  geom_point()+
  # geom_boxplot()+
  # geom_line()+
  geom_errorbar(aes(ymin=low, ymax=high, width=0.5), position="dodge")+
  facet_wrap(~Category.x, scales="free")+
  ylab("Z-score")+
  scale_colour_manual("Mutation Type", values=cols[1:6])+
  theme_classic()+
  theme(legend.position="none",
    axis.title.x=element_blank())
ggsave("/net/bipolar/jedidiah/mutation/images/qcomp3.png", width=12, height=6)

chrp2 %>%
  filter(OBS==1) %>%
  group_by(Category.x, Model) %>%
  ungroup() %>%
  mutate(Model = plyr::mapvalues(Model, c("MU_7P", "MU_7S"), c("Common variants", "Singletons"))) %>%
  # summarise(zscore_m=mean(zscore),
  #   low=t.test(zscore)$conf.int[1],
  #   high=t.test(zscore)$conf.int[2]) %>%
ggplot(aes(x=Model, y=MU, colour=Category.x))+
  geom_boxplot()+#aes(alpha=factor(OBS)))+
  scale_y_log10()+
  # geom_boxplot()+
  # geom_line()+
  # geom_errorbar(aes(ymin=low, ymax=high, width=0.5), position="dodge")+
  facet_wrap(~Category.x)+
  ylab("Relative mutation rate")+
  scale_colour_manual("Mutation Type", values=cols[c(4,1,7,2,5,8)])+
  theme_bw()+
  theme(#legend.position="none",
    axis.title.x=element_blank())
ggsave("/net/bipolar/jedidiah/mutation/images/qcomp3a.png", width=12, height=6)

qdat1 <- chrp %>%
  filter(OBS==1) %>%
  # group_by(Category.x) %>%
  summarise(n=n(),
    MU_3=n/sum(MU_3),
    MU_5=n/sum(MU_5),
    MU_7=n/sum(MU_S),
    MU_L=n/sum(MU)) %>%
    gather(Model, Q, MU_3:MU_L)

ggplot(qdat1, aes(x=Model, y=Q, group=n))+
  geom_point()+
  geom_line()+
  # facet_wrap(~Category.x, scales="free")+
  theme_bw()+
  theme(legend.position="none")
ggsave("/net/bipolar/jedidiah/mutation/images/qcomp1.png")
