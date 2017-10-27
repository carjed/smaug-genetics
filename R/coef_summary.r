# Read table of coefficients output from logit models
coefs <- read.table(paste0(parentdir, "/output/logmod_data/coefs/coefs_full.txt"),
  header=F, stringsAsFactors=F)
names(coefs) <- c("Cov", "Est", "SE", "Z", "pval", "Sequence", "Category")
coefs$Category <- ifelse(substr(coefs$Sequence, 4, 5)=="CG",
  paste0("cpg_", coefs$Category), coefs$Category)

coefsout <- coefs %>%
  mutate(Category = plyr::mapvalues(Category, orderedcats1, orderedcats2))
coefsout$Category <- factor(coefsout$Category, levels=orderedcats2)

write.table(coefsout, paste0(parentdir,
  "/output/logmod_data/coefs/supplementary_table_7.txt"),
  quote=F, col.names=T, row.names=F, sep="\t")

# ratefile <- paste0(parentdir, "/output/7bp_1000k_rates.txt")
# rates <- read.table(ratefile, header=T, stringsAsFactors=F)
rates <- ratelist[[4]]

rates2 <- rates %>%
  mutate(Category=Type) %>%
  mutate(Sequence=substr(Motif, 0, 7))

rates2$Category <- ifelse(substr(rates2$Sequence, 4, 5)=="CG",
  paste0("cpg_", rates2$Category), rates2$Category)

rates <- rates %>%
  mutate(Category=gsub("cpg_", "", Type)) %>%
  mutate(Sequence=substr(Motif, 0, 7))
motifdat <- rates %>%
  dplyr::select(Category, Sequence, rel_prop=ERV_rel_rate)

orderedcovs <- c("H3K9me3", "RR", "TIME", "H3K27me3",
  "H3K36me3", "DHS", "GC", "CpGI")



sigcoefs <- merge(coefs, rates2, by=c("Category", "Sequence")) %>%
  filter(nERVs>=20) %>%
  filter(!(Cov %in% c("(Intercept)", "DP"))) %>%
  group_by(Cov, Category) %>%
  # group_by(Category, Sequence) %>%
  mutate(qval=p.adjust(pval, method="fdr")) %>%
  filter(qval<0.05) %>%
  filter(!grepl("Intercept|DP", Cov)) %>%
  ungroup() %>%
  mutate(dir=ifelse(Est>0, "Up", "Down"),
    Est=ifelse(Cov=="RR", Est*10, Est),
    Est=ifelse(Cov=="GC", Est/10, Est)) %>%
  ungroup() %>%
  group_by(Cov, Category, dir) %>%
  mutate(n=n())

sigcoefs$Category2 <- plyr::mapvalues(sigcoefs$Category, orderedcats1, orderedcats2)
sigcoefs$Category2 <- factor(sigcoefs$Category2, levels=orderedcats2)

sig10txt <- sigcoefs %>%
  mutate(Estmax=ifelse(dir=="Up", 3, 0.5)) %>%
  mutate(Estmax=ifelse(Cov=="CpGI", ifelse(dir=="Up", 10, 0.125), Estmax),
    label=n) %>%
  filter(n>=10)

##############################################################################
# Theme for plotting coefficient distributions
##############################################################################
theme_coef <- function(base_size = 12, base_family = ""){
  theme_bw(base_size = base_size) %+replace%
  theme(legend.position="bottom",
    legend.title=element_blank(),
    strip.text=element_text(size=12),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.text.y=element_text(size=10),
    axis.title.y=element_text(size=12, angle=90))
}

##############################################################################
# Plot all effects (incl. NS)
##############################################################################
coefs <- coefs %>%
  mutate(Category = plyr::mapvalues(Category, orderedcats1, orderedcats2))

coefs$Category <- factor(coefs$Category, levels=orderedcats2)
coefs %>%
  filter(!grepl("Intercept", Cov)) %>%
  filter(exp(Est)<10) %>%
  ggplot(aes(fill=Category))+
    geom_hline(yintercept=1, linetype="dashed")+
    geom_violin(aes(x=Category, y=exp(Est), fill=Category))+
    scale_fill_manual(values=gp_cols, drop=FALSE)+
    scale_colour_manual(values=gp_cols, drop=FALSE)+
    facet_wrap(~Cov, ncol=4, scales="free_y")+
    ylab("odds ratio for mutability")+
    guides(fill = guide_legend(nrow = 3))+
    theme_coef()
ggsave(paste0(parentdir, "/images/coef_violin_full.pdf"), width=9, height=8)

##############################################################################
# Plot all significant effects
##############################################################################
sigcoefs %>%
  ggplot(aes(x=Category2, y=exp(Est), alpha=factor(dir), fill=Category2))+
    geom_hline(yintercept=1, linetype="dashed")+
    geom_violin(position="identity", scale="area")+
    scale_fill_manual(values=gp_cols, drop=FALSE)+
    scale_colour_manual(values=gp_cols, drop=FALSE)+
    scale_x_discrete(drop=FALSE)+
    scale_y_log10(breaks=c(0.5, 1, 2, 4))+
    scale_alpha_discrete(range = c(0.95, 0.96), guide=F)+
    facet_wrap(~Cov, ncol=4, drop=F)+
    ylab("odds ratio for mutability")+
    guides(fill = guide_legend(nrow = 3))+
    theme_coef()
ggsave(paste0(parentdir, "/images/coef_violin.pdf"), width=6.5, height=6.5)

##############################################################################
# Plot non-CpG effects (w/ >10 subtypes)
##############################################################################
non_cpg_covs <- c("H3K9me3", "RR", "TIME",
  "H3K36me3", "DHS", "GC")

non_cpg_plotdat <- sigcoefs %>%
  filter(Cov %in% non_cpg_covs) %>%
  filter(n>=10) %>%
  mutate(Est=ifelse(Cov=="TIME", -Est, Est))
non_cpg_plotdat$Cov <- factor(non_cpg_plotdat$Cov, levels=non_cpg_covs)
sig10txt1 <- sig10txt %>%
  group_by(Cov, dir, Category2) %>%
  summarise(label=max(label), Estmax=max(Estmax)) %>%
  ungroup() %>%
  mutate(Estmax=ifelse(Cov=="TIME", 1/Estmax, Estmax)) %>%
  filter(Cov %in% non_cpg_covs)
sig10txt1$Cov <- factor(sig10txt1$Cov, levels=non_cpg_covs)

non_cpg_plotdat %>%
  ggplot(aes(x=Category2, y=exp(Est), alpha=factor(dir), fill=Category2))+
    geom_hline(yintercept=1, linetype="dashed")+
    geom_text(data=sig10txt1,
      aes(x=Category2, y=Estmax, label=label, colour=Category2),
      vjust=1, size=4, angle=90)+
    geom_violin(position="identity", scale="area")+
    scale_fill_manual(values=gp_cols, drop=FALSE)+
    scale_colour_manual(values=gp_cols, drop=FALSE)+
    scale_x_discrete(drop=FALSE)+
    # scale_y_log10(breaks=c(.125, .25, 0.5, 1, 2, 4, 8))+
    scale_y_log10(breaks=c(0.5, 1, 2, 4))+
    scale_alpha_discrete(range = c(0.95, 0.96), guide=F)+
    facet_wrap(~Cov, ncol=3, drop=F)+
    ylab("odds ratio for mutability")+
    guides(fill = guide_legend(nrow = 3))+
    theme_coef()
ggsave(paste0(parentdir, "/images/coef_violin_10.pdf"), width=6.5, height=6.5)

non_cpg_plotdat %>%
  ggplot(aes(x=exp(Est), y=Category2, alpha=factor(dir), colour=Category2, fill=Category2))+
    geom_vline(yintercept=1, linetype="dashed")+
    geom_text(data=sig10txt1,
      aes(x=Estmax, y=Category2, label=label, colour=Category2),
      vjust=1, size=4, angle=90)+
    geom_bar(position="identity", scale="area")+
    scale_fill_manual(values=gp_cols, drop=FALSE)+
    scale_colour_manual(values=gp_cols, drop=FALSE)+
    scale_x_discrete(drop=FALSE)+
    # scale_y_log10(breaks=c(.125, .25, 0.5, 1, 2, 4, 8))+
    scale_y_log10(breaks=c(0.5, 1, 2, 4))+
    scale_alpha_discrete(range = c(0.95, 0.96), guide=F)+
    facet_wrap(~Cov, ncol=3, drop=F)+
    ylab("odds ratio for mutability")+
    guides(fill = guide_legend(nrow = 3))+
    theme_coef()
ggsave(paste0(parentdir, "/images/coef_hist_10.pdf"), width=6.5, height=6.5)


# non_cpg_plotdat %>%
#   ggplot(aes(x=exp(Est), y=Category, alpha=factor(dir), fill=Category))+
#     # geom_hline(yintercept=1, linetype="dashed")+
#     # geom_text(data=sig10txt1,
#     #   aes(x=Category2, y=Estmax, label=label, colour=Category2),
#     #   vjust=1, size=4, angle=90)+
#     geom_joy()+
#     scale_fill_manual(values=gp_cols, drop=FALSE)+
#     scale_colour_manual(values=gp_cols, drop=FALSE)+
#     # scale_x_discrete(drop=FALSE)+
#     # scale_y_log10(breaks=c(.125, .25, 0.5, 1, 2, 4, 8))+
#     # scale_x_log10(breaks=c(0.5, 1, 2, 4))+
#     # scale_alpha_discrete(range = c(0.95, 0.96), guide=F)+
#     facet_wrap(~Cov, ncol=3, drop=F)#+
#     # ylab("odds ratio for mutability")+
#     # guides(fill = guide_legend(nrow = 3))+
#     # theme_coef()
# ggsave(paste0(parentdir, "/images/coef_joy_10.pdf"), width=6.5, height=6.5)


##############################################################################
# Plot CpG island effects (w/ >10 subtypes)
##############################################################################
cpgi_covs <- c("CpGI", paste0("BLANK", 1:5))

cpgi_plotdat <- sigcoefs %>%
  filter(Cov %in% cpgi_covs) %>%
  filter(n>=10)
cpgi_plotdat$Cov <- factor(cpgi_plotdat$Cov, levels=cpgi_covs)
sig10txt2 <- sig10txt %>%
  group_by(Cov, dir, Category2) %>%
  filter(Cov %in% cpgi_covs) %>%
  summarise(label=max(label), Estmax=max(Estmax))
sig10txt2$Cov <- factor(sig10txt2$Cov, levels=cpgi_covs)

cpgi_plotdat %>%
  ggplot(aes(x=Category2, y=exp(Est), alpha=factor(dir), fill=Category2))+
    geom_hline(yintercept=1, linetype="dashed")+
    geom_text(data=sig10txt2,
      aes(x=Category2, y=Estmax, label=label, colour=Category2),
      vjust=1, size=4, angle=90)+
    geom_violin(position="identity", scale="area")+
    scale_fill_manual(values=gp_cols, drop=FALSE)+
    scale_colour_manual(values=gp_cols, drop=FALSE)+
    scale_x_discrete(drop=FALSE)+
    # scale_y_log10(breaks=c(.125, .25, 0.5, 1, 2, 4, 8))+
    scale_y_log10(breaks=c(.125, .25, 0.5, 1, 2, 4, 8), labels=c(.125, .25, 0.5, 1, 2, 4, 8))+
    scale_alpha_discrete(range = c(0.95, 0.96), guide=F)+
    facet_wrap(~Cov, ncol=3, drop=F)+
    ylab("odds ratio for mutability")+
    guides(fill = guide_legend(nrow = 3))+
    theme_coef()
ggsave(paste0(parentdir, "/images/coef_violin_cpgi.pdf"), width=6.5, height=6.5)

##############################################################################
# Plot CpG-specific effects
##############################################################################

cpg_covs <- c("LAMIN", "EXON", "H3K4me1", "H3K4me3", "H3K27ac")

cpg_plotdat <- sigcoefs %>%
  filter(Cov %in% cpg_covs) %>%
  filter(n>=10)

cpg_plotdat$Cov <- factor(cpg_plotdat$Cov, levels=cpg_covs)
sig10txt3 <- sig10txt %>%
  filter(Cov %in% cpg_covs)
sig10txt3$Cov <- factor(sig10txt3$Cov, levels=cpg_covs)

cpg_plotdat %>%
  filter(exp(Est) > 0.4) %>%
  ggplot(aes(x=Cov, y=exp(Est), alpha=factor(dir), fill=Category2))+
    geom_hline(yintercept=1, linetype="dashed")+
    geom_text(data=sig10txt3,
      aes(x=Cov, y=Estmax, label=label, colour=Category2),
      vjust=1, size=6, angle=90)+
    geom_violin(position="identity", scale="area")+
    scale_fill_manual(values=gp_cols, drop=FALSE)+
    scale_colour_manual(values=gp_cols, drop=FALSE)+
    scale_x_discrete(drop=FALSE)+
    # scale_y_log10(breaks=c(.125, .25, 0.5, 1, 2, 4, 8))+
    scale_y_log10(breaks=c(0.5, 1, 2, 4), limits=c(0.25,5))+
    scale_alpha_discrete(range = c(0.95, 0.96), guide=F)+
    # facet_wrap(~Cov, ncol=4, drop=F)+
    ylab("odds ratio for mutability")+
    guides(fill = FALSE, colour=FALSE)+
    theme_coef()+
    theme(axis.text.x=element_text(angle=90))
ggsave(paste0(parentdir, "/images/coef_violin_cpg_effects.pdf"), width=6.5, height=6.5)

##############################################################################
# Code below is used to validate specific results of feature-associated subtypes
# -Must already have read in DNM data
#
# runTest() function counts the total observed number of DNMs within
# feature-associated subtypes, the number of these within each feature, and the
# expected number assuming no effect of genomic features, then tests for
# difference between expected/observed
##############################################################################
getfeats <- function(sitedf){
  # Add histone marks to site data
  hists <- c("H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3",
    "H3K27ac", "H3K27me3", "H3K36me3")
  dflist <- list()
  for(i in 1:length(hists)){
    mark <- hists[i]
    file <- paste0(parentdir, "/reference_data/sort.E062-",
      mark, ".bed")
    hist <- binaryCol(sitedf, file)
    dflist[[i]] <- hist
  }

  df <- as.data.frame(do.call(cbind, dflist))
  names(df) <- hists
  sitedf <- cbind(sitedf, df)

  # Add other features
  sitedf$EXON <- binaryCol(sitedf,
    paste0(parentdir, "/reference_data/GRCh37_RefSeq_sorted.bed"))
  sitedf$CpGI <- binaryCol(sitedf,
    paste0(parentdir, "/reference_data/cpg_islands_sorted.bed"))
  sitedf$RR <- rcrCol(sitedf,
    paste0(parentdir, "/reference_data/recomb_rate.bed"))
  sitedf$LAMIN <- binaryCol(sitedf,
    paste0(parentdir, "/reference_data/lamin_B1_LADS2.bed"))
  sitedf$DHS <- binaryCol(sitedf,
    paste0(parentdir, "/reference_data/DHS.bed"))
  sitedf$TIME <- repCol(sitedf,
    paste0(parentdir, "/reference_data/lymph_rep_time.txt"))
  sitedf$GC <- gcCol(sitedf,
    paste0(parentdir, "/reference_data/gc10kb.bed"))

  return(sitedf)
}

runTest <- function(cov, dir){
  covtmp <- covdat %>% filter(Cov==cov)
  if(dir=="Up"){
    covdir <- covtmp %>% filter(Est>0)
    diralt <- "greater"
  } else {
    covdir <- covtmp %>% filter(Est<0)
    diralt <- "less"
  }

  covdirsub <- paste0(covdir$Category, "_", covdir$Sequence)

  # subset DNMs to those in associated subtypes
  dnmstmp <- input_dnms %>%
    mutate(Category2 = ifelse(substr(SEQ, 4, 5)=="CG",
      paste0("cpg_", Category), Category)) %>%
    mutate(subtype=paste0(Category, "_", substr(SEQ, 1, 7))) %>%
    # filter(subtype %in% paste0(covdir$Category, "_", covdir$Sequence))
    mutate(insub=ifelse(subtype %in% covdirsub, 1, 0))
    # filter(Category2 == "cpg_GC_AT")
    # filter(Category %in% covdir$Category &
    #   substr(SEQ, 1, 7) %in% covdir$Sequence)

  dnms2 <- dnmstmp %>%
    filter(insub==1) %>%
    mutate(DNM=1)
  dirlist <- sapply(unlist(covdir$Sequence),
    function(x)
      paste0(parentdir, "/output/logmod_data/motifs3/", x, ".txt"))

  motif_positions <- do.call(rbind,
    lapply(dirlist,
      function(x) read.table(x, header=F, stringsAsFactors=F)))

  mut_cats <- c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA")
  names(motif_positions) <- c("CHR", "POS", "Sequence", mut_cats, "DP")
  sites <- merge(motif_positions, dnms2[,c("CHR", "POS", "DNM")],
    by=c("CHR", "POS"), all.x=TRUE) %>%
    mutate(DNM=ifelse(DNM==1, 1, 0))

  sites2 <- sites %>%
    sample_n(100000)

  sites2 <- rbind(sites2, sites[sites$DNM==1,])



  sites3 <- getfeats(sites2)
  dnm_mod_formula <- as.formula(paste("DNM ~",
    paste(names(sites3)[-(1:12)], collapse="+")))
  dnm_mod <- glm(dnm_mod_formula, data=sites3, family="binomial")

  result <- tidy(dnm_mod) %>% filter(term==cov)
  row <- c(dir, result)
  return(row)
  ###
}

##############################################################################
# Test for mutation rate differences in dnms
##############################################################################
covdat <- sigcoefs
covdat$Category <- gsub("cpg_", "", covdat$Category)
covs <- unique(covdat$Cov)
# covs <- covs[grepl("H3", covs)]
testdat <- data.frame()
dnmdat <- data.frame()
covdir <- covdat %>%
  mutate(Dir=ifelse(Est>=0, "Up", "Down")) %>%
  group_by(Cov, Dir) %>%
  summarise(n=n())
for(i in 1:nrow(covdir)){
# for(i in 1:5){
  cov <- covdir[i,]$Cov
  dir <- covdir[i,]$Dir
  cat("Running tests on", cov, dir, "\n")
  if(covdir[i,]$n >= 10){
    tmpdata <- runTest(cov, dir)
    # testrow <- tmpdata$test
    testdat <- rbind(testdat, tmpdata)

    # tmpdnms <- tmpdata$dnms
    # dnmdat <- rbind(dnmdat, tmpdnms)
  }
}

dnmdat %>%
  mutate(cov=ifelse(dir=="Down" & cov=="TIME", "TIME_LATE", cov)) %>%
  mutate(cov=ifelse(dir=="Up" & cov=="TIME", "TIME_EARLY", cov)) %>%
  mutate(dir=ifelse(dir=="Down" & cov=="TIME_LATE", "Up", dir)) %>%
  filter(inside %in% c(0,1)) %>%
  ggplot(aes(x=cov, y=log(MU), colour=factor(inside)))+
    geom_violin()+
    geom_boxplot(fill=NA)+
    # geom_point(aes(colour=factor(inside)), alpha=0.1,
    #   position=position_jitterdodge())+
    facet_wrap(~dir, scales="free", nrow=2)+
    ylab("log(relative mutation rate)")
ggsave(paste0(parentdir, "/images/covbars_full.png"), width=12, height=8)

cov_ord <- c("H3K9me3", "RR", "H3K27me3",
  "H3K27ac", "EXON", "H3K4me1", "H3K4me3", "LAMIN",
  "GC", "H3K36me3", "CpGI", "TIME", "DHS")

testdat$cov <- factor(testdat$cov, levels=cov_ord)

testdat %>%
  filter(cov %in% cov_ord) %>%
  dplyr::select(cov, dir,
    Inside=estimate2, Outside=estimate1,
    statistic, df=parameter, p.value) %>%
    arrange(cov)

##############################################################################
# Test for mutation rate differences in dnms, excluding CpGs
##############################################################################
covdat <- sigcoefs %>%
  filter(!grepl("cpg", Category))
covdat$Category <- gsub("cpg_", "", covdat$Category)
covs <- unique(covdat$Cov)
# covs <- covs[grepl("H3", covs)]
testdat3 <- data.frame()
covdir <- covdat %>%
  mutate(Dir=ifelse(Est>=0, "Up", "Down")) %>%
  group_by(Cov, Dir) %>%
  summarise(n=n())
for(i in 1:nrow(covdir)){
  cov <- covdir[i,]$Cov
  dir <- covdir[i,]$Dir
  if(covdir[i,]$n > 10){
    row <- runTest(cov, dir)
    testdat3 <- rbind(testdat3, row)
  }
}

##############################################################################
# Test for mutation rate differences in DHS CCAAT-box motifs
##############################################################################
covdat <- sigcoefs %>%
  filter(Category == "AT_GC" & substr(Sequence, 1, 5)=="CCAAT")
covdat$Category <- gsub("cpg_", "", covdat$Category)
covs <- unique(covdat$Cov)

row <- runTest("DHS", "Up")

# covs <- covs[grepl("H3", covs)]
testdat2 <- data.frame()
covdir <- covdat %>%
  mutate(Dir=ifelse(Est>=0, "Up", "Down")) %>%
  group_by(Cov, Dir) %>%
  summarise(n=n())
for(i in 1:nrow(covdir)){
  cov <- covdir[i,]$Cov
  dir <- covdir[i,]$Dir
  if(covdir[i,]$n > 10){
    row <- runTest(cov, dir)
    testdat2 <- rbind(testdat2, row)
  }
}

###
testdat[testdat$cov=="TIME",]$dir <- "Up"
names(testdat)[3:4]<-c("Observed", "Expected")

td2 <- testdat %>%
  gather(group, value, c(Observed:Expected))
td3 <- data.frame()
for(i in 1:nrow(td2)){
  row <- td2[i,]
  ci <- prop.test(row$value, row$n)$conf.int[1:2]*row$n
  names(ci) <- c("lo", "hi")
  if(row$dir=="Up"){
    binom.pval <- binom.test(x=as.integer(row$value), n=row$n, p=row$propexp,
      alternative="greater")$p.value
  } else {
    binom.pval <- binom.test(x=as.integer(row$value), n=row$n, p=row$propexp,
      alternative="less")$p.value
  }

  newrow <- data.frame(c(row, ci, binom.pval=binom.pval), stringsAsFactors=F)
  td3 <- rbind(td3, newrow)
}

limits <- aes(ymax=hi, ymin=lo)

dirlabels <- c(Down="Subtypes with predicted depletion",
  Up="Subtypes with predicted enrichment")

td3$gp <- paste0(td3$dir, " (", td3$group, ")")
# td3$gp <- as.factor(paste0(td3$dir, " (", td3$group, ")"))
# levels(td3$group) <- c("Expected", "Observed")
# levels(td3$gp) <- unique(td3$gp)[c(3:4, 1:2)]

ggplot(td3[td3$value>60,],
    aes(x=cov, y=value, colour=gp, fill=gp, group=group))+
  geom_bar(stat="identity", position="dodge")+
  geom_errorbar(limits, position="dodge", colour="black")+
  # scale_y_log10(limits=c(60, 2000))+
  facet_wrap(~dir, scales="free", nrow=2, labeller=labeller(dir=dirlabels))+
  scale_colour_manual(values=rb[c(1:2,4:5)])+
  scale_fill_manual(values=rb[c(1:2,4:5)],
    labels=rep(c("Expected #DNMs under null\n(features have no effect)",
      "Observed #DNMs"), 2))+
  ylab("#DNMs")+
  theme_classic()+
  theme(legend.title=element_blank(),
    legend.key.size = unit(1.1, "cm"),
    strip.text.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.x = element_text(size=14, angle=45, vjust=1, hjust=1),
    axis.title.x=element_blank())+
  guides(colour=FALSE)
ggsave(paste0(parentdir, "/images/covbars_full.png"), width=12, height=8)
