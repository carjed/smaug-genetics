###############################################################################
# Analysis by sample
###############################################################################
i <- 1
cbp <- 5

bindir3 <- paste0(parentdir, "/motif_counts/3-mers/full")
p2 <- "1000kb_full.txt"
bins1Mb <- get_bins(bindir3, p2)
bins1Mb$CHR <- gsub("chr", "", bins1Mb$CHR)

get_mct_b <- function(bins){
  out <- bins %>%
  	group_by(CHR, BIN, Motif) %>%
  	summarise(nMotifs=sum(nMotifs))
  return(out)
}

# Count singletons per 3-mer subtype per individual
ind_counts <- full_data$sites %>%
  mutate(Type=gsub("cpg_", "", Category2),
    SEQA=substr(Motif, cbp-i, cbp+i),
    SEQB=substr(Motif, cbp*3-i, cbp*3+i),
    Motif=paste0(SEQA, "(", SEQB, ")")) %>%
  group_by(ID, Type, Motif) %>%
  summarise(n=n())

# count 3-mer motifs genome-wide
mct3 <- get_mct_b(bins1Mb) %>%
  group_by(Motif) %>%
  summarise(nMotifs=sum(nMotifs))

# get ids in ped file
ped <- read.table(pedfile, header=F, stringsAsFactors=F)
ped2 <- read.table(qcped, header=T, stringsAsFactors=F)

# get per-person rates and filter IDs
ind_counts <- merge(ind_counts, mct3, by="Motif") %>%
  filter(ID %in% ped$V1) %>%
  mutate(ERV_rel_rate=n/nMotifs,
    subtype=paste0(Type, "_", Motif))

# Coerce from long to wide format
# ind_wide <- ind_counts  %>%
#   dplyr::select(ID, subtype, n) %>%
#   spread(subtype, n)

# # Alternative to ERV rates--use
# iw2 <- sweep(ind_wide[,-1], 2, colSums(ind_wide[,-1]),`/`)
# ind_wide <- data.frame(ID=ind_wide[,1], iw2)

# # Run PCA on sample matrix
# ind_pca <- prcomp(scale(ind_wide[,-c(1)]))
#
# # Cumulative variance explained by 96 PCs
# ind_cumvar <- cumsum((ind_pca$sdev)^2) / sum(ind_pca$sdev^2)
#
# # Data frame containing sample IDs and components
# ind_wide_c_pcs <- cbind(ID=ind_wide$ID, as.data.frame(ind_pca$x))
#
# # Run NMF on sample matrix; get top signature for each sample
# # set.seed(36087318)
#
# # ner<-nmfEstimateRank(as.matrix(ind_wide[,-c(1)]), seq(2,5), nrun=10)

###############################################################################
# get signature loadings
###############################################################################
get_loads <- function(widedat, nmfdat){
  sigloads <- data.frame(subtype=names(widedat),
    sig1=coef(nmfdat)[1,],
    sig2=coef(nmfdat)[2,],
    sig3=coef(nmfdat)[3,]) %>%
    mutate(sig1=sig1/sum(sig1), sig2=sig2/sum(sig2), sig3=sig3/sum(sig3)) %>%
    gather(subtype, value)

  names(sigloads) <- c("subtype", "sig", "value")
  sigloads <- sigloads %>%
    mutate(Category=substr(subtype, 1, 5),
      Sequence=substr(subtype, 7, 14))

  return(sigloads)
}

###############################################################################
# plot signature loadings
###############################################################################
plot_loads <- function(sigloads){
  p <- ggplot(sigloads, aes(x=Sequence, y=value))+
    geom_bar(stat="identity")+
    facet_grid(sig~Category, scales="free_x")+
    # geom_label_repel(data=sigloads[sigloads$sig3>0.005,], aes(x=sig2, y=sig3, label=subtype))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90, hjust=1))
  return(p)
}

###############################################################################
# plot signature contribution across individuals
###############################################################################
plot_ind_sigs <- function(sigdat){
  p <- ggplot(sigdat,
      aes(x=ID, y=prob, group=cluster, colour=cluster)) +
    geom_line()+
    scale_x_discrete(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0), limits=c(0,1))+
    facet_wrap(~Study, scales="free_x", nrow=1)+
    ylab("signature contribution")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
      legend.position="bottom")

  return(p)
}

###############################################################################
# Prepare data and run NMF
###############################################################################
ind_wide <- ind_counts  %>%
  dplyr::select(ID, subtype, ERV_rel_rate) %>%
  spread(subtype, ERV_rel_rate)

ind_wide[is.na(ind_wide)] <- 0

# ind_wide <- ind_wide[,-c(64,60,56,52)]

ind_nmf <- nmf(as.matrix(ind_wide[,-c(1)]), 3, nrun=10)
ind_pred <- predict(ind_nmf, "rows", prob=T)

sigloads <- get_loads(ind_wide[,-c(1)], ind_nmf)
plot_loads(sigloads)
ggsave(paste0(parentdir, "/images/sigloads.png"))

# Add NMF basis and predicted class to data frame
ind_nmf_plotdat <- data.frame(ID=ind_wide$ID, basis(ind_nmf), sig=ind_pred[[1]])

# head(ind_nmf_plotdat)
# summary(ind_nmf_plotdat$sig)

samples <- ped %>%
  dplyr::select(ID=V1, Study=V7, BP=V6, PC1=V14, PC2=V15, PC3=V16, PC4=V17)

ind_nmf_plotdat <- merge(ind_nmf_plotdat, samples, by="ID") %>%
  mutate(sum=X1+X2+X3,
    X1a=X1/sum,
    X2a=X2/sum,
    X3a=X3/sum)

ind_nmf_plotdat$ID <- factor(ind_nmf_plotdat$ID, levels = unique(ind_nmf_plotdat$ID))

ind_nmf_long <- ind_nmf_plotdat %>%
  gather(cluster, prob, X1a:X3a) %>%
  arrange(Study, sum)

plot_ind_sigs(ind_nmf_long)
ggsave(paste0(parentdir, "/images/by_ind_all.png"), width=12, height=8)

# cbind(ind_wide, g2=ind_nmf_plotdat$sig) %>%
#   mutate(g2=ifelse(g2==1, TRUE, FALSE)) %>%
#   mutate(key=paste0(ID, "+", g2)) %>%
#   dplyr::select(-c(ID, g2)) %>%
#   gather(key, val) %>% setNames(c("IDf", "key2", "val")) %>%
#   tidyr::separate(IDf, into=c("ID", "g2"), sep="[+]") %>% #head
#   mutate(Category=substr(key2, 1, 5),
#     Sequence=substr(key2, 7, 14)) %>%
#   group_by(Category, Sequence) %>%
#   do(tidy(t.test(val~g2, data=.))) %>% #head %>% data.frame
#   dplyr::select(Category, Sequence, estimate1, estimate2, p.value) %>%
#   # filter(Category=="GC_TA")
#   filter(p.value<0.05/96)

###############################################################################
# IDs to keep
###############################################################################
keep_cis <- ind_nmf_long %>%
  dplyr::select(cluster, prob) %>%
  group_by(cluster) %>%
  summarise_each(funs(mean,sd)) %>%
  mutate(conf.low=mean-1.96*sd, conf.high=mean+1.96*sd)
  # do(tidy(t.test(.$prob))) %>%
  # dplyr::select(cluster, conf.low, conf.high)

keep_ids <- ind_nmf_plotdat %>%
  filter(X1a > keep_cis[1,]$conf.low & X1a < keep_cis[1,]$conf.high) %>%
  filter(X2a > keep_cis[2,]$conf.low & X2a < keep_cis[2,]$conf.high) %>%
  filter(X3a > keep_cis[3,]$conf.low & X3a < keep_cis[3,]$conf.high) %>%
  dplyr::select(ID)

drop_ids <- ind_nmf_plotdat %>%
  filter(!(ID %in% keep_ids$ID)) %>%
  dplyr::select(ID)

ped2$drop <- (ped2$HAID %in% drop_ids$ID)
ped2 %>% do(tidy(t.test(gcbias~drop, data=.)))
ped2 %>%
  group_by(drop) %>%
  summarise(cov=mean(coverage),
    sing=mean(Singletons),
    gcbias=mean(gcbias),
    insmedian=mean(insmedian))

# plates <- read.table(plates, header=F, stringsAsFactors=F)
# names(plates) <- c("ID", "PLATE")
# ind_nmf_plotdat <- merge(ind_nmf_plotdat, plates, by="ID")
#
# ped2 %>%
#   filter(HAID %in% sig1$ID) %>%
#   dplyr::select(HAID, N_SNPs, Singletons, Study)

###############################################################################
# Repeat with kept IDs
###############################################################################
ind_wide2 <- ind_wide %>%
  filter(ID %in% keep_ids$ID)
ind_nmf2 <- nmf(ind_wide2[,-c(1)], 3, nrun=10)
ind_pred2 <- predict(ind_nmf2, "rows", prob=T)

sigloads <- get_loads(ind_wide2[,-c(1)], ind_nmf2)
plot_loads(sigloads)
ggsave(paste0(parentdir, "/images/sigloads2.png"))

ind_nmf_plotdat <- data.frame(ID=ind_wide2$ID, basis(ind_nmf2), sig=ind_pred2[[1]])

ind_nmf_plotdat <- merge(ind_nmf_plotdat, samples, by="ID") %>%
  mutate(sum=X1+X2+X3,
    X1a=X1/sum,
    X2a=X2/sum,
    X3a=X3/sum)

ind_nmf_plotdat$ID <- factor(ind_nmf_plotdat$ID, levels = unique(ind_nmf_plotdat$ID))

ind_nmf_long <- ind_nmf_plotdat %>%
  gather(cluster, prob, X1a:X3a) %>%
  arrange(Study, sum)

plot_ind_sigs(ind_nmf_long)
ggsave(paste0(parentdir, "/images/by_ind_keep.png"), width=12, height=8)

# test for significant differences between groups
cbind(ind_wide2, g2=ind_nmf_plotdat$sig) %>%
  mutate(g2=ifelse(g2==2, TRUE, FALSE)) %>%
  # mutate(g2=(ID %in% s2ids$ID)) %>% # dplyr::select(ID, AT_CG_AAA.TTT., g2) %>% group_by(g2) %>% summarise(val=mean(AT_CG_AAA.TTT.))
  # group_by(g2) %>%
  mutate(key=paste0(ID, "+", g2)) %>%
  dplyr::select(-c(ID, g2)) %>%
  gather(key, val) %>% setNames(c("IDf", "key2", "val")) %>%
  tidyr::separate(IDf, into=c("ID", "g2"), sep="[+]") %>% #head
  # tidyr::separate(key2, into=c("ID", "g2"), sep="[+]") %>%
  mutate(Category=substr(key2, 1, 5),
    Sequence=substr(key2, 7, 14)) %>%
  group_by(Category, Sequence) %>%
  do(tidy(t.test(val~g2, data=.))) %>% #head %>% data.frame
  dplyr::select(Category, Sequence, estimate1, estimate2, p.value) %>%
  # filter(Category=="GC_TA")
  filter(p.value<0.05/96)

ggplot()+
  geom_point(data=ind_nmf_plotdat[ind_nmf_plotdat$sig==3,],
    aes(x=PC4, y=PC2), colour="gray30", alpha=0.3)+
  geom_point(data=ind_nmf_plotdat[ind_nmf_plotdat$sig==1,],
    aes(x=PC4, y=PC2), colour="blue", alpha=0.3)+
  geom_point(data=ind_nmf_plotdat[ind_nmf_plotdat$sig==2,],
    aes(x=PC4, y=PC2), colour="yellow", alpha=0.3)+
  theme_bw()
ggsave("/net/bipolar/jedidiah/mutation/images/sig_snp_pcs.png", width=8, height=8)

###############################################################################
# Correlate with Alexandrov's somatic signatures
###############################################################################

somatic_corr <- 0

if(somatic_corr){
  cs <- read.table(paste0(parentdir, "/cancer_sigs.txt"), header=T, sep="\t")

  cs <- cs[,1:33]

  snames <- paste0("S", 1:30)

  names(cs) <- c("Category", "Seq3", "Type", snames)

  cs$Category <- as.factor(cs$Category)
  levels(cs$Category) <- c("GC_TA", "GC_CG", "GC_AT", "AT_TA", "AT_GC", "AT_CG")

  cs$Seq3a <- as.character(reverse(complement(DNAStringSet(cs$Seq3))))
  #cs<-cs %>% mutate(Seq3b=revcomp(Seq3))
  #cs$Seq3a<-apply(cs$Seq3, 2, function(x) revcomp(x))

  cs$Sequence <- ifelse(
    substr(cs$Seq3,2,2)<substr(cs$Seq3a,2,2),
    paste0(cs$Seq3,"(",cs$Seq3a,")"),
    paste0(cs$Seq3a,"(",cs$Seq3,")")
  )

  cs2 <- cs %>%
    arrange(Category, Sequence) %>%
    mutate(subtype=paste0(Category, "_", Sequence))

  bin_coefs <- data.frame(subtype=colnames(coef(bins_nmf)), t(coef(bins_nmf)))
  names(bin_coefs) <- c("subtype", "binsig1", "binsig2", "binsig3")

  ind_coefs <- data.frame(subtype=colnames(coef(ind_nmf)), t(coef(ind_nmf)))
  names(ind_coefs) <- c("subtype", "indsig1", "indsig2", "indsig3")

  cs2 <- merge(cs2, bin_coefs, by="subtype")
  cs2 <- merge(cs2, ind_coefs, by="subtype")

  corrheat <- cs2 %>%
    dplyr::select(c(S1:S30, binsig1:binsig3, indsig1:indsig3)) %>%
    correlate() %>%
    focus(c(binsig1:binsig3, indsig1:indsig3)) %>%
    gather(key, correlation, -rowname)

  corrheat$rowname <- as.factor(corrheat$rowname)
  levels(corrheat$rowname) <- paste0("S", 1:30)

  ggplot(corrheat, aes(x=key, y=rowname, fill=correlation))+
    geom_tile()+
    xlab("ERV Signature")+
    ylab("Somatic Signature")+
    scale_fill_gradient(low="white",high="darkgreen")
  ggsave(paste0(parentdir, "/images/sigheat.png"))

}

###############################################################################
# Additional QC
###############################################################################
ind_extra <- 0
if(ind_extra){
  sig1 <- ind_wide_c %>%
    filter(sig==1)

  sig1 <- ind_nmf_plotdat %>% filter(sig==1)
  sig2 <- ind_nmf_plotdat %>% filter(sig==2)

  ped2 %>%
    filter(HAID %in% sig1$ID)

  ped2 %>%
    dplyr::select(HAID, N_SNPs, Singletons, Study, coverage, qcfail) %>%
    filter(HAID %in% sig2$ID)

  gcb_hi<-ped2 %>%
    arrange(desc(gcbias)) %>%
    dplyr::select(ID=HAID, gcbias, Study) %>%
    head(300) %>%
    arrange(ID)

  sig1_2<-ind_nmf_plotdat %>%
    filter(sig==1) %>%
    arrange(desc(prob)) %>%
    dplyr::select(ID, Study, sig, PLATE) %>%
    arrange(ID) %>%
    distinct(ID, .keep_all = TRUE)

  sig1_2 %>% group_by(PLATE) %>% summarise(n=n())

  sig1_2$gcb_hi <- sig1_2$ID %in% gcb_hi$ID
  sig1_2 <- merge(sig1_2, gcb_hi, by="ID", all.x=T)
  sig1_2 %>% arrange(desc(gcbias)) %>% head(50)

}