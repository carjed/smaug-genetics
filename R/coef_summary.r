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
  filter(n>=10)
non_cpg_plotdat$Cov <- factor(non_cpg_plotdat$Cov, levels=non_cpg_covs)
sig10txt1 <- sig10txt %>%
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

##############################################################################
# Plot CpG island effects (w/ >10 subtypes)
##############################################################################
cpgi_covs <- c("CpGI", paste0("BLANK", 1:5))

cpgi_plotdat <- sigcoefs %>%
  filter(Cov %in% cpgi_covs) %>%
  filter(n>=10)
cpgi_plotdat$Cov <- factor(cpgi_plotdat$Cov, levels=cpgi_covs)
sig10txt2 <- sig10txt %>%
  filter(Cov %in% cpgi_covs)
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
    scale_y_log10(breaks=c(.125, .25, 0.5, 1, 2, 4, 8))+
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
    scale_y_log10(breaks=c(0.5, 1, 2, 4))+
    scale_alpha_discrete(range = c(0.95, 0.96), guide=F)+
    # facet_wrap(~Cov, ncol=4, drop=F)+
    ylab("odds ratio for mutability")+
    guides(fill = FALSE, colour=FALSE)+
    theme_coef()+
    theme(axis.text.x=element_text(angle=90))
ggsave(paste0(parentdir, "/images/coef_violin_cpg_effects.pdf"), width=6.5, height=3)

##############################################################################
# Code below is used to validate specific results of feature-associated subtypes
# -Must already have read in DNM data
#
# runTest() function counts the total observed number of DNMs within
# feature-associated subtypes, the number of these within each feature, and the
# expected number assuming no effect of genomic features, then tests for
# difference between expected/observed
##############################################################################
runTest <- function(cov, dir){
  covtmp <- covdat %>% filter(Cov==cov)
  if(dir=="Up"){
    covdir <- covtmp %>% filter(Est>0)
    diralt <- "greater"
  } else {
    covdir <- covtmp %>% filter(Est<0)
    diralt <- "less"
  }

  # subset DNMs to those in associated subtypes
  dnmstmp <- input_dnms %>%
    mutate(Category2 = ifelse(substr(SEQ, 4, 5)=="CG",
      paste0("cpg_", Category), Category)) %>%
    mutate(subtype=paste0(Category, "_", substr(SEQ, 1, 7))) %>%
    filter(subtype %in% paste0(covdir$Category, "_", covdir$Sequence)) %>%
    filter(Category2 == "cpg_GC_AT")
    # filter(Category %in% covdir$Category &
    #   substr(SEQ, 1, 7) %in% covdir$Sequence)

  if(nrow(dnmstmp)>=10){
    if(cov=="GC"){
      covbase <- paste0(parentdir, "/reference_data/high_gc")
      covbed <- paste0(parentdir, "/reference_data/gc10kb.bed")
      dnmstmp$GC <- gcCol(dnmstmp, covbed)
      dnmstmp$inside <- ifelse(dnmstmp$GC>=0.55, 1, 0)
      awkstr <- "awk -F \"\\t\" '{if($4>=0.55 && NR>1) print }'"
    } else if(cov=="TIME"){
      covbed <- paste0(parentdir, "/reference_data/lymph_rep_time.txt")
      dnmstmp$TIME <- repCol(dnmstmp, covbed)
      if(dir=="Down"){
        covbase <- paste0(parentdir, "/reference_data/late_rt")
        dnmstmp$inside <- ifelse(dnmstmp$TIME<=-1.25, 1, 0)
        awkstr <- "awk -F \"\\t\" '{if($3<=-1.25 && NR>1) print $1 \"\\t\"$2-500+1\"\\t\"$2+500}'"
      } else if(dir=="Up"){
        covbase <- paste0(parentdir, "/reference_data/early_rt")
        dnmstmp$inside <- ifelse(dnmstmp$TIME>1.25, 1, 0)
        awkstr <- "awk -F \"\\t\" '{if($3>1.25 && NR>1) print $1 \"\\t\"$2-500+1\"\\t\"$2+500}'"
      }
    } else if(cov=="RR"){
      covbase <- paste0(parentdir, "/reference_data/high_rr")
      covbed <- paste0(parentdir, "/reference_data/recomb_rate.bed")
      dnmstmp$RR <- rcrCol(dnmstmp, covbed)
      dnmstmp$RR[is.na(dnmstmp$RR)] <- 0
      dnmstmp$inside <- ifelse(dnmstmp$RR>=2, 1, 0)
      awkstr <- "awk -F\"\\t\" '{if($4>=2 && NR>1) print $1\"\\t\"$2\"\\t\"$3}'"
    } else if(cov=="DHS"){
      covbase <- paste0(parentdir, "/reference_data/DHS")
      covbed <- paste0(covbase, ".bed")
      dnmstmp$inside <- binaryCol(dnmstmp, covbed)
    } else if(cov=="CpGI"){
      covbase <- paste0(parentdir, "/reference_data/cpg_islands_sorted")
      covbed <- paste0(covbase, ".bed")
      dnmstmp$inside <- binaryCol(dnmstmp, covbed)
    } else if(cov=="LAMIN"){
      covbase <- paste0(parentdir, "/reference_data/lamin_B1_LADS2")
      covbed <- paste0(covbase, ".bed")
      dnmstmp$inside <- binaryCol(dnmstmp, covbed)
    } else if(cov=="EXON"){
      covbase <- paste0(parentdir, "/reference_data/GRCh37_RefSeq_sorted")
      covbed <- paste0(covbase, ".bed")
      dnmstmp$inside <- binaryCol(dnmstmp, covbed)
    } else {
      covbase <- paste0(parentdir,
        "/reference_data/sort.E062-", cov)
      covbed <- paste0(covbase, ".bed")
      dnmstmp$inside <- binaryCol(dnmstmp, covbed)
    }

    # seqs <- unique(unlist(c(covdir$Sequence, lapply(covdir$Sequence, revcomp))))
    # write out list of sequence motifs (+ reverse complements)
    # that contain a DNM
    din <- dnmstmp #%>% filter(inside==1)
    seqs <- unique(unlist(c(substr(din$SEQ, 1, 7), lapply(substr(din$SEQ, 1, 7), revcomp))))
    write.table(seqs, paste0(parentdir, "/seqs.txt"),
      col.names=F, row.names=F, quote=F, sep="\t")

    # Create fasta file from feature file
    if(!file.exists(paste0(covbase, ".fa"))){

      if(cov=="RR"){
        getfastacmd <- paste0("bedtools getfasta -fi ", parentdir,
          "/reference_data/human_g1k_v37.fasta -bed <(cut -f1-4 ", covbed,
          " | sed 's/chr//g' | ", awkstr,
          " | sed /^X/d | sed /^Y/d) -fo ", covbase, ".fa")

      } else if (cov=="TIME"){
        getfastacmd <- paste0("bedtools getfasta -fi ", parentdir,
          "/reference_data/human_g1k_v37.fasta -bed <(cut -f1-3 ", covbed,
          " | sed 's/chr//g' | ", awkstr,
          " | sed /^X/d | sed /^Y/d) -fo ", covbase, ".fa")
      } else if (cov=="GC"){
        getfastacmd <- paste0("bedtools getfasta -fi ", parentdir,
          "/reference_data/human_g1k_v37.fasta -bed <(cut -f1-5 ", covbed,
          " | sed 's/chr//g' | ", awkstr,
          " | sed /^X/d | sed /^Y/d) -fo ", covbase, ".fa")
      } else {
        getfastacmd <- paste0("bedtools getfasta -fi ", parentdir,
          "/reference_data/human_g1k_v37.fasta -bed <(sed 's/chr//g' ", covbed,
          " | sed /^X/d | sed /^Y/d) -fo ", covbase, ".fa")
      }

      get_fasta_file <- paste0(parentdir, "/", cov, "_", dir, "get_fasta.sh")
      fileConn<-file(get_fasta_file)
      writeLines(c("#!/bin/bash", getfastacmd), fileConn)
      close(fileConn)
      # system(getfastacmd)
      system(paste0("bash ", get_fasta_file), ignore.stderr=TRUE)
    }

    # grepcmd <- paste0("grep -o -Ff ",
    #   parentdir, "/seqs.txt ", covbase, ".fa | sort | uniq -c > ",
    #   parentdir, "/testcounts.txt")
    # system(grepcmd)

    if(!file.exists(paste0(covbase, "_testcounts.txt"))){
      count_motifs_cmd <- paste0("python data_mgmt/lib/motif_count.py",
        " -i ", covbase, ".fa",
        " -m ", parentdir, "/seqs.txt",
        " -o ", covbase, "_testcounts.txt")

      system(count_motifs_cmd)
    }

    motifcts <- read.table(paste0(covbase, "_testcounts.txt"),
      header=F, stringsAsFactors=F)

    names(motifcts) <- c("SEQ", "Count")
    motifcts$REVSEQ <- unlist(lapply(motifcts$SEQ, revcomp))
    motifcts$Sequence <- ifelse(substr(motifcts$SEQ,4,4) %in% c("A", "C"),
      motifcts$SEQ, motifcts$REVSEQ)

    motifcts2 <- motifcts %>%
      group_by(Sequence) %>%
      summarise(Count=sum(Count))
    covdir2 <- merge(motifcts2, covdir[,c("Sequence", "nMotifs")], by=c("Sequence")) %>%
      group_by(Sequence, Count) #%>%
      # filter(row_number(nMotifs) == 1)
      # filter(paste0(Category, "_", Sequence) %in% dnmstmp$subtype)

    ndnms_in <- sum(dnmstmp$inside, na.rm=T)
    ndnms_out <- nrow(dnmstmp)-ndnms_in

    nmotifs_in <- sum(covdir2$Count)
    nmotifs_out <- sum(covdir2$nMotifs) - nmotifs_in

    # outdat <- list()
    if(ndnms_in > 5){
      test_result <- tidy(
        prop.test(c(ndnms_in, ndnms_out),
                  c(nmotifs_in, nmotifs_out), alternative=diralt))

      # test_result <- dnmstmp %>% do(tidy(t.test(MU~inside, data=., alternative=diralt)))
      nseqs <- length(seqs)
      outdat <- cbind(cov, dir, test_result)
      # outdat$dnms <- dnmstmp %>%
      #   dplyr::select(CHR, POS, MU, inside) %>%
      #   mutate(cov=cov, dir=dir)
      return(outdat)
     } # else {
    #   return(rep(0,12))
    # }

  }
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
