# Read table of coefficients output from logit models
coefs <- read.table(paste0(parentdir, "/output/logmod_data/coefs/coefs_full.txt"), header=F, stringsAsFactors=F)
names(coefs) <- c("Cov", "Est", "SE", "Z", "pval", "Sequence", "Category")
coefs$Category <- ifelse(substr(coefs$Sequence, 4, 5)=="CG", paste0("cpg_", coefs$Category), coefs$Category)

coefsout <- coefs %>%
  mutate(Category = plyr::mapvalues(Category, orderedcats1, orderedcats2))
coefsout$Category <- factor(coefsout$Category, levels=orderedcats2)

write.table(coefsout, paste0(parentdir, "/output/logmod_data/coefs/supplementary_table_7.txt"),
  quote=F, col.names=T, row.names=F, sep="\t")
# coefs$Cov <- factor(plotcts$Cov, levels=orderedcovs)

ratefile <- paste0(parentdir, "/output/7bp_1000k_rates.txt")
rates <- read.table(ratefile, header=T, stringsAsFactors=F)

rates <- rates %>%
  mutate(Category=gsub("cpg_", "", Category2)) %>%
  mutate(Sequence=substr(Sequence, 0, 7))
motifdat <- rates %>%
  dplyr::select(Category, Sequence, rel_prop)

orderedcovs <- c("H3K9me3", "H3K27me3", "LAMIN", "BLANK",
  "RR", "TIME", "H3K4me1", "H3K4me3",
  "H3K36me3", "DHS", "GC", "CpGI")

sigcoefs <- coefs %>%
  group_by(Cov, Category) %>%
  # group_by(Category, Sequence) %>%
  mutate(qval=p.adjust(pval, method="fdr")) %>%
  filter(qval<0.05) %>%
  filter(!grepl("Intercept|DP", Cov)) %>%
  mutate(n=n()) %>%
  filter(n>=10) %>%
  ungroup()

# Update data for plotting
plotdat <- sigcoefs %>%
  mutate(Category = plyr::mapvalues(Category, orderedcats1, orderedcats2))
plotdat$Category <- factor(plotdat$Category, levels=orderedcats2)
plotdat$Cov <- factor(plotdat$Cov, levels=orderedcovs)
plotdat$dir <- ifelse(plotdat$Est>0, "Up", "Down")
plotdat$Est <- ifelse(plotdat$Cov=="RR", plotdat$Est*10, plotdat$Est)
plotdat$Est <- ifelse(plotdat$Cov=="GC", plotdat$Est/10, plotdat$Est)

plotcts <- plotdat %>%
  group_by(Cov, Category, dir) %>%
  summarise(n=n()) %>%
  mutate(Estmax=ifelse(dir=="Up", 3, 0.5)) %>%
  mutate(Estmax=ifelse(Cov=="CpGI", ifelse(dir=="Up", 10, 0.125), Estmax), label=n) %>%
  filter(n>4)
plotcts$Category <- factor(plotcts$Category, levels=orderedcats2)
plotcts$Cov <- factor(plotcts$Cov, levels=orderedcovs)

# ggplot()+
svglite(paste0(parentdir, "/images/coef_violin2.svg"), width=7, height=7)
plotdat %>%
  filter(Cov!="CpGI") %>%
  ggplot(aes(x=Category, y=exp(Est), alpha=factor(dir), fill=Category))+
    geom_hline(yintercept=1, linetype="dashed")+
    geom_text(data=plotcts[plotcts$Cov!="CpGI",],
      aes(x=Category, y=Estmax, label=label, colour=Category), vjust=1, size=4, angle=90)+
    geom_violin(position="identity", scale="area")+
    scale_fill_manual(values=cols, drop=FALSE)+
    scale_colour_manual(values=cols, drop=FALSE)+
    scale_x_discrete(drop=FALSE)+
    # scale_y_log10(breaks=c(.125, .25, 0.5, 1, 2, 4, 8))+
    scale_y_log10(breaks=c(0.5, 1, 2, 4))+
    scale_alpha_discrete(range = c(0.95, 0.96), guide=F)+
    facet_wrap(~Cov, ncol=4, drop=F)+
    ylab("odds ratio for mutability")+
    theme_bw()+
    guides(fill = guide_legend(nrow = 3))+
    theme(legend.position="bottom",
    legend.title=element_blank(),
    strip.text=element_text(size=12),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.text.y=element_text(size=10),
    axis.title.y=element_text(size=12))
dev.off()
# ggsave(paste0(parentdir, "/images/coef_violin2.svg"), width=7, height=7)

# my_svg <- function(file, width, height) {
#    library(RSvgDevice)
#    devSVG(file = file, width = width, height = height, bg = "white", fg = "black",
#           onefile = TRUE, xmlHeader = TRUE)
# }

svglite(paste0(parentdir, "/images/coef_violin2_cpgi2.svg"), width=7, height=7)
plotdat %>%
  ggplot(aes(x=Category, y=exp(Est), alpha=factor(dir), fill=Category))+
    geom_hline(yintercept=1, linetype="dashed")+
    geom_text(data=plotcts[plotcts$Cov=="CpGI",],
      aes(x=Category, y=Estmax, label=label, colour=Category), vjust=1, angle=90, size=4)+
    geom_violin(position="identity", scale="area")+
    scale_fill_manual(values=cols, drop=FALSE)+
    scale_colour_manual(values=cols, drop=FALSE)+
    scale_x_discrete(drop=FALSE)+
    scale_y_log10(breaks=c(.125, .25, 0.5, 1, 2, 4, 8), labels=c(.125, .25, 0.5, 1, 2, 4, 8))+
    # scale_y_log10(breaks=c(0.5, 1, 2))+
    scale_alpha_discrete(range = c(0.95, 0.96), guide=F)+
    facet_wrap(~Cov, ncol=4, drop=F)+
    # facet_wrap(~Cov, ncol=4, drop=T)+
    ylab("odds ratio for mutability")+
    theme_bw()+
    guides(fill = guide_legend(nrow = 3))+
    theme(legend.position="bottom",
    legend.title=element_blank(),
    strip.text=element_text(size=12),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.text.y=element_text(size=10),
    axis.title.y=element_text(size=12))
dev.off()
# ggsave(paste0(parentdir, "/images/coef_violin2_cpgi2.svg"), width=7, height=7)

coefs <- coefs %>%
  mutate(Category = plyr::mapvalues(Category, orderedcats1, orderedcats2))

coefs$Category <- factor(coefs$Category, levels=orderedcats2)
coefs %>%
  filter(!grepl("Intercept", Cov)) %>%
  filter(exp(Est)<10) %>%
  ggplot(aes(fill=Category))+
    geom_hline(yintercept=1, linetype="dashed")+
    geom_violin(aes(x=Category, y=exp(Est), fill=Category))+
    scale_fill_manual(values=cols, drop=FALSE)+
    scale_colour_manual(values=cols, drop=FALSE)+
    facet_wrap(~Cov, ncol=4, scales="free_y")+
    ylab("odds ratio for mutability")+
    theme_bw()+
    guides(fill = guide_legend(nrow = 3))+
    theme(legend.position="bottom",
    strip.text=element_text(size=16),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.text.y=element_text(size=14),
    axis.title.y=element_text(size=16))
ggsave(paste0(parentdir, "/images/coef_violin_full.png"), width=9, height=8)

##############################################################################
# Code below is used to validate specific results of feature-associated subtypes
# -Must already have read in DNM data
#
# runTest() function counts the total observed number of de novo mutations within
# feature-associated subtypes, the number of these within each feature, and the
# expected number assuming no effect of genomic features, then tests for
# difference between expected/observed
##############################################################################
runTest <- function(cov, dir){
  covtmp <- covdat %>% filter(Cov==cov)
  if(dir=="Up"){
    covdir <- covtmp %>% filter(Est>0)
  } else {
    covdir <- covtmp %>% filter(Est<0)
  }

  dnmstmp <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", substr(chrpfdnm$SEQ, 0, 7)) %in%
    paste0(covdir$Category, "_", covdir$Sequence),]
  if(nrow(dnmstmp)>=5){
    if(cov=="GC"){
      covbase <- paste0(parentdir, "/reference_data/high_gc")
      dnmstmp$GC <- gcCol(dnmstmp,
        paste0(parentdir, "/output/3bp_10k/full_bin.txt"))
      dnmstmp$inside <- ifelse(dnmstmp$GC>=0.55, 1, 0)
    } else if(cov=="TIME"){
      dnmstmp$TIME <- repCol(dnmstmp,
        paste0(parentdir, "/reference_data/lymph_rep_time.txt"))
      if(dir=="Down"){
        covbase <- paste0(parentdir, "/reference_data/late_rt")
        dnmstmp$inside <- ifelse(dnmstmp$TIME<=-1.25, 1, 0)
      } else if(dir=="Up"){
        covbase <- paste0(parentdir, "/reference_data/early_rt")
        dnmstmp$inside <- ifelse(dnmstmp$TIME>=1.25, 1, 0)
      }
    } else if(cov=="RR"){
      covbase <- paste0(parentdir, "/reference_data/high_rr")
      dnmstmp$RR <- rcrCol(dnmstmp,
        paste0(parentdir, "/reference_data/recomb_rate.bed"))
      dnmstmp$inside <- ifelse(dnmstmp$RR>=2, 1, 0)
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
    } else {
      covbase <- paste0(parentdir, "/reference_data/histone_marks/broad/sort.E062-", cov)
      covbed <- paste0(covbase, ".bed")
      dnmstmp$inside <- binaryCol(dnmstmp, covbed)
    }

    obs <- sum(dnmstmp$inside)

    seqs <- unlist(c(covdir$Sequence, lapply(covdir$Sequence, revcomp)))
    write.table(seqs, paste0(parentdir, "/seqs.txt"), col.names=F, row.names=F, quote=F, sep="\t")

    grepcmd <- paste0("grep -o -Ff ", parentdir, "/seqs.txt ", covbase, ".fa | sort | uniq -c > ", parentdir, "/testcounts.txt")
    system(grepcmd)

    motifcts <- read.table(paste0(parentdir, "/testcounts.txt"), header=F, stringsAsFactors=F)

    names(motifcts) <- c("Count", "SEQ")
    motifcts$REVSEQ <- unlist(lapply(motifcts$SEQ, revcomp))
    motifcts$Sequence <- ifelse(substr(motifcts$SEQ,4,4) %in% c("A", "C"),
      motifcts$SEQ, motifcts$REVSEQ)
      # paste0(motifcts$SEQ, "(", motifcts$REVSEQ, ")"),
      # paste0(motifcts$REVSEQ, "(", motifcts$SEQ, ")"))
    motifcts2 <- motifcts %>%
      group_by(Sequence) %>%
      summarise(Count=sum(Count))
    covdir2 <- merge(covdir, motifcts2, by=c("Sequence"))

    covdir3 <- merge(covdir2, motifdat, by=c("Category", "Sequence"))
    covdir3$exp <- covdir3$Count*covdir3$rel_prop*1.67e-06*1074
    exp <- sum(covdir3$exp)

    total <- nrow(dnmstmp)
    if(dir=="Up"){
      test <- prop.test(c(obs, exp), c(total, total), alternative="greater")
    } else {
      test <- prop.test(c(obs, exp), c(total, total), alternative="less")
    }
    newrow <- data.frame(cov, dir=dir, obs=obs, exp=exp, n=total,
      propobs=test$estimate[1], propexp=test$estimate[2], pval=test$p.value)
    # testdat <- rbind(testdat, newrow)
    newrow
  }
}

covdat <- sigcoefs
covdat$Category <- gsub("cpg_", "", covdat$Category)
covs <- unique(covdat$Cov)
# covs <- covs[grepl("H3", covs)]
testdat <- data.frame()
covdir <- covdat %>%
  mutate(Dir=ifelse(Est>=0, "Up", "Down")) %>%
  group_by(Cov, Dir) %>%
  summarise(n=n())
for(i in 1:nrow(covdir)){
  cov <- covdir[i,]$Cov
  dir <- covdir[i,]$Dir
  if(covdir[i,]$n > 10){
    row <- runTest(cov, dir)
    testdat <- rbind(testdat, row)
  }
}

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
