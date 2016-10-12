suppressMessages(usePackage(ggplot2))
suppressMessages(usePackage(dplyr))
suppressMessages(usePackage(tidyr))
suppressMessages(usePackage(reshape2))
suppressMessages(usePackage(RColorBrewer))
suppressMessages(usePackage(MASS))
suppressMessages(usePackage(speedglm))
suppressMessages(usePackage(boot))
suppressMessages(usePackage(devtools))
suppressMessages(usePackage(psych))

source("./R/get_functions.r")

# Function counts the total observed number of de novo mutations within
# feature-associated subtypes, the number of these within each feature, and the
# expected number assuming no effect of genomic features, then tests for
# difference between expected/observed
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
      covbase <- "/net/bipolar/jedidiah/mutation/reference_data/high_gc"
      dnmstmp$GC <- gcCol(dnmstmp,
        "/net/bipolar/jedidiah/mutation/output/3bp_10k/full_bin.txt")
      dnmstmp$inside <- ifelse(dnmstmp$GC>=0.55, 1, 0)
    } else if(cov=="TIME"){
      dnmstmp$TIME <- repCol(dnmstmp,
        "/net/bipolar/jedidiah/mutation/reference_data/lymph_rep_time.txt")
      if(dir=="Down"){
        covbase <- "/net/bipolar/jedidiah/mutation/reference_data/late_rt"
        dnmstmp$inside <- ifelse(dnmstmp$TIME<=-1.25, 1, 0)
      } else if(dir=="Up"){
        covbase <- "/net/bipolar/jedidiah/mutation/reference_data/early_rt"
        dnmstmp$inside <- ifelse(dnmstmp$TIME>=1.25, 1, 0)
      }
    } else if(cov=="RR"){
      covbase <- "/net/bipolar/jedidiah/mutation/reference_data/high_rr"
      dnmstmp$RR <- rcrCol(dnmstmp,
        "/net/bipolar/jedidiah/mutation/reference_data/recomb_rate.bed")
      dnmstmp$inside <- ifelse(dnmstmp$RR>=2, 1, 0)
    } else if(cov=="DHS"){
      covbase <- "/net/bipolar/jedidiah/mutation/reference_data/DHS"
      covbed <- paste0(covbase, ".bed")
      dnmstmp$inside <- binaryCol(dnmstmp, covbed)
    } else if(cov=="CpGI"){
      covbase <- "/net/bipolar/jedidiah/mutation/reference_data/cpg_islands_sorted"
      covbed <- paste0(covbase, ".bed")
      dnmstmp$inside <- binaryCol(dnmstmp, covbed)
    } else if(cov=="LAMIN"){
      covbase <- "/net/bipolar/jedidiah/mutation/reference_data/lamin_B1_LADS2"
      covbed <- paste0(covbase, ".bed")
      dnmstmp$inside <- binaryCol(dnmstmp, covbed)
    } else {
      covbase <- paste0("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/broad/sort.E062-", cov)
      covbed <- paste0(covbase, ".bed")
      dnmstmp$inside <- binaryCol(dnmstmp, covbed)
    }

    obs <- sum(dnmstmp$inside)

    seqs <- unlist(c(covdir$Sequence, lapply(covdir$Sequence, revcomp)))
    write.table(seqs, "/net/bipolar/jedidiah/mutation/seqs.txt", col.names=F, row.names=F, quote=F, sep="\t")

    grepcmd <- paste0("grep -o -Ff /net/bipolar/jedidiah/mutation/seqs.txt ", covbase, ".fa | sort | uniq -c > /net/bipolar/jedidiah/mutation/testcounts.txt")
    system(grepcmd)

    motifcts <- read.table("/net/bipolar/jedidiah/mutation/testcounts.txt", header=F, stringsAsFactors=F)

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

# motiffile <- "/net/bipolar/jedidiah/mutation/output/7bp_1000k_rates.txt"
# motifdat <- read.table(motiffile, header=T, stringsAsFactors=F)
# motifdat <- motifdat %>%
#   mutate(Category=gsub("cpg_", "", Category2)) %>%
#   mutate(Sequence=substr(Sequence, 0, 7)) %>%
#   dplyr::select(Category, Sequence, rel_prop)

covdat <- read.table("/net/bipolar/jedidiah/mutation/fa_motifs.txt", header=T, stringsAsFactors=F)
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

td2<-testdat %>% gather(group, value, c(Observed:Expected))
td3 <- data.frame()
for(i in 1:nrow(td2)){
  row <- td2[i,]
  ci <- prop.test(row$value, row$n)$conf.int[1:2]*row$n
  names(ci) <- c("lo", "hi")
  if(row$dir=="Up"){
    binom.pval <- binom.test(x=row$value, n=row$n, p=row$propexp, alternative="greater")$p.value
  } else {
    binom.pval <- binom.test(x=row$value, n=row$n, p=row$propexp, alternative="less")$p.value
  }

  newrow <- data.frame(c(row, ci, binom.pval=binom.pval), stringsAsFactors=F)
  td3 <- rbind(td3, newrow)
}

limits <- aes(ymax=hi, ymin=lo)

dirlabels <- c(Down="Subtypes with predicted depletion", Up="Subtypes with predicted enrichment")

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
  scale_fill_manual(values=rb[c(1:2,4:5)], labels=rep(c("Expected #DNMs under null\n(features have no effect)", "Observed #DNMs"), 2))+
  ylab("#DNMs")+
  theme_classic()+
  theme(legend.title=element_blank(),
    legend.key.size = unit(1.1, "cm"),
    strip.text.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.x = element_text(size=14, angle=45, vjust=1, hjust=1),
    axis.title.x=element_blank())+
  guides(colour=FALSE)
ggsave("/net/bipolar/jedidiah/mutation/images/covbars_full.png", width=12, height=8)

# testdat2 <- testdat %>%
#   filter(cov != "CpGI" & cov != "GC") %>%
#   gather("gp", "val", obs:exp) %>%
#   mutate(gp=ifelse(gp=="exp", "Expected", "Observed"),
#     dir=ifelse(dir=="Up", "Predicted enrichment", "Predicted depletion"),
#     sigalpha=ifelse(pval<0.1, ifelse(pval<0.05, 1, 0.75), 0.6))
#
# ggplot(testdat2, aes(x=gp, y=val, fill=gp, group=dir, alpha=factor(sigalpha)))+
#   geom_bar(stat="identity", position="dodge")+
#   facet_grid(dir~cov, space="free", drop=TRUE)+
#   scale_fill_brewer("Number of mutations in feature-associated subtypes", palette="Dark2")+
#   scale_alpha_discrete("p-value", labels=c("n.s.", "<0.1", "<0.05"))+
#   theme_bw()+
#   theme(axis.title.x=element_blank(),
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank(),
#     axis.title.y=element_blank(),
#     # legend.title=element_text(),
#     legend.position="bottom")
# ggsave("/net/bipolar/jedidiah/mutation/images/covdat.png", width=10, height=5)
