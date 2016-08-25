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

runTest <- function(cov, dir){
  covtmp <- covdat %>% filter(Cov==cov)
  if(dir=="Up"){
    covdir <- covtmp %>% filter(Est>0)
  } else {
    covdir <- covtmp %>% filter(Est<0)
  }

  dnmstmp <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in%
    paste0(covdir$Category, "_", covdir$Sequence),]
  if(nrow(dnmstmpup)>=5){
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
    test <- prop.test(c(obs, exp), c(total, total))
    newrow <- data.frame(cov, dir=dir, obs=obs, exp=exp, n=total,
      propobs=test$estimate[1], propexp=test$estimate[2], pval=test$p.value)
    # testdat <- rbind(testdat, newrow)
    newrow
  }
}

motiffile <- "/net/bipolar/jedidiah/mutation/output/7bp_1000k_rates.txt"
motifdat <- read.table(motiffile, header=T, stringsAsFactors=F)
motifdat <- motifdat %>%
  mutate(Category=gsub("cpg_", "", Category2)) %>%
  mutate(Sequence=substr(Sequence, 0, 7)) %>%
  dplyr::select(Category, Sequence, rel_prop)

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
