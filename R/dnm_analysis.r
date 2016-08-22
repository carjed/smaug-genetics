motiffile <- "/net/bipolar/jedidiah/mutation/output/7bp_1000k_rates.txt"
motifdat <- read.table(motiffile, header=T, stringsAsFactors=F)
motifdat <- motifdat %>%
  mutate(Category=gsub("cpg_", "", Category2)) %>%
  mutate(Sequence=substr(Sequence, 0, 7)) %>%
  dplyr::select(Category, Sequence, rel_prop)

covdat <- read.table("/net/bipolar/jedidiah/mutation/fa_motifs.txt", header=T, stringsAsFactors=F)
# covdat$Category <- gsub("cpg_", "", covdat$Category)
covs <- unique(covdat$Cov)
covs <- covs[grepl("H3", covs)]
testdat <- data.frame()
for(i in 1:length(covs)){
# for(i in 1:2){
  cov <- covs[i]
  covtmp <- covdat %>% filter(Cov==cov)
  covup <- covtmp %>% filter(Est > 0)
  covdown <- covtmp %>% filter(Est < 0)

  dnmstmpup <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in%
    paste0(covup$Category, "_", covup$Sequence),]
  dnmstmpdown <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in%
    paste0(covdown$Category, "_", covdown$Sequence),]

  # fastacmd <- paste0("bash bedtools getfasta -fi /net/bipolar/jedidiah/mutation/reference_data/human_g1k_v37.fasta -bed <(sed \'s/chr//g\' /net/bipolar/jedidiah/mutation/reference_data/histone_marks/broad/sort.E062-", cov, ".bed | sed /^X/d | sed /^Y/d) -fo /net/bipolar/jedidiah/mutation/reference_data/histone_marks/broad/sort.E062-", cov, ".fa")
  # system(fastacmd)

  covfile <- paste0("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/broad/sort.E062-", cov, ".bed")
  covbed <- read.table(covfile, header=F)
  covsum <- covbed %>%
    filter(!grepl("X|Y", V1)) %>%
    mutate(length=V3-V2) %>%
    summarise(total=sum(length), prop=total/3e9)

  if(nrow(dnmstmpup)>=5){
    dnmstmpup$inside <- binaryCol(dnmstmpup, covfile)

    upobs <- sum(dnmstmpup$inside)
    # upexp <- covsum$prop*nrow(dnmstmpup)

    seqsup <- unlist(c(covup$Sequence, lapply(covup$Sequence, revcomp)))
    write.table(seqsup, "/net/bipolar/jedidiah/mutation/seqsup.txt", col.names=F, row.names=F, quote=F, sep="\t")

    grepcmd <- paste0("grep -o -Ff /net/bipolar/jedidiah/mutation/seqsup.txt /net/bipolar/jedidiah/mutation/reference_data/histone_marks/broad/sort.E062-", cov, ".fa | sort | uniq -c > /net/bipolar/jedidiah/mutation/testcounts.txt")
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
    covup2 <- merge(covup, motifcts2, by=c("Sequence"))

    covup3 <- merge(covup2, motifdat, by=c("Category", "Sequence"))
    covup3$exp <- covup3$Count*covup3$rel_prop*1.67e-06*1074
    upexp <- sum(covup3$exp)

    uptotal <- nrow(dnmstmpup)
    uptest <- prop.test(c(upobs, upexp), c(uptotal, uptotal))
    uprow <- data.frame(cov, dir="up", obs=upobs, exp=upexp, n=uptotal,
      propobs=uptest$estimate[1], propexp=uptest$estimate[2], pval=uptest$p.value)
    testdat <- rbind(testdat, uprow)
  }

  if(nrow(dnmstmpdown)>=5){
    dnmstmpdown$inside <- binaryCol(dnmstmpdown, covfile)

    downobs <- sum(dnmstmpdown$inside)
    seqsdown <- unlist(c(covdown$Sequence, lapply(covdown$Sequence, revcomp)))
    write.table(seqsdown, "/net/bipolar/jedidiah/mutation/seqsdown.txt", col.names=F, row.names=F, quote=F, sep="\t")

    grepcmd <- paste0("grep -o -Ff /net/bipolar/jedidiah/mutation/seqsdown.txt /net/bipolar/jedidiah/mutation/reference_data/histone_marks/broad/sort.E062-", cov, ".fa | sort | uniq -c > /net/bipolar/jedidiah/mutation/testcounts.txt")
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
    covdown2 <- merge(covdown, motifcts2, by=c("Sequence"))

    covdown3 <- merge(covdown2, motifdat, by=c("Category", "Sequence"))
    covdown3$exp <- covdown3$Count*covdown3$rel_prop*1.67e-06*1074
    downexp <- sum(covdown3$exp)

    downtotal <- nrow(dnmstmpdown)
    downtest <- prop.test(c(downobs, downexp), c(downtotal, downtotal))
    downrow <- data.frame(cov, dir="down", obs=downobs, exp=downexp, n=downtotal,
      propobs=downtest$estimate[1], propexp=downtest$estimate[2], pval=downtest$p.value)
    testdat <- rbind(testdat, downrow)
  }
}

# Repeat with CpGI
cov <- "CpGI"
covtmp <- covdat %>% filter(Cov==cov)
covup <- covtmp %>% filter(Est > 0)
covdown <- covtmp %>% filter(Est < 0)

dnmstmpup <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in% paste0(covup$Category, "_", covup$Sequence),]
dnmstmpdown <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in% paste0(covdown$Category, "_", covdown$Sequence),]

covfile <- paste0("/net/bipolar/jedidiah/mutation/reference_data/cpg_islands_sorted.bed")
covbed <- read.table(covfile, header=F)
covsum <- covbed %>%
  filter(!grepl("X|Y", V1)) %>%
  mutate(length=V4) %>%
  summarise(total=sum(length), prop=total/3e9)

if(nrow(dnmstmpup)>=5){
  dnmstmpup$inside <- binaryCol(dnmstmpup, covfile)

  upobs <- sum(dnmstmpup$inside)
  seqsup <- unlist(c(covup$Sequence, lapply(covup$Sequence, revcomp)))
  write.table(seqsup, "/net/bipolar/jedidiah/mutation/seqsup.txt", col.names=F, row.names=F, quote=F, sep="\t")

  grepcmd <- paste0("grep -o -Ff /net/bipolar/jedidiah/mutation/seqsup.txt /net/bipolar/jedidiah/mutation/reference_data/cpg_islands_sorted.fa | sort | uniq -c > /net/bipolar/jedidiah/mutation/testcounts.txt")
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
  covup2 <- merge(covup, motifcts2, by=c("Sequence"))

  covup3 <- merge(covup2, motifdat, by=c("Category", "Sequence"))
  covup3$exp <- covup3$Count*covup3$rel_prop*1.67e-06*1074
  upexp <- sum(covup3$exp)
  uptotal <- nrow(dnmstmpup)
  uptest <- prop.test(c(upobs, upexp), c(uptotal, uptotal))
  uprow <- data.frame(cov, dir="up", obs=upobs, exp=upexp, n=uptotal,
    propobs=uptest$estimate[1], propexp=uptest$estimate[2], pval=uptest$p.value)
  testdat <- rbind(testdat, uprow)
}

if(nrow(dnmstmpdown)>=5){
  dnmstmpdown$inside <- binaryCol(dnmstmpdown, covfile)

  downobs <- sum(dnmstmpdown$inside)
  seqsdown <- unlist(c(covdown$Sequence, lapply(covdown$Sequence, revcomp)))
  write.table(seqsdown, "/net/bipolar/jedidiah/mutation/seqsdown.txt", col.names=F, row.names=F, quote=F, sep="\t")

  grepcmd <- paste0("grep -o -Ff /net/bipolar/jedidiah/mutation/seqsdown.txt /net/bipolar/jedidiah/mutation/reference_data/cpg_islands_sorted.fa | sort | uniq -c > /net/bipolar/jedidiah/mutation/testcounts.txt")
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
  covdown2 <- merge(covdown, motifcts2, by=c("Sequence"))

  covdown3 <- merge(covdown2, motifdat, by=c("Category", "Sequence"))
  covdown3$exp <- covdown3$Count*covdown3$rel_prop*1.67e-06*1074
  downexp <- sum(covdown3$exp)
  downtotal <- nrow(dnmstmpdown)
  downtest <- prop.test(c(downobs, downexp), c(downtotal, downtotal))
  downrow <- data.frame(cov, dir="down", obs=downobs, exp=downexp, n=downtotal,
    propobs=downtest$estimate[1], propexp=downtest$estimate[2], pval=downtest$p.value)
  testdat <- rbind(testdat, downrow)
}

# Repeat with DHS
cov <- "DHS"
covtmp <- covdat %>% filter(Cov==cov)
covup <- covtmp %>% filter(Est > 0)
covdown <- covtmp %>% filter(Est < 0)

dnmstmpup <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in% paste0(covup$Category, "_", covup$Sequence),]
dnmstmpdown <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in% paste0(covdown$Category, "_", covdown$Sequence),]

covfile <- paste0("/net/bipolar/jedidiah/mutation/reference_data/DHS.bed")
covbed <- read.table(covfile, header=F)
covsum <- covbed %>%
  filter(!grepl("X|Y", V1)) %>%
  mutate(length=V3-V2) %>%
  summarise(total=sum(length), prop=total/3e9)

if(nrow(dnmstmpup)>=5){
  dnmstmpup$inside <- binaryCol(dnmstmpup, covfile)

  upobs <- sum(dnmstmpup$inside)
  seqsup <- unlist(c(covup$Sequence, lapply(covup$Sequence, revcomp)))
  write.table(seqsup, "/net/bipolar/jedidiah/mutation/seqsup.txt", col.names=F, row.names=F, quote=F, sep="\t")

  grepcmd <- paste0("grep -o -Ff /net/bipolar/jedidiah/mutation/seqsup.txt /net/bipolar/jedidiah/mutation/reference_data/DHS.fa | sort | uniq -c > /net/bipolar/jedidiah/mutation/testcounts.txt")
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
  covup2 <- merge(covup, motifcts2, by=c("Sequence"))

  covup3 <- merge(covup2, motifdat, by=c("Category", "Sequence"))
  covup3$exp <- covup3$Count*covup3$rel_prop*1.67e-06*1074
  upexp <- sum(covup3$exp)
  uptotal <- nrow(dnmstmpup)
  uptest <- prop.test(c(upobs, upexp), c(uptotal, uptotal))
  uprow <- data.frame(cov, dir="up", obs=upobs, exp=upexp, n=uptotal,
    propobs=uptest$estimate[1], propexp=uptest$estimate[2], pval=uptest$p.value)
  testdat <- rbind(testdat, uprow)
}

if(nrow(dnmstmpdown)>=5){
  dnmstmpdown$inside <- binaryCol(dnmstmpdown, covfile)

  downobs <- sum(dnmstmpdown$inside)
  seqsdown <- unlist(c(covdown$Sequence, lapply(covdown$Sequence, revcomp)))
  write.table(seqsdown, "/net/bipolar/jedidiah/mutation/seqsdown.txt", col.names=F, row.names=F, quote=F, sep="\t")

  grepcmd <- paste0("grep -o -Ff /net/bipolar/jedidiah/mutation/seqsdown.txt /net/bipolar/jedidiah/mutation/reference_data/DHS.fa | sort | uniq -c > /net/bipolar/jedidiah/mutation/testcounts.txt")
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
  covdown2 <- merge(covdown, motifcts2, by=c("Sequence"))

  covdown3 <- merge(covdown2, motifdat, by=c("Category", "Sequence"))
  covdown3$exp <- covdown3$Count*covdown3$rel_prop*1.67e-06*1074
  downexp <- sum(covdown3$exp)
  downtotal <- nrow(dnmstmpdown)
  downtest <- prop.test(c(downobs, downexp), c(downtotal, downtotal))
  downrow <- data.frame(cov, dir="down", obs=downobs, exp=downexp, n=downtotal,
    propobs=downtest$estimate[1], propexp=downtest$estimate[2], pval=downtest$p.value)
  testdat <- rbind(testdat, downrow)
}

# Repeat with LAD
cov <- "LAMIN"
covtmp <- covdat %>% filter(Cov==cov)
covup <- covtmp %>% filter(Est > 0)
covdown <- covtmp %>% filter(Est < 0)

dnmstmpup <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in% paste0(covup$Category, "_", covup$Sequence),]
dnmstmpdown <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in% paste0(covdown$Category, "_", covdown$Sequence),]

covfile <- paste0("/net/bipolar/jedidiah/mutation/reference_data/lamin_B1_LADS2.bed")
covbed <- read.table(covfile, header=F)
covsum <- covbed %>%
  filter(!grepl("X|Y", V1)) %>%
  mutate(length=V3-V2) %>%
  summarise(total=sum(length), prop=total/3e9)

if(nrow(dnmstmpup)>=5){
  dnmstmpup$inside <- binaryCol(dnmstmpup, covfile)

  upobs <- sum(dnmstmpup$inside)
  seqsup <- unlist(c(covup$Sequence, lapply(covup$Sequence, revcomp)))
  write.table(seqsup, "/net/bipolar/jedidiah/mutation/seqsup.txt", col.names=F, row.names=F, quote=F, sep="\t")

  grepcmd <- paste0("grep -o -Ff /net/bipolar/jedidiah/mutation/seqsup.txt /net/bipolar/jedidiah/mutation/reference_data/lamin_B1_LADS2.fa | sort | uniq -c > /net/bipolar/jedidiah/mutation/testcounts.txt")
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
  covup2 <- merge(covup, motifcts2, by=c("Sequence"))

  covup3 <- merge(covup2, motifdat, by=c("Category", "Sequence"))
  covup3$exp <- covup3$Count*covup3$rel_prop*1.67e-06*1074
  upexp <- sum(covup3$exp)
  uptotal <- nrow(dnmstmpup)
  uptest <- prop.test(c(upobs, upexp), c(uptotal, uptotal))
  uprow <- data.frame(cov, dir="up", obs=upobs, exp=upexp, n=uptotal,
    propobs=uptest$estimate[1], propexp=uptest$estimate[2], pval=uptest$p.value)
  testdat <- rbind(testdat, uprow)
}

if(nrow(dnmstmpdown)>=5){
  dnmstmpdown$inside <- binaryCol(dnmstmpdown, covfile)

  downobs <- sum(dnmstmpdown$inside)
  seqsdown <- unlist(c(covdown$Sequence, lapply(covdown$Sequence, revcomp)))
  write.table(seqsdown, "/net/bipolar/jedidiah/mutation/seqsdown.txt", col.names=F, row.names=F, quote=F, sep="\t")

  grepcmd <- paste0("grep -o -Ff /net/bipolar/jedidiah/mutation/seqsdown.txt /net/bipolar/jedidiah/mutation/reference_data/lamin_B1_LADS2.fa | sort | uniq -c > /net/bipolar/jedidiah/mutation/testcounts.txt")
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
  covdown2 <- merge(covdown, motifcts2, by=c("Sequence"))

  covdown3 <- merge(covdown2, motifdat, by=c("Category", "Sequence"))
  covdown3$exp <- covdown3$Count*covdown3$rel_prop*1.67e-06*1074
  downexp <- sum(covdown3$exp)
  downtotal <- nrow(dnmstmpdown)
  downtest <- prop.test(c(downobs, downexp), c(downtotal, downtotal))
  downrow <- data.frame(cov, dir="down", obs=downobs, exp=downexp, n=downtotal,
    propobs=downtest$estimate[1], propexp=downtest$estimate[2], pval=downtest$p.value)
  testdat <- rbind(testdat, downrow)
}

# recombination rate
{
  cov <- "RR"
  covtmp <- covdat %>% filter(Cov==cov)
  covup <- covtmp %>% filter(Est > 0)
  covdown <- covtmp %>% filter(Est < 0)

  dnmstmpup <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in%
    paste0(covup$Category, "_", covup$Sequence),]
  dnmstmpupinv <- chrpfdnm[!(paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in%
    paste0(covup$Category, "_", covup$Sequence)),]

  dnmstmpup$RR <- rcrCol(dnmstmpup,
    "/net/bipolar/jedidiah/mutation/reference_data/recomb_rate.bed")

  dnmstmpupinv$RR <- rcrCol(dnmstmpupinv,
    "/net/bipolar/jedidiah/mutation/reference_data/recomb_rate.bed")
  rcruptest <- t.test(dnmstmpup$RR, dnmstmpupinv$RR)
}

# replication timing
{
  cov <- "TIME"
  covtmp <- covdat %>% filter(Cov==cov)
  covup <- covtmp %>% filter(Est > 0)
  covdown <- covtmp %>% filter(Est < 0)

  dnmstmpup <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in%
    paste0(covup$Category, "_", covup$Sequence),]
  dnmstmpupinv <- chrpfdnm[!(paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in%
    paste0(covup$Category, "_", covup$Sequence)),]

  dnmstmpup$TIME <- repCol(dnmstmpup,
    "/net/bipolar/jedidiah/mutation/reference_data/lymph_rep_time.txt")

  dnmstmpupinv$TIME <- repCol(dnmstmpupinv,
    "/net/bipolar/jedidiah/mutation/reference_data/lymph_rep_time.txt")

  reptimeuptest <- t.test(dnmstmpup$TIME, dnmstmpupinv$TIME)

  dnmstmpdown <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in%
    paste0(covdown$Category, "_", covdown$Sequence),]
  dnmstmpdowninv <- chrpfdnm[!(paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in%
    paste0(covdown$Category, "_", covdown$Sequence)),]

  dnmstmpdown$TIME <- repCol(dnmstmpdown,
    "/net/bipolar/jedidiah/mutation/reference_data/lymph_rep_time.txt")

  dnmstmpdowninv$TIME <- repCol(dnmstmpdowninv,
    "/net/bipolar/jedidiah/mutation/reference_data/lymph_rep_time.txt")

  reptimedowntest <- t.test(dnmstmpdown$TIME, dnmstmpdowninv$TIME)
}

# gc content
{
  cov <- "GC"
  covtmp <- covdat %>% filter(Cov==cov)
  covup <- covtmp %>% filter(Est > 0)
  covdown <- covtmp %>% filter(Est < 0)

  dnmstmpup <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in%
    paste0(covup$Category, "_", covup$Sequence),]
  dnmstmpupinv <- chrpfdnm[!(paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in%
    paste0(covup$Category, "_", covup$Sequence)),]

  dnmstmpup$GC <- gcCol(dnmstmpup,
    "/net/bipolar/jedidiah/mutation/output/3bp_10k/full_bin.txt")

  dnmstmpupinv$GC <- gcCol(dnmstmpupinv,
    "/net/bipolar/jedidiah/mutation/output/3bp_10k/full_bin.txt")

  gcuptest <- t.test(dnmstmpup$GC, dnmstmpupinv$GC)

  dnmstmpdown <- chrpfdnm[paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in%
    paste0(covdown$Category, "_", covdown$Sequence),]
  dnmstmpdowninv <- chrpfdnm[!(paste0(chrpfdnm$Category.x, "_", chrpfdnm$Sequence) %in%
    paste0(covdown$Category, "_", covdown$Sequence)),]

  dnmstmpdown$GC <- gcCol(dnmstmpdown,
    "/net/bipolar/jedidiah/mutation/output/3bp_10k/full_bin.txt")

  dnmstmpdowninv$GC <- gcCol(dnmstmpdowninv,
    "/net/bipolar/jedidiah/mutation/output/3bp_10k/full_bin.txt")

  gcdowntest <- t.test(dnmstmpdown$GC, dnmstmpdowninv$GC)
}
