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
  mutate(Estmax=ifelse(Cov=="CpGI", ifelse(dir=="Up", 10, 0.125), Estmax),
    label=n) %>%
  filter(n>4)
plotcts$Category <- factor(plotcts$Category, levels=orderedcats2)
plotcts$Cov <- factor(plotcts$Cov, levels=orderedcovs)

theme_coef <- function(base_size = 12, base_family = ""){
  theme_bw(base_size = base_size) %+replace%
  theme(legend.position="bottom",
    legend.title=element_blank(),
    strip.text=element_text(size=12),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.text.y=element_text(size=10),
    axis.title.y=element_text(size=12))
}

svglite(paste0(parentdir, "/images/coef_violin2.svg"), width=7, height=7)
plotdat %>%
  filter(Cov!="CpGI") %>%
  ggplot(aes(x=Category, y=exp(Est), alpha=factor(dir), fill=Category))+
    geom_hline(yintercept=1, linetype="dashed")+
    geom_text(data=plotcts[plotcts$Cov!="CpGI",],
      aes(x=Category, y=Estmax, label=label, colour=Category),
      vjust=1, size=4, angle=90)+
    geom_violin(position="identity", scale="area")+
    scale_fill_manual(values=cols, drop=FALSE)+
    scale_colour_manual(values=cols, drop=FALSE)+
    scale_x_discrete(drop=FALSE)+
    # scale_y_log10(breaks=c(.125, .25, 0.5, 1, 2, 4, 8))+
    scale_y_log10(breaks=c(0.5, 1, 2, 4))+
    scale_alpha_discrete(range = c(0.95, 0.96), guide=F)+
    facet_wrap(~Cov, ncol=4, drop=F)+
    ylab("odds ratio for mutability")+
    guides(fill = guide_legend(nrow = 3))+
    theme_coef()
dev.off()

svglite(paste0(parentdir, "/images/coef_violin2_cpgi2.svg"), width=7, height=7)
plotdat %>%
  ggplot(aes(x=Category, y=exp(Est), alpha=factor(dir), fill=Category))+
    geom_hline(yintercept=1, linetype="dashed")+
    geom_text(data=plotcts[plotcts$Cov=="CpGI",],
      aes(x=Category, y=Estmax, label=label, colour=Category),
      vjust=1, angle=90, size=4)+
    geom_violin(position="identity", scale="area")+
    scale_fill_manual(values=cols, drop=FALSE)+
    scale_colour_manual(values=cols, drop=FALSE)+
    scale_x_discrete(drop=FALSE)+
    scale_y_log10(breaks=c(.125, .25, 0.5, 1, 2, 4, 8),
      labels=c(.125, .25, 0.5, 1, 2, 4, 8))+
    scale_alpha_discrete(range = c(0.95, 0.96), guide=F)+
    facet_wrap(~Cov, ncol=4, drop=F)+
    ylab("odds ratio for mutability")+
    guides(fill = guide_legend(nrow = 3))+
    theme_coef()
dev.off()

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
    guides(fill = guide_legend(nrow = 3))+
    theme_coef()
ggsave(paste0(parentdir, "/images/coef_violin_full.png"), width=9, height=8)

##############################################################################
# Interpolate GC effects for Jun
##############################################################################
plotdat_gc <- plotdat %>%
  filter(Cov=="GC" & Category=="A>G") %>%
  mutate(Odds=exp(Est), subtype=paste0(Category, "_", Sequence)) %>%
  dplyr::select(Category, Sequence, subtype, Odds) %>%
  filter(Odds<2)

plotdat_exp <- plotdat_gc %>%
  tidyr::expand(subtype, ind=seq(0:10)-1)

plotdat_gc <- merge(plotdat_exp, plotdat_gc, by="subtype") %>%
  mutate(Odds=Odds^ind, pctGC=ind*10, seq3=substr(Sequence, 3, 5))

# Get top 6 positive/negative effects
top <- plotdat_gc %>%
  ungroup() %>%
  filter(pctGC==100) %>%
  top_n(6, Odds) %>%
  dplyr::select(subtype)

bottom <- plotdat_gc %>%
  ungroup() %>%
  filter(pctGC==100) %>%
  top_n(-6, Odds) %>%
  dplyr::select(Sequence)

tb <- rbind(top, bottom)

plotdat_gc <- plotdat_gc %>%
  filter(subtype %in% tb$subtype)

ggplot(plotdat_gc, aes(x=pctGC, y=log(Odds), colour=subtype, group=subtype))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=0)+
  scale_x_continuous(expand=c(0,0.4),
    breaks=seq(0, 100, 10),
    labels=seq(0, 100, 10))+
  scale_colour_manual(values=iwhPalette[1:12])+
  ylab("log-odds of mutability")+
  xlab("%GC")+
  theme_bw()+
  theme(legend.position="bottom")

ggsave(paste0(parentdir, "/images/gc_effect.png"))

sitefile1 <- paste0(parentdir, "/output/logmod_data/motifs3/CGTAATT.txt")
sitefile2 <- paste0(parentdir, "/output/logmod_data/motifs3/TCGATTG.txt")

sites1 <- read.table(sitefile1, header=F, stringsAsFactors=F)
names(sites1) <- c("CHR", "POS", "Sequence", mut_cats, "DP")

sites1$GC <- gcCol(sites1,
  paste0(parentdir, "/reference_data/gc10kb.bed"))

sites2 <- read.table(sitefile2, header=F, stringsAsFactors=F)
names(sites2) <- c("CHR", "POS", "Sequence", mut_cats, "DP")

sites2$GC <- gcCol(sites2,
  paste0(parentdir, "/reference_data/gc10kb.bed"))

sites1$GP <- "CGT[A>G]ATT (+)"
sites2$GP <- "TCG[A>G]TTG (-)"

mround <- function(x,base){
  base*round(x/base)
}

# sites2 <- sites2 %>%
rbind(sites1, sites2) %>%
mutate(GC=mround(GC*100, 0.5)) %>%
group_by(GP, GC, AT_GC) %>%
tally() %>%
spread(AT_GC, n) %>%
setNames(c("GP", "GC", "nonmut", "nERVs")) %>%
na.omit() %>%
mutate(tot=nERVs+nonmut, pctmut=nERVs/tot) %>% #data.frame
filter(nERVs>2) %>%
ggplot(aes(x=GC, y=pctmut, colour=GP))+
  geom_point(aes(group=GP, size=nERVs))+
  geom_smooth(method="lm", aes(weight=tot))+
  facet_wrap(~GP)+
  # geom_line()+
  scale_colour_manual(values=iwhPalette[c(6,10)])+
  xlab("%GC")+
  ylab("ERV fraction")+
  theme_bw()
ggsave(paste0(parentdir, "/images/gc_effect_raw.png"))

# rbind(sites1, sites2) %>%
# # mutate(GC=mround(GC*100, 2)) %>%
# group_by(GP, GC, AT_GC) %>%
# filter(GC>0.3 & GC<0.6) %>%
# # tally() %>%
# # spread(AT_GC, n) %>%
# # setNames(c("GP", "GC", "nonmut", "mut")) %>%
# # na.omit() %>%
# # mutate(tot=mut+nonmut, pctmut=mut/tot) %>% #data.frame
# # filter(mut>3) %>%
# ggplot(aes(x=AT_GC, y=GC, fill=GP, group=AT_GC))+
#   # geom_point(alpha=0.3, shape=21)+
#   geom_violin()+
#   facet_wrap(~GP)+
#   # geom_line()+
#   scale_colour_manual(values=iwhPalette[c(6,10)])+
#   xlab("%GC")+
#   ylab("mutated")+
#   theme_bw()
# ggsave(paste0(parentdir, "/images/gc_effect_raw2.png"), width=8, height=6)

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
  } else {
    covdir <- covtmp %>% filter(Est<0)
  }

  dnmstmp <- chrpfdnm[paste0(chrpdnm$Category, "_", substr(chrpfdnm$SEQ, 0, 7)) %in%
    paste0(covdir$Category, "_", covdir$Sequence),]
  if(nrow(dnmstmp)>=5){
    if(cov=="GC"){
      covbase <- paste0(parentdir, "/reference_data/high_gc")
      covbed <- paste0(parentdir, "/reference_data/gc10kb.bed")
      dnmstmp$GC <- gcCol(dnmstmp, covbed)
      dnmstmp$inside <- ifelse(dnmstmp$GC>=0.55, 1, 0)
      awkstr <- "awk -F\"\\t\" '{if($4>=0.55 && NR>1) print }'"
    } else if(cov=="TIME"){
      covbed <- paste0(parentdir, "/reference_data/lymph_rep_time.txt")
      dnmstmp$TIME <- repCol(dnmstmp, covbed)
      if(dir=="Down"){
        covbase <- paste0(parentdir, "/reference_data/late_rt")
        dnmstmp$inside <- ifelse(dnmstmp$TIME<=-1.25, 1, 0)
        awkstr <- "awk -F\"\\t\" '{if($3<=-1.25 && NR>1) print $1\"\\t\"$2-500+1\"\\t\"$2+500}'"
      } else if(dir=="Up"){
        covbase <- paste0(parentdir, "/reference_data/early_rt")
        dnmstmp$inside <- ifelse(dnmstmp$TIME>1.25, 1, 0)
        awkstr <- "awk -F\"\\t\" '{if($3>1.25 && NR>1) print $1\"\\t\"$2-500+1\"\\t\"$2+500}'"
      }
    } else if(cov=="RR"){
      covbase <- paste0(parentdir, "/reference_data/high_rr")
      covbed <- paste0(parentdir, "/reference_data/recomb_rate.bed")
      dnmstmp$RR <- rcrCol(dnmstmp, covbed)
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
    } else {
      covbase <- paste0(parentdir,
        "/reference_data/sort.E062-", cov)
      covbed <- paste0(covbase, ".bed")
      dnmstmp$inside <- binaryCol(dnmstmp, covbed)
    }

    obs <- sum(dnmstmp$inside)

    seqs <- unlist(c(covdir$Sequence, lapply(covdir$Sequence, revcomp)))
    write.table(seqs, paste0(parentdir, "/seqs.txt"),
      col.names=F, row.names=F, quote=F, sep="\t")

    # Create fasta file from feature file
    if(!file.exists(paste0(covbase, ".fa"))){

      if(cov=="RR"){
        getfastacmd <- paste0("bedtools getfasta -fi ", parentdir,
          "/reference_data/human_g1k_v37.fasta -bed <(cut -f1-4", covbed,
          " | sed 's/chr//g' | ", awkstr,
          " | sed /^X/d | sed /^Y/d) -fo ", covbase, ".fa")

      } else if (cov=="TIME"){
        getfastacmd <- paste0("bedtools getfasta -fi ", parentdir,
          "/reference_data/human_g1k_v37.fasta -bed <(cut -f1-3", covbed,
          " | sed 's/chr//g' | ", awkstr,
          " | sed /^X/d | sed /^Y/d) -fo ", covbase, ".fa")
      } else if (cov=="GC"){
        getfastacmd <- paste0("bedtools getfasta -fi ", parentdir,
          "/reference_data/human_g1k_v37.fasta -bed <(cut -f1-5", covbed,
          " | sed 's/chr//g' | ", awkstr,
          " | sed /^X/d | sed /^Y/d) -fo ", covbase, ".fa")
      } else {
        getfastacmd <- paste0("bedtools getfasta -fi ", parentdir,
          "/reference_data/human_g1k_v37.fasta -bed <(sed 's/chr//g' ", covbed,
          " | sed /^X/d | sed /^Y/d) -fo ", covbase, ".fa")
      }

      system(getfastacmd)
    }

    # Search feature fasta for motifs in list
    grepcmd <- paste0("grep -o -Ff ",
      parentdir, "/seqs.txt ", covbase, ".fa | sort | uniq -c > ",
      parentdir, "/testcounts.txt")
    system(grepcmd)

    motifcts <- read.table(paste0(parentdir, "/testcounts.txt"),
      header=F, stringsAsFactors=F)

    names(motifcts) <- c("Count", "SEQ")
    motifcts$REVSEQ <- unlist(lapply(motifcts$SEQ, revcomp))
    motifcts$Sequence <- ifelse(substr(motifcts$SEQ,4,4) %in% c("A", "C"),
      motifcts$SEQ, motifcts$REVSEQ)

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
