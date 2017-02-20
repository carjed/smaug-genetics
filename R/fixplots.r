      axis.title.y=element_text(size=16),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=14))+
  ylab(expression(paste("Fraction of variance explained (Nagelkerke ", R^2, ")", sep="")))
ggsave(paste0(parentdir, "/images/EvC_rsq_full.png"), width=8, height=6)
ptm <- proc.time()
cat("Analyzing genomic features...\n")
source("./R/coef_summary.r")
tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
require(svglite)
ptm <- proc.time()
cat("Analyzing genomic features...\n")
source("./R/coef_summary.r")
tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
head(coefs)
dim(coefs)
length(unique(paste0(coefs$Sequence, "_", coefs$Category)))
head(rates_full)
ls()
head(r5m)
r5m1<-r5m %>% mutate(nCG=str_count(Sequence, c("C","G")))
r5m1<-r5m %>% mutate(nCG=str_count(Sequence, "C"))
head(r5m1)
r5m1<-r5m %>% mutate(nC=str_count(paste0(v2, v3), "C"), gp=ifelse(nCG<4, 1, 0))
head(r5m1)
head(r5m1, 20)
r5m1<-r5m %>%
  mutate(nC=str_count(paste0(v2, v3), "C"),
    nG=str_count(paste0(v2, v3), "G"),
    nCG=nC+nG,
    gp=ifelse(nCG<4, 1, 0)) %>%
    group_by(Category2) %>%
    spread(gp, prop_diff5)
head(r5m1)
dim(r5m1)
data.frame(head(r5m1))
r5m1<-r5m %>%
  mutate(nC=str_count(paste0(v2, v3), "C"),
    nG=str_count(paste0(v2, v3), "G"),
    nCG=nC+nG,
    gp=ifelse(nCG<4, "low", "high")) %>%
    group_by(Category2) %>%
    spread(gp, prop_diff5)
mean(r5m1$low)
mean(r5m1$low, na.rm=T)
mean(r5m1$high, na.rm=T)
r5m1<-r5m %>%
  mutate(nC=str_count(paste0(v2, v3), "C"),
    nG=str_count(paste0(v2, v3), "G"),
    nCG=nC+nG,
    gp=ifelse(nCG<4, "low", "high")) %>%
    group_by(Category2) %>%
    spread(gp, prop_diff5) %>%
    summarise(pval=t.test(low, high, na.omit)$p.value)
r5m1<-r5m %>%
  mutate(nC=str_count(paste0(v2, v3), "C"),
    nG=str_count(paste0(v2, v3), "G"),
    nCG=nC+nG,
    gp=ifelse(nCG<4, "low", "high")) %>%
    group_by(Category2) %>%
    spread(gp, prop_diff5) %>%
    summarise(pval=t.test(low, high, na.action=na.omit)$p.value)
dim(r5m1)
r5m1
r5m1<-r5m %>%
  mutate(nC=str_count(paste0(v2, v3), "C"),
    nG=str_count(paste0(v2, v3), "G"),
    nCG=nC+nG,
    gp=ifelse(nCG<4, "low", "high")) %>%
    group_by(Category2) %>%
    spread(gp, prop_diff5) %>%
    summarise(pval=t.test(log(low), log(high), na.action=na.omit)$p.value)
r5m1
r5m1<-r5m %>%
  mutate(nC=str_count(paste0(v2, v3), "C"),
    nG=str_count(paste0(v2, v3), "G"),
    nCG=nC+nG,
    gp=ifelse(nCG<4, "low", "high")) %>%
    group_by(Category2) %>%
    spread(gp, prop_diff5) %>%
    summarise(low=mean(log(low)),
      high=mean(log(high)),
      pval=t.test(log(low), log(high), na.action=na.omit)$p.value)
r5m1<-r5m %>%
  mutate(nC=str_count(paste0(v2, v3), "C"),
    nG=str_count(paste0(v2, v3), "G"),
    nCG=nC+nG,
    gp=ifelse(nCG<4, "low", "high")) %>%
    group_by(Category2) %>%
    spread(gp, prop_diff5) %>%
    summarise(low=mean(log(low), na.rm=T),
      high=mean(log(high), na.rm=T),
      pval=t.test(log(low), log(high), na.action=na.omit)$p.value)
r5m1<-r5m %>%
  mutate(nC=str_count(paste0(v2, v3), "C"),
    nG=str_count(paste0(v2, v3), "G"),
    nCG=nC+nG,
    gp=ifelse(nCG<4, "low", "high")) %>%
    group_by(Category2) %>%
    spread(gp, prop_diff5)
mean(r5m1$low, na.rm=T)
mean(log(r5m1$low), na.rm=T)
r5m1<-r5m %>%
  mutate(nC=str_count(paste0(v2, v3), "C"),
    nG=str_count(paste0(v2, v3), "G"),
    nCG=nC+nG,
    gp=ifelse(nCG<4, "low", "high")) %>%
    group_by(Category2) %>%
    spread(gp, prop_diff5) %>%
    summarise(low=mean(low, na.rm=T),
      high=mean(high, na.rm=T),
      pval=t.test(log(low), log(high), na.action=na.omit)$p.value)
r5m1
r5m1<-r5m %>%
  mutate(nC=str_count(paste0(v2, v3), "C"),
    nG=str_count(paste0(v2, v3), "G"),
    nCG=nC+nG,
    gp=ifelse(nCG<4, "low", "high")) %>%
    group_by(Category2) %>%
    spread(gp, prop_diff5) %>%
    summarise(lowm=mean(log(low), na.rm=T),
      highm=mean(log(high), na.rm=T),
      pval=t.test(log(low), log(high), na.action=na.omit)$p.value)
r5m1
exp(-0.001)
exp(.24347)
r5m1
r5m1$lowm<-exp(r5m1$lowm)
r5m1$highm<-exp(r5m1$highm)
r5m1
head(r5m)
ls()
head(rates1)
rates1
head(rates)
dim(rates)
summary(rates$num)
dim(rates[rates$num<5,])
dim(rates[rates$num<2,])
rates %>% group_by(Category2) %>% summarise(n=sum(num))
ls()
summfile
summ_5bp_100k <- read.table(summfile, header=F, stringsAsFactors=F, skip=1)
names(summ_5bp_100k)<-c(
"CHR", "POS", "REF", "ALT", "DP", "AN", "SEQ", "ALTSEQ", "GC")
gc()


##
# compare singleton counts
##
getBinomP <- function(x1,x2){
    prop <- 1789/3716
    pval <- binom.test(c(x1, x2), p=prop)$p.value
    return(pval)
}

# Gene-wise tests for significant excess of singletons in cases
gene_wise <- function(data){
  outdat <- data %>%
    group_by(BP, GENE) %>%
    # filter(CADD>20) %>%
    summarise(n=n()) %>%
    spread(BP, n)

  nIND.adj <- 1789/3716

  names(outdat)<-c("GENE", "CONTROL", "CASE")


  # outdat <- outdat[complete.cases(outdat),]
  outdat[is.na(outdat)] <- 0
  outdat$TOTAL <- outdat$CASE+outdat$CONTROL


  outdat <- outdat %>%
    filter(TOTAL>5) %>%
    rowwise() %>%
    mutate(pval=getBinomP(CASE, CONTROL)) %>%
    arrange(pval)

  outdat$qval <- p.adjust(outdat$pval, method="fdr")

  return(outdat)
}

# Read in annotated BRIDGES positions
bridges.sites <- read.table("/exports/home/jedidiah/BRIDGES.anno.sub.sites", header=F, stringsAsFactors=F)
names(bridges.sites) <- c("CHR", "POS", "ANNO", "GENE", "CADD")

# Subset BRIDGES positions to singletons
bridges.sing <- merge(bridges.sites, dat_5bp_100k$summ, c("CHR", "POS"))
names(bridges.sing)[13] <- "ID"

bridges.sing$ind <- paste0(bridges.sing$CHR, "_", bridges.sing$POS)

# Read in .ped file with bipolar status
ped <- read.table("/net/wonderland/home/mattz/BRIDGES/Assoc/BRIDGES_Assoc.20160226.ped", header=F, stringsAsFactors=F)
ped2 <- ped[,1:6]

names(ped2) <- c("FAM", "ID", "F", "M", "SEX", "BP")

# Merge bipolar status with main file
bridges.sing <- merge(bridges.sing, ped2, by="ID")

# Read in ExAC positions
exac.sites <- read.table("/exports/home/jedidiah/ExAC.sites", header=F)
names(exac.sites) <- c("CHR", "POS")
exac.sites$ind <- paste0(exac.sites$CHR, "_", exac.sites$POS)

# Read in constrained gene data (Samocha et al., 2014)
constrained <- read.table("/exports/home/jedidiah/constrained.txt", header=T, stringsAsFactors=F)
names(constrained)[2]<-"GENE"

# Read in psychiatric genes (Purcell et al., 2013)
psych.genes <- read.table("~/psych.genes", header=F, stringsAsFactors=F)
names(psych.genes) <- c("GENE")

# Indicator for status in ExAC
bridges.sing$exac <- bridges.sing$ind %in% exac.sites$ind

# Indicator for constrained genes
bridges.sing$constrained <- bridges.sing$GENE %in% constrained$GENE

# Indicator for psych genes
bridges.sing$psych <- bridges.sing$GENE %in% psych.genes$GENE

# indicator for high CADD scores
bridges.sing$CADD15 <- ifelse(as.numeric(bridges.sing$CADD)>15, TRUE, FALSE)
bridges.sing$CADD20 <- ifelse(as.numeric(bridges.sing$CADD)>20, TRUE, FALSE)

# all singletons (row 1)
all_gw <- bridges.sing %>%
  group_by(BP) %>%
  mutate(GENE="none") %>%
  gene_wise()

all <- bridges.sing %>%
  gene_wise()

all_hc_gw <- bridges.sing %>%
  group_by(BP) %>%
  filter(constrained==TRUE) %>%
  mutate(GENE="none") %>%
  gene_wise()

all_hc <- bridges.sing %>%
  filter(constrained==TRUE) %>%
  gene_wise()

all_psych_gw <- bridges.sing %>%
  group_by(BP) %>%
  filter(psych==TRUE) %>%
  mutate(GENE="none") %>%
  gene_wise()

all_psych <- bridges.sing %>%
  filter(psych==TRUE) %>%
  gene_wise()

# rare singletons (row 2)
rare_gw <- bridges.sing %>%
  group_by(BP) %>%
  filter(exac==FALSE) %>%
  mutate(GENE="none") %>%
  gene_wise()

rare_all <- bridges.sing %>%
  filter(exac==FALSE) %>%
  gene_wise()

rare_hc_gw <- bridges.sing %>%
  group_by(BP) %>%
  filter(exac==FALSE & constrained==TRUE) %>%
  mutate(GENE="none") %>%
  gene_wise()

rare_hc <- bridges.sing %>%
  filter(exac==FALSE & constrained==TRUE) %>%
  gene_wise()

rare_psych_gw <- bridges.sing %>%
  group_by(BP) %>%
  filter(exac==FALSE & psych==TRUE) %>%
  mutate(GENE="none") %>%
  gene_wise()

rare_psych <- bridges.sing %>%
  filter(exac==FALSE & psych==TRUE) %>%
  gene_wise()

# rare high-CADD singletons (row 3)

rare_c_gw <- bridges.sing %>%
  group_by(BP) %>%
  filter(exac==FALSE & CADD15==TRUE) %>%
  mutate(GENE="none") %>%
  gene_wise()

rare_c_all <- bridges.sing %>%
  # filter(exac==FALSE & CADD20==TRUE) %>%
  filter(exac==FALSE & CADD) %>%
  # filter(CADD15==TRUE) %>%
  gene_wise()

rare_c_hc_gw <- bridges.sing %>%
  group_by(BP) %>%
  filter(exac==FALSE & CADD15==TRUE & constrained==TRUE) %>%
  mutate(GENE="none") %>%
  gene_wise()

rare_c_hc <- bridges.sing %>%
  # filter(exac==FALSE & CADD15==TRUE & constrained==TRUE) %>%
  filter(as.numeric(CADD)>20 & constrained==TRUE) %>%
  gene_wise()

rare_c_psych_gw <- bridges.sing %>%
  group_by(BP) %>%
  filter(exac==FALSE & CADD15==TRUE & psych==TRUE) %>%
  mutate(GENE="none") %>%
  gene_wise()

rare_c_psych <- bridges.sing %>%
  # filter(exac==FALSE & as.numeric(CADD)>25 & psych==TRUE) %>%
  # filter(as.numeric(CADD)>30 & psych==TRUE) %>%
  filter(CADD15==TRUE & psych==TRUE & exac==FALSE) %>%
  gene_wise()


bridges.sing.rare.cadd.genes <- gene_wise(bridges.sing.rare.cadd)
bridges.sing.constrained.rare.cadd.genes <- gene_wise(bridges.sing.constrained.rare.cadd)

bridges.sing.psych.rare.cadd.genes <- gene_wise(bridges.sing.psych.rare.cadd)

spectra_3_1Mb <- dat_5bp_100k$summ %>%
  dplyr::select(CHR, BIN, Category2, Sequence) %>%
  mutate(SEQ3=substr(Sequence,2,6),
    subtype=paste0(Category2, "_", SEQ3),
    chrbin=paste0("chr", CHR, "_", BIN)) %>%
  group_by(subtype, chrbin) %>%
  summarise(n=n())

  spectra <- spectra_3_1Mb %>%
    spread(chrbin, n)

  spectra[is.na(spectra)]<-0

dat<- as.matrix(spectra[,2:2735])

sigs_nmf2<-nmf(dat, 2)
evar(sigs_nmf2, dat)

sigs_nmf3<-nmf(dat, 3)
evar(sigs_nmf3, dat)

sigs_nmf4<-nmf(dat, 4)
evar(sigs_nmf4, dat)

sigs_nmf5<-nmf(dat, 5)
evar(sigs_nmf5, dat)


  sigs_nmf<-nmf(as.matrix(spectra[,2:800]), 3)
  sigs_pca<-princomp(as.matrix(spectra[,2:800]))


anno2b.plot<-anno2b %>% group_by(GENE, BP) %>% summarise(CADD=mean(CADD))
warnings()
head(anno2b)
class(anno2b$CADD)
anno2b$CADD<-as.numeric(anno2b$CADD)
anno2b.plot<-anno2b %>% group_by(GENE, BP) %>% summarise(CADD=mean(CADD))
dim(anno2b.plot)
head(anno2b.plot)
anno2b.plot<-anno2b %>% group_by(GENE, BP) %>% summarise(CADD=mean(CADD), CONS=max(mis_z_sign))
head(anno2b.plot)
ggplot(anno2b.plot, aes(x=CONS, y=CADD, group=BP, colour=BP)+geom_point()
ggplot(anno2b.plot, aes(x=CONS, y=CADD, group=BP, colour=BP))+geom_point()
ggsave("/net/bipolar/jedidiah/mutation/images/cadd_constr.png")
ggplot(anno2b.plot, aes(x=CONS, y=CADD, group=factor(BP), colour=factor(BP)))+geom_point()
ggsave("/net/bipolar/jedidiah/mutation/images/cadd_constr.png")
ggplot(anno2b.plot, aes(x=log(CONS), y=CADD, group=factor(BP), colour=factor(BP)))+geom_point()
ggsave("/net/bipolar/jedidiah/mutation/images/cadd_constr.png")
head(anno2b.plot)
anno2b.plot<-anno2b %>% group_by(GENE, BP) %>% summarise(CADDmean=mean(CADD), CADDmax=max(CADD), CONS=max(mis_z_sign))
ggplot(anno2b.plot, aes(x=log(CONS), y=CADDmax, group=factor(BP), colour=factor(BP)))+geom_point()
ggsave("/net/bipolar/jedidiah/mutation/images/cadd_constr_max.png")
ggplot(anno2b.plot, aes(x=log(CONS), y=CADDmean, group=factor(BP), colour=factor(BP)))+geom_point()
ggsave("/net/bipolar/jedidiah/mutation/images/cadd_constr_max.png")
anno2b.plot1<-anno2b.plot %>% spread(BP)
anno2b.plot1<-anno2b.plot %>% spread(BP, "case", "control")
anno2b.plot1<-anno2b.plot %>% spread(BP, CADD)
anno2b.plot1<-anno2b.plot %>% spread(CADD, BP)
anno2b.plot %>% arrange(desc(CONS))
anno2b.plot<-anno2b %>% group_by(GENE, BP) %>% summarise(CADDmean=mean(CADD), CADDmax=max(CADD), CONS=max(mis_z_sign, n=n()))
anno2b.plot %>% arrange(desc(CONS))
anno2b.plot<-anno2b %>% group_by(GENE, BP) %>% summarise(CADDmean=mean(CADD), CADDmax=max(CADD), CONS=max(mis_z_sign), n=n())
anno2b.plot %>% arrange(desc(CONS))
8
anno2b.plot %>% group_by(BP) %>% summarise(cor=cor(CADDmean, CONS))
ggplot(anno2b.plot, aes(x=CONS, y=CADDmean, group=factor(BP), colour=factor(BP)))+geom_point()+xlab("Conservation Score")+ylab("mean CADD score")
ggsave("/net/bipolar/jedidiah/mutation/images/cadd_constr.png")
ggplot(anno2b.plot, aes(x=CONS, y=CADDmean, group=factor(BP), colour=factor(BP)))+geom_point(alpha=0.5)+xlab("Conservation Score")+ylab("mean CADD score")
ggsave("/net/bipolar/jedidiah/mutation/images/cadd_constr.png")
ggplot(anno2b.plot, aes(x=CONS, y=CADDmean, group=factor(BP), colour=factor(BP)))+geom_point(alpha=0.5)+xlab("Conservation Score")+ylab("mean CADD score")+geom_smooth()
ggsave("/net/bipolar/jedidiah/mutation/images/cadd_constr.png")
head(anno2b.plot)
anno2b.plot<-annot2b.plot %>% filter(n>5)
anno2b.plot <- anno2b.plot %>% filter(n>5)
ggplot(anno2b.plot, aes(x=CONS, y=CADDmean, group=factor(BP), colour=factor(BP)))+geom_point(alpha=0.5)+xlab("Conservation Score")+ylab("mean CADD score")+geom_smooth()
ggsave("/net/bipolar/jedidiah/mutation/images/cadd_constr.png")
ggplot(anno2b.plot, aes(x=factor(, y=CADDmean, group=factor(BP), colour=factor(BP)))+xlab("Conservation Score")+ylab("mean CADD score per gene")+geom_box()
ggplot(anno2b.plot, aes(x=CONS, y=CADDmean, group=factor(BP), colour=factor(BP)))+geom_point(alpha=0.5)+xlab("Conservation Score")+ylab("mean CADD score")+geom_smooth()
ggsave("/net/bipolar/jedidiah/mutation/images/cadd_constr.png")
anno2b.plot.hc <- anno2b.plot %>% filter(CADD>20)
anno2b.plot.hc <- anno2b.plot %>% filter(CADDmean>20)
ggplot(anno2b.plot.hc, aes(x=CONS, y=CADDmean, group=factor(BP), colour=factor(BP)))+geom_point(alpha=0.5)+xlab("Conservation Score")+ylab("mean CADD score")+geom_smooth()
ggsave("/net/bipolar/jedidiah/mutation/images/cadd_constr_hc.png")
dim(anno2b.plot)
1262/2
ggplot(anno2b.plot.hc, aes(x=CONS, y=CADDmean, group=factor(BP), colour=factor(BP)))+geom_point(alpha=0.5)+xlab("Constraint Score")+ylab("mean CADD score")+geom_smooth()
ggsave("/net/bipolar/jedidiah/mutation/images/cadd_constr_hc.png")
ggplot(anno2b.plot, aes(x=CONS, y=CADDmean, group=factor(BP), colour=factor(BP)))+geom_point(alpha=0.5)+xlab("Constraint Score")+ylab("mean CADD score")+geom_smooth()
ggsave("/net/bipolar/jedidiah/mutation/images/cadd_constr.png")
rates1
head(rates)
summary(rates$rel_prop)
summary(rates5$rel_prop)
summary(rates3$rel_prop)
summary(rates$num)
rates %>% filter(Category2=="AT_GC") %>% filter(substr(Sequence,2,5)=="CAAT") %>% summarise(num=sum(num), COUNT=sum(COUNT), prop=num/COUNT)
rates %>% filter(Category2=="AT_GC") %>% filter(substr(Sequence,2,5)%in%c("CAAT", "CTAT")) %>% summarise(num=sum(num), COUNT=sum(COUNT), prop=num/COUNT)
rates %>% filter(Category2=="AT_GC") %>% summarise(num=sum(num), COUNT=sum(COUNT), prop=num/COUNT)
0.024/0.0062
rates %>% filter(Category2=="GC_AT") %>% filter(substr(Sequence,3,4)%in%c("GC")) %>% summarise(num=sum(num), COUNT=sum(COUNT), prop=num/COUNT)
rates %>% filter(Category2=="GC_AT") %>% filter(substr(Sequence,4,5)%in%c("CG")) %>% summarise(num=sum(num), COUNT=sum(COUNT), prop=num/COUNT)
rates %>% filter(Category2=="cpg_GC_AT") %>% filter(substr(Sequence,4,5)%in%c("CG")) %>% summarise(num=sum(num), COUNT=sum(COUNT), prop=num/COUNT)
0.0965/0.013
rates %>% filter(Category2=="AT_TA") %>% filter(substr(Sequence,2,7)%in%c("TTAAAA")) %>% summarise(num=sum(num), COUNT=sum(COUNT), prop=num/COUNT)
0.0098/0.0016
24489*14
ls()
head(plotdat)
dim(plotdat)
head(rates)
head(sigcoefs)
head(sigcoefs) %>% data.frame()
unique(sigcoefs$Category)
rates2<-rates
names(rates2)
names(rates2)[3]<-"C2"
names(rates2)[3]<-"Category2"
names(rates2)[9]<-"C2"
names(rates2)[3]<-"Category"
r3<-merge(sigcoefs, rates2, by=c("Category", "Sequence"))
head(r3)
summary(r3$num)
summary(rates$num)
dim(r3[r3$num>1000,])
r3a<-r3[r3$num>1000,]
dim(r3a)
head(r3a)
r3a<-r3[r3$num>741,]
dim(r3a)
6030/6514
summary(rates$num)
0.026*36087319
36087319-938270
35149049/63566013
0.6*63566013
0.6*63566013-35149049
2990559/(2990559+35149049)
19596845/(19596845+4890329+9719435+2332114+5787895)
head(r7s)
head(chrp)
rsq
rsq<-unlist(lapply(overall_models, function(x) NagelkerkeR2(x)))[c(2,4,6,8,10,12,14,16,18)]
rsq
runDNMLogit<-function(data, group){
if(nrow(data)>1e6){
logmod1a<-glm(OBS~Category.x, data=data, family=binomial())
logmod1b<-glm(OBS~Category, data=data, family=binomial())
logmod3<-glm(OBS~MU_1+resid3, data=data, family=binomial())
logmod5<-glm(OBS~MU_1+resid3+resid5, data=data, family=binomial())
logmod7<-glm(OBS~MU_1+resid3+resid5+resid7, data=data, family=binomial())
logmodL<-glm(OBS~MU_1+resid3+resid5+resid7+residL, data=data, family=binomial())
# logmod3<-glm(OBS~Category.x+MU_3, data=data, family=binomial())
# logmod5<-glm(OBS~MU_1+resid5, data=data, family=binomial())
# logmod5<-glm(OBS~Category.x+MU_5, data=data, family=binomial())
# logmod7<-glm(OBS~MU_1+resid7, data=data, family=binomial())
# logmod7<-glm(OBS~Category.x+MU_S, data=data, family=binomial())
# logmodL<-glm(OBS~MU_1+residLa, data=data, family=binomial())
  logmodS<-glm(OBS~Category+MU_7S, data=data, family=binomial())
  logmodP<-glm(OBS~Category+MU_7P, data=data, family=binomial())
  logmodA<-glm(OBS~Category+MU_A, data=data, family=binomial())
#
# logmod1a<-glm(OBS~Category+MU_1, data=data, family=binomial())
logmod3a<-glm(OBS~Category+MU_3, data=data, family=binomial())
logmod5a<-glm(OBS~Category+MU_5, data=data, family=binomial())
logmod7a<-glm(OBS~Category+MU_S, data=data, family=binomial())
outdat<-list(logmod1a, logmod1b, logmod3, logmod5, logmod7, logmodL,
logmodS, logmodP, logmodA, logmod1a, logmod3a, logmod5a, logmod7a)
} else {
logmod3<-glm(OBS~MU_3, data=data, family=binomial())
logmod5<-glm(OBS~MU_3+resid5, data=data, family=binomial())
logmod7<-glm(OBS~MU_3+resid5+resid7, data=data, family=binomial())
logmodS<-glm(OBS~MU_7S, data=data, family=binomial())
logmodP<-glm(OBS~MU_7P, data=data, family=binomial())
logmodA<-glm(OBS~MU_A, data=data, family=binomial())
logmodL<-glm(OBS~MU_3+resid5+resid7+residL, data=data, family=binomial())
outdat<-list(logmod3, logmod5, logmod7, logmodL,
logmodS, logmodP, logmodA)
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
test1a1b <- lrtest(overall_models[[1]], overall_models[[2]])
test1b3 <- lrtest(overall_models[[2]], overall_models[[3]])
test35 <- lrtest(overall_models[[3]], overall_models[[4]])
test57 <- lrtest(overall_models[[4]], overall_models[[5]])
test7L <- lrtest(overall_models[[5]], overall_models[[6]])
fulllist<-list(test1a1b, test1b3, test35, test57, test7L)
pvals <- lapply(fulllist, function(x) x[[5]][2])
rsq<-unlist(lapply(overall_models, function(x) NagelkerkeR2(x)))[c(2,4,6,8,10,12,14,16,18)]
rsq
runDNMLogit<-function(data, group){
if(nrow(data)>1e6){
logmod1a<-glm(OBS~Category.x, data=data, family=binomial())
logmod1b<-glm(OBS~Category, data=data, family=binomial())
logmod3<-glm(OBS~Category+MU_1+resid3, data=data, family=binomial())
logmod5<-glm(OBS~Category+MU_1+resid3+resid5, data=data, family=binomial())
logmod7<-glm(OBS~Category+MU_1+resid3+resid5+resid7, data=data, family=binomial())
logmodL<-glm(OBS~Category+MU_1+resid3+resid5+resid7+residL, data=data, family=binomial())
# logmod3<-glm(OBS~Category.x+MU_3, data=data, family=binomial())
# logmod5<-glm(OBS~MU_1+resid5, data=data, family=binomial())
# logmod5<-glm(OBS~Category.x+MU_5, data=data, family=binomial())
# logmod7<-glm(OBS~MU_1+resid7, data=data, family=binomial())
# logmod7<-glm(OBS~Category.x+MU_S, data=data, family=binomial())
# logmodL<-glm(OBS~MU_1+residLa, data=data, family=binomial())
  logmodS<-glm(OBS~Category+MU_7S, data=data, family=binomial())
  logmodP<-glm(OBS~Category+MU_7P, data=data, family=binomial())
  logmodA<-glm(OBS~Category+MU_A, data=data, family=binomial())
#
# logmod1a<-glm(OBS~Category+MU_1, data=data, family=binomial())
logmod3a<-glm(OBS~Category+MU_3, data=data, family=binomial())
logmod5a<-glm(OBS~Category+MU_5, data=data, family=binomial())
logmod7a<-glm(OBS~Category+MU_S, data=data, family=binomial())
outdat<-list(logmod1a, logmod1b, logmod3, logmod5, logmod7, logmodL,
logmodS, logmodP, logmodA, logmod1a, logmod3a, logmod5a, logmod7a)
} else {
logmod3<-glm(OBS~MU_3, data=data, family=binomial())
logmod5<-glm(OBS~MU_3+resid5, data=data, family=binomial())
logmod7<-glm(OBS~MU_3+resid5+resid7, data=data, family=binomial())
logmodS<-glm(OBS~MU_7S, data=data, family=binomial())
logmodP<-glm(OBS~MU_7P, data=data, family=binomial())
logmodA<-glm(OBS~MU_A, data=data, family=binomial())
logmodL<-glm(OBS~MU_3+resid5+resid7+residL, data=data, family=binomial())
outdat<-list(logmod3, logmod5, logmod7, logmodL,
logmodS, logmodP, logmodA)
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
test1a1b <- lrtest(overall_models[[1]], overall_models[[2]])
test1b3 <- lrtest(overall_models[[2]], overall_models[[3]])
test35 <- lrtest(overall_models[[3]], overall_models[[4]])
test57 <- lrtest(overall_models[[4]], overall_models[[5]])
test7L <- lrtest(overall_models[[5]], overall_models[[6]])
fulllist<-list(test1a1b, test1b3, test35, test57, test7L)
pvals <- lapply(fulllist, function(x) x[[5]][2])
rsq<-unlist(lapply(overall_models, function(x) NagelkerkeR2(x)))[c(2,4,6,8,10,12,14,16,18)]
rsq
overall_models[[1]]
exp(-4.21)
overall_models[[2]]
overall_models[[3]]
inv.logit(-4.21)
rsq
head(overall_dat)
head(sigcoefs)
max(sigcoefs$Est)
exp(8.68)
Lv7v5v3_full<-ggplot()+
  geom_bar(aes(x=category, y=rsq, fill=mod), stat="identity", position="dodge")+
  # scale_fill_manual("Model", values=cbbPalette[c(1,5,4,6,7,8)])+
scale_fill_manual("Model", values=c(iwhPalette[c(1,2,3,4,5,9)]))+
  theme_bw()+
  theme(
      legend.text=element_text(size=14),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=16),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=14))+
  ylab(expression(paste("Fraction of variance explained (Nagelkerke ", R^2, ")", sep="")))
Lv7v5v3<-ggplot(rsqdat)+
  geom_bar(aes(x=Category, y=rsq, fill=mod), stat="identity", position="dodge")+
  # scale_fill_manual("Model", values=cbbPalette[c(4,6,7,8)])+
scale_fill_manual("Model", values=c(iwhPalette[c(3,4,5,9)]))+
  geom_segment(data=corplot, aes(x=xst, xend=xend, y=yst, yend=yend))+
  geom_segment(data=corplot2, aes(x=xst, xend=xend, y=yst, yend=yend))+
  geom_text(data=corplot2, aes(x=xst+1/9, y=yst+.005, label=pval, angle=90), size=6)+
  scale_y_continuous(limits=c(0,0.1))+
  theme_bw()+
  theme(
      legend.text=element_text(size=14),
      axis.title.x=element_text(size=16),
      axis.title.y=element_text(size=16),
    axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
      axis.text.y=element_text(size=14))+
  xlab("Mutation Type")+
  ylab(expression(paste("Fraction of variance explained (Nagelkerke ", R^2, ")", sep="")))
grid.arrange(grobs=c(Lv7v5v3_full, Lv7v5v3), widths=c(1,2))
require(grid)
grid.arrange(grobs=c(Lv7v5v3_full, Lv7v5v3), widths=c(1,2))
require(gridExtra)
grid.arrange(grobs=c(Lv7v5v3_full, Lv7v5v3), widths=c(1,2))
grid.arrange(Lv7v5v3_full, Lv7v5v3, widths=c(1,2))
corplot
rsq_full_dat<- combineddat %>%
  filter(mod %in% mod[1:6])
Lv7v5v3_full<-ggplot(rsq_full_dat)+
  geom_bar(aes(x=category, y=rsq, fill=mod), stat="identity", position="dodge")+
  # scale_fill_manual("Model", values=cbbPalette[c(1,5,4,6,7,8)])+
scale_fill_manual("Model", values=c(iwhPalette[c(1,2,3,4,5,9)]))+
  theme_bw()+
  theme(
      legend.text=element_text(size=14),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=16),
      axis.text.x=element_blank(),
      axis.text.y=element_text(size=14))+
  ylab(expression(paste("Fraction of variance explained (Nagelkerke ", R^2, ")", sep="")))
grid.arrange(Lv7v5v3_full, Lv7v5v3, widths=c(1,2))
ggsave(paste0(parentdir, "/images/Lv7v5v3_rsq_combined.png"), width=4, height=6)
grid.arrange(grobs=c(Lv7v5v3_full, Lv7v5v3), widths=c(1,2))
head(rsqdat)
Lv7v5v3<-ggplot(rsqdat)+
  geom_bar(aes(x=Category, y=rsq, fill=mod), stat="identity", position="dodge")+
  # scale_fill_manual("Model", values=cbbPalette[c(4,6,7,8)])+
scale_fill_manual("Model", values=c(iwhPalette[c(3,4,5,9)]))+
  geom_segment(data=corplot, aes(x=xst, xend=xend, y=yst, yend=yend))+
  geom_segment(data=corplot2, aes(x=xst, xend=xend, y=yst, yend=yend))+
  geom_text(data=corplot2, aes(x=xst+1/9, y=yst+.005, label=pval, angle=90), size=6)+
  scale_y_continuous(limits=c(0,0.1))+
  theme_bw()+
  theme(
      legend.text=element_text(size=14),
      axis.title.x=element_text(size=16),
      axis.title.y=element_text(size=16),
    axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
      axis.text.y=element_text(size=14))+
  xlab("Mutation Type")+
  ylab(expression(paste("Fraction of variance explained (Nagelkerke ", R^2, ")", sep="")))
grid.arrange(grobs=c(Lv7v5v3_full, Lv7v5v3), widths=c(1,2))

require(gtable)
Lv7v5v3_full<-ggplotGrob(Lv7v5v3_full)
Lv7v5v3<-ggplotGrob(Lv7v5v3)

g<-arrangeGrob(Lv7v5v3_full, Lv7v5v3, nrow=1, widths=c(1,2))
ggsave(paste0(parentdir, "/images/Lv7v5v3_rsq_combined.png"), width=12, height=6, g)

EvC_full<-ggplotGrob(EvC_full)
EvC<-ggplotGrob(EvC)
g<-arrangeGrob(EvC_full, EvC, nrow=1, widths=c(1,2))
ggsave(paste0(parentdir, "/images/EvC_rsq_combined.png"), width=12, height=6, g)

png(paste0(parentdir, "/images/Lv7v5v3_rsq_combined.png"), width=12, height=6)
grid::grid.newpage()
grid::grid.draw(cbind(Lv7v5v3_full, Lv7v5v3))
dev.off()
install_github("baptiste/egg")
require(egg)
require(scales)
require(egg)
ggsave(paste0(parentdir, "/images/Lv7v5v3_rsq_combined.png"), width=12, height=6, g)
ggsave(paste0(parentdir, "/images/Lv7v5v3_rsq_combined.svg"), width=12, height=6, g)
savehistory("/net/bipolar/jedidiah/fixplots.r")
