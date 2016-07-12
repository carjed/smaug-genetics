##############################################################################
# Merge aggregate data with covariate data
##############################################################################
cat("Merging bin-wise data with covariate data...\n")
agg_cov <- merge(agg_5bp_100k, mut_cov, by=c("CHR", "BIN"))

ratiofilter <- 0
if(ratiofilter){
  agg_cov$ratio <- agg_cov$exp/agg_cov$obs
  agg_cov <- filter(agg_cov, ratio<5) %>% mutate(med=(maxn+minn)/2)
}

##############################################################################
# Get relative rates for 1bp and 3bp motifs and check for homogeneity
# among constituent motifs
##############################################################################
cat("Getting relative rates for submotifs...\n")
# Get submotifs
rates5$Seq3 <- substr(rates5$Sequence, adj, nbp-adj+1)
rates5$Seq5 <- substr(rates5$Sequence, adj-1, nbp-1)

# Recalculate rates for submotifs
rates1 <- rates5 %>%
  dplyr::select(Sequence, Category2, COUNT, num) %>%
  group_by(Category2) %>%
  summarise(COUNT1=sum(COUNT), num1=sum(num), rel_prop1=num1/COUNT1)

rates3 <- rates5 %>%
  dplyr::select(Seq3, Category2, COUNT, num) %>%
  group_by(Seq3, Category2) %>%
  summarise(COUNT3=sum(COUNT), num3=sum(num), rel_prop3=num3/COUNT3)

rates5a <- rates5 %>%
  dplyr::select(Seq5, Category2, COUNT, num) %>%
  group_by(Seq5, Category2) %>%
  summarise(COUNT5=sum(COUNT), num5=sum(num), rel_prop5=num5/COUNT5)

# Merge to single data frame
rates_full <- merge(rates5, rates1, by="Category2")
rates_full <- merge(rates_full, rates3, by=c("Seq3", "Category2"))
rates_full <- merge(rates_full, rates5a, by=c("Seq5", "Category2"))

# Identify motifs whose components have most similar relative rates
cat("Checking for homogeneity of similar motifs...\n")
count<-0
results<-data.frame()

for(i in unique(rates5$Category2)){
  aggcat <- rates5[rates5$Category2==i,]
  # names(aggcat)[1] <- "Sequence"
  nbpt <- 7

  tmpresults <- data.frame()
  for(j in unique(substr(aggcat$Sequence,2,nbpt-1))){

    ref <- unique(substr(aggcat$Sequence,adj+1,adj+1))
    dat <- aggcat[substr(aggcat$Sequence,2,nbpt-1)==j,]
    dat$exp <- dat$COUNT*(sum(dat$n)/sum(dat$COUNT))

    test <- chisq.test(dat$n, p=dat$exp/sum(dat$exp))

    datrow <- data.frame(Category2=i, Sequence=j, pval=test$p.value)
    tmpresults <- rbind(tmpresults, datrow)
  }

  tmpresults$Q <- ntile(tmpresults$pval, 10)
  results <- rbind(results, tmpresults)
}

##############################################################################
# Run likelihood analysis on full data
##############################################################################
# ra1a %>% filter(model=="L5")
logit_curves<-0
if(logit_curves){
  source("./logit_curves.r")
}

##############################################################################
# Initialize data of motif counts to use as covariates in negbin model
##############################################################################
cat("Initializing data for negbin models...\n")
dat_5bp_100k$bin$CHR <- as.integer(gsub("chr", "", dat_5bp_100k$bin$CHR))

names(dat_5bp_100k$bin) <- gsub('\\(', '_', names(dat_5bp_100k$bin))
names(dat_5bp_100k$bin) <- gsub('\\)', '_', names(dat_5bp_100k$bin))

atcols <- c(names(dat_5bp_100k$bin)[1:5],
  names(dat_5bp_100k$bin)[which(substr(names(dat_5bp_100k$bin),
    adj+1, adj+1)=="A")])

# For binsAT, etc., must subselect variables

binsAT <- dat_5bp_100k$bin %>%
  select_(.dots = atcols) %>%
	arrange(CHR, BIN)

gcdn <- c("CA", "CC", "CT")
gccols <- c(names(dat_5bp_100k$bin)[1:5],
  names(dat_5bp_100k$bin)[which(substr(names(dat_5bp_100k$bin),
    adj+1, adj+2) %in% gcdn)])

binsGC <- dat_5bp_100k$bin %>%
  select_(.dots = gccols) %>%
	arrange(CHR, BIN)

cpggccols <- c(names(dat_5bp_100k$bin)[1:5],
  names(dat_5bp_100k$bin)[which(substr(names(dat_5bp_100k$bin),
   adj+1, adj+2)=="CG")])

binscpgGC <- dat_5bp_100k$bin %>%
  select_(.dots = cpggccols) %>%
	arrange(CHR, BIN)

##############################################################################
# Run combined negbin model
##############################################################################
# Add expected #singletons per bin using 1bp motifs
# -used in comparing category-specific models with weighted mean method
a3 <- merge(agg_cov, rates1[,c(1,2,4)], by="Category2")
# a3$exp1 <- a3$nmotifs*a3$rel_prop1

logit_gcta<-read.table("/net/bipolar/jedidiah/mutation/output/predicted/binned/GC_TA_binned.txt", header=F, stringsAsFactors=F)
names(logit_gcta)<-c("CHR", "BIN", "GC_TA", "cpg_GC_TA")
lgcta<-logit_gcta %>% gather(Category2, LSUM, GC_TA:cpg_GC_TA)

logit_gccg<-read.table("/net/bipolar/jedidiah/mutation/output/predicted/binned/GC_CG_binned.txt", header=F, stringsAsFactors=F)
names(logit_gccg)<-c("CHR", "BIN", "GC_CG", "cpg_GC_CG")
lgccg<-logit_gccg %>% gather(Category2, LSUM, GC_CG:cpg_GC_CG)

logit_gcat<-read.table("/net/bipolar/jedidiah/mutation/output/predicted/binned/GC_AT_binned.txt", header=F, stringsAsFactors=F)
names(logit_gcat)<-c("CHR", "BIN", "GC_AT", "cpg_GC_AT")
lgcat<-logit_gcat %>% gather(Category2, LSUM, GC_AT:cpg_GC_AT)

logit_atcg<-read.table("/net/bipolar/jedidiah/mutation/output/predicted/binned/AT_CG_binned.txt", header=F, stringsAsFactors=F)
names(logit_atcg)<-c("CHR", "BIN", "LSUM")
latcg<-logit_atcg %>% mutate(Category2="AT_CG") %>% dplyr::select(CHR, BIN, Category2, LSUM)

logit_atgc<-read.table("/net/bipolar/jedidiah/mutation/output/predicted/binned/AT_GC_binned.txt", header=F, stringsAsFactors=F)
names(logit_atgc)<-c("CHR", "BIN", "LSUM")
latgc<-logit_atgc %>% mutate(Category2="AT_GC") %>% dplyr::select(CHR, BIN, Category2, LSUM)

logit_atta<-read.table("/net/bipolar/jedidiah/mutation/output/predicted/binned/AT_TA_binned.txt", header=F, stringsAsFactors=F)
names(logit_atta)<-c("CHR", "BIN", "LSUM")
latta<-logit_atta %>% mutate(Category2="AT_TA") %>% dplyr::select(CHR, BIN, Category2, LSUM)

ldat<-rbind(latcg, latgc, latta, lgcat, lgccg, lgcta)

a3<-merge(a3, ldat, by=c("Category2", "CHR", "BIN"))

# Aggregate across categories for total #singletons per bin
a3a <- a3 %>%
  # group_by_(.dots=lapply(names(a3)[c(2,3,5:17)], as.symbol)) %>%
  group_by_(.dots=lapply(danames, as.symbol)) %>%
  summarise(obs=sum(obs))

# Add motif counts per bin
a3a1 <- merge(a3a, dat_5bp_100k$bin, by=c("CHR", "BIN", "prop_GC"))

bases <- c("A", "C", "G", "T")
nts <- c("A", "C")

# Uses poisson regression instead of negbin, since negbin fails to converge
# with default parameters of glm.nb()
overall<-0
if(overall){
  cat("Running overall models...\n")

  cat("Formatting data...\n")
  a3a2 <- getSubMotifs(a3a1, nts, bases)

  # Create model formulas
  # Must fix to use filtered motifs
  motif_mod_form <- as.formula(paste("obs~",
  	paste(names(a3a2)[(ncol(a3a2)-511):ncol(a3a2)], collapse="+")))

  feat_mod_form <- as.formula(paste("obs~",
    paste(c(covnames, "prop_GC"), collapse="+")))

  full_mod_form <- as.formula(paste("obs~",
  	paste(names(a3a2)[c(3:15, (ncol(a3a2)-511):ncol(a3a2))], collapse="+")))

  mut_lm_m_all <- glm(motif_mod_form, data=a3a2, family="poisson")
  mut_lm_fe_all <- glm(feat_mod_form, data=a3a2, family="poisson")
  mut_lm_f_all <- glm(full_mod_form, data=a3a2, family="poisson")

  # Calculate McFadden's pseudo R-squared (unadjusted)
  mrsq <- 1-mut_lm_m_all$deviance/mut_lm_m_all$null.deviance
  fersq <- 1-mut_lm_fe_all$deviance/mut_lm_fe_all$null.deviance
  frsq <- 1-mut_lm_f_all$deviance/mut_lm_f_all$null.deviance
  maic <- AIC(mut_lm_m_all)
  feaic <- AIC(mut_lm_fe_all)
  faic <- AIC(mut_lm_f_all)

  cat("Adjusted R-squared of combined model (motifs): ", mrsq, "\n")
  cat("AIC of combined model (motifs): ", maic, "\n")

  cat("Adjusted R-squared of combined model (features): ", fersq, "\n")
  cat("AIC of combined model (features): ", feaic, "\n")

  cat("Adjusted R-squared of combined model (full): ", frsq, "\n")
  cat("AIC of combined model (full): ", faic, "\n")
}

##############################################################################
# Model each category independently
##############################################################################
ptm <- proc.time()
cat("Running per-category models...\n")
source("./negbin_run.r")
tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

# Specify category as factor
compare.all$Category2 <- factor(compare.all$Category2)
# levels(compare.all$Category2) <- c(
	# "AT>CG", "AT>GC", "AT>TA",
	# "(CpG)GC>AT", "(CpG)GC>CG", "(CpG)GC>TA",
	# "GC>AT", "GC>CG", "GC>TA")

overallcor<-compare.all %>%
  group_by(res) %>%
  summarise(cor=cor(exp,obs), corsq=cor^2)

##############################################################################
# Plot barcharts comparing obs/exp correlation for different models
##############################################################################
mod.corr <- compare.all %>%
	# filter(exp > 100) %>%
	group_by(Category2, res) %>%
	summarise(num=length(exp),
		cor=cor(exp, obs, method="pearson"),
		cor.p=cor.test(exp, obs, method="pearson")$p.value)
mod.corr$SE <- corSE(mod.corr$cor, mod.corr$num)
limits <- aes(ymax = mod.corr$cor + mod.corr$SE,
	ymin=mod.corr$cor - mod.corr$SE)
dodge <- position_dodge(width=0.9)

mc2<-mod.corr %>%
  filter(res %in% c("full", "logit", "marginal7", "motifs5")) %>%
  mutate(gp=ifelse(res %in% c("full", "logit"), "Sequence+Features", "Sequence Only"))
mc2$res<-as.factor(mc2$res)
levels(mc2$res)<-c("K-mers+Features", "Positional+Features", "Positional", "K-mers")
mc2$res<-factor(mc2$res,levels(mc2$res)[c(3,2,4,1)])

customPalette <- brewer.pal(4, "RdBu")[c(3,4,2,1)]
# "#92C5DE" "#F4A582" "#0571B0" "#CA0020"

ggplot(mc2, aes(x=Category2, y=cor, colour=res, fill=res, group=res))+
	geom_bar(stat="identity", position=dodge)+
  # geom_point(position=dodge)+
  # scale_colour_brewer("Predictor",palette="RdBu")+
  # scale_fill_brewer("Predictor", palette="RdBu")+
  scale_colour_manual("Model", values=customPalette)+
  scale_fill_manual("Model", values=customPalette)+
	xlab("Category")+
	ylab("Correlation")+
	# geom_errorbar(limits, position=dodge, width=0.25)+
  # facet_wrap(~gp, nrow=1)+
	theme_bw()+
	theme(legend.title = element_text(size=18),
		legend.text = element_text(size=14),
    legend.position="bottom",
		# axis.title.x = element_text(size=20),
    axis.title.x = element_blank(),
		axis.title.y = element_text(size=20),
		axis.text.y = element_text(size=16, angle=90, hjust=0.5),
    axis.ticks.x = element_blank(),
    # axis.text.x = element_blank())
		axis.text.x = element_text(size=16, angle = 45,  vjust=1, hjust=1.01))

modelbar <- paste0(parentdir, "/images/gw_5bp_vs_mod_3.png")
ggsave(modelbar, width=8, height=6)

# Get AIC for each model/category
compare.aic <- compare.aic %>% spread(Category2, AIC)

mc2b<-mod.corr %>%
  filter(res %in% c("full", "logit", "marginal7", "motifs5")) %>%
  ungroup() %>%
  dplyr::select(Category2, num, res, cor) %>%
  spread(res, cor) %>%
  mutate(corpf=r.test(num, full, logit)$p,
    corps=r.test(num, marginal7, motifs5)$p)


comp.mods<-compare.all %>%
  filter(res %in% c("full", "logit", "marginal7", "motifs5")) %>%
  mutate(gp=ifelse(res %in% c("full", "logit"), "Sequence+Features", "Sequence Only"))
comp.mods$res<-as.factor(comp.mods$res)
levels(comp.mods$res)<-c("K-mers+Features", "Positional+Features", "Positional", "K-mers")
comp.mods$res<-factor(comp.mods$res,levels(comp.mods$res)[c(3,2,4,1)])

ggplot(comp.mods, aes(x=obs, y=exp, group=res, colour=res))+
  geom_point(alpha=0.3)+
  scale_colour_manual("Model", values=customPalette)+
  xlab("Observed Singletons")+
  ylab("Predicted Singletons")+
  geom_smooth(method=lm, se=F, colour="black")+
  facet_wrap(~res, scales="free")+
  theme_bw()+
  theme(legend.position="none")
ggsave("/net/bipolar/jedidiah/mutation/images/1mb_scatter.png", height=4, width=5.5)


##############################################################################
# Plot per-chromosome variation
##############################################################################
cat("Plotting per-chromosome variation...\n")
require(ggbio)
data(hg19IdeogramCyto, package = "biovizBase")
ic2<-as.data.frame(hg19IdeogramCyto)

# Read in chromosome lengths for specifying plot range
chrlen<-read.table("/net/bipolar/jedidiah/mutation/reference_data/hg19.genome", header=T, stringsAsFactors=F)

plotlist<-list()
for(i in 1:22){
  # Specify chromosome
  chrname<-paste0("chr",i)

  #Subset and update data
  d2<-compare.all %>%
    filter(CHR==i, res=="full") %>%
    mutate(ratio=exp/obs, diff=exp-obs) %>%
    filter(abs(diff)<1200)
  d2$BIN<-d2$BIN*1000000

  # Get chromosome band
  ic2chr<-ic2[ic2$seqnames==chrname,]
  ic2chr$mean<-(ic2chr$start+ic2chr$end)/2
  oetest<-t.test(d2$exp, d2$obs, conf.level=0.99)
  testmin<-oetest$conf.int[[1]]
  testmax<-oetest$conf.int[[2]]
  ic2a<-ic2chr %>%
    group_by(seqnames) %>%
    summarise(min=min(start), max=max(end))

  # Plot scatterplot & loess curve, with banding overlay
  p2<-ggplot()+
    geom_rect(data=ic2chr,
      aes(xmin=start, xmax=end,
        # ymin=min(d2$ratio)-0.3, ymax=min(d2$ratio)-0.1, fill=gieStain),
        ymin=testmin, ymax=testmax, fill=gieStain),
      alpha=0.6)+
    geom_rect(data=ic2a,
      aes(xmin=min, xmax=max,
        ymin=testmin, ymax=testmax),
      colour="black", alpha=0)+
    scale_fill_manual(values=getOption("biovizBase")$cytobandColor)+
    geom_point(data=d2[d2$ratio<2,],
        aes(x=BIN, y=diff, colour=Category2, group=Category2), alpha=0.3)+
    geom_smooth(data=d2[d2$ratio<2,],
        aes(x=BIN, y=diff, colour=Category2, group=Category2), span=0.2, se=FALSE)+
    scale_colour_manual("Category", values=myPaletteCat(9))+
    scale_x_continuous(limits=c(0,chrlen[chrlen$chrom==chrname,]$size))+
    # scale_y_continuous(limits=c(0.25,2))+
    # geom_hline(yintercept=0)+
    # geom_hline(yintercept=testmin, linetype="dashed")+
    # geom_hline(yintercept=testmax, linetype="dashed")+
    # geom_text(data=ic2chr,
    #   aes(x=mean, y=min(d2$ratio)-0.3, label=name),
    #   hjust=0, angle=90, size=3)+
    # facet_wrap(~Category2, scales="free", ncol=1)+
    ylab("Predicted-Observed")+
    xlab(NULL)+
    theme_bw()+
    theme(legend.justification = c(0, 1),
      legend.position = c(0, 1),
      legend.title=element_blank(),
      legend.text=element_text(size=14),
      axis.title.y=element_text(size=16),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
    guides(fill=FALSE, colour = guide_legend(nrow = 3))
  plotlist[[i]]<-p2
  p2
  ratiofile<-paste0("/net/bipolar/jedidiah/mutation/images/chr",i,"_ratio.png")
  ggsave(ratiofile, width=10, height=5)
}

pdf("/net/bipolar/jedidiah/mutation/images/allchr.pdf", onefile=T,
  width=8, height=88, paper="letter")
multiplot(plotlist=plotlist[1:4], cols=1)
multiplot(plotlist=plotlist[5:8], cols=1)
multiplot(plotlist=plotlist[9:12], cols=1)
multiplot(plotlist=plotlist[13:16], cols=1)
multiplot(plotlist=plotlist[17:20], cols=1)
multiplot(plotlist=plotlist[21:22], cols=1)
dev.off()
