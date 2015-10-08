if(negbin_model){
cat("Initializing negbin regression model...\n")
ptm <- proc.time()
source("negbin_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Finished in", tottime, "seconds \n")
}
if(negbin_model){
cat("Initializing negbin regression model...\n")
ptm <- proc.time()
source("negbin_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Finished in", tottime, "seconds \n")
}
gc()
head(compare.all)
c2<-summarise(compare.all, me=mean(exp), mo=mean(obs))
head(c2)
c2
data.frame(c2)
head(summagg2)
head(s2)
head(agg_cov)
if(run_agg){
ptm <- proc.time()
cat("Aggregating data...\n")
source("agg_dat.r")
aggV <- aggData(dat_5bp_100k, 2) #<-modify the adj value for 3bp data
agg_5bp_100k <- aggV$oe
rates1 <- aggV$agg
ratefile <- paste0(parentdir, "/output/", nbp, "bp_", bink, "k_rates.txt")
write.table(rates1, ratefile, col.names=T, row.names=F, quote=F, sep="\t")
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
}
##############################################################################
# Build file of covariates and run PCA
##############################################################################
ptm <- proc.time()
mutcov2file<-paste0(parentdir, "/output/logmod_data/", bink, "kb_mut_cov2.txt")
if(!file.exists(mutcov2file)){
cat("Building covariate data...\n")
source("get_covs.r")
} else {
cat("Reading existing covariate datafile:", mutcov2file, "...\n")
mut_cov<-read.table(mutcov2file, header=F, stringsAsFactors=F)
}
if(pcs==1){
names(mut_cov)<-c("CHR", "BIN", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13")
} else {
names(mut_cov)<-c("CHR", "BIN", "H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3", "H3K27ac", "H3K27me3", "H3K36me3", "CPGI", "EXON", "TIME", "RATE", "prop_GC", "LAMIN")
}
danames<-names(mut_cov)
covnames<-danames[-c(1:2, 14)]
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
##############################################################################
# Run negative binomial regression models and plot
##############################################################################
if(negbin_model){
cat("Initializing negbin regression model...\n")
ptm <- proc.time()
source("negbin_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Finished in", tottime, "seconds \n")
}
ggplot(agg_cov, aes(x=prop_GC, y=motifvar))+
geom_point(alpha=0.5)+
scale_colour_manual(values=myPaletteCat(8)[3])+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
ggplot(agg_cov, aes(x=prop_GC, y=motifvar))+
geom_point(alpha=0.5)+
scale_colour_manual(values=myPaletteCat(8)[3:4])+
facet_wrap(~Category2, scales="free")+
theme_bw()
ggplot(agg_cov, aes(x=prop_GC, y=motifvar, colour=myPaletteCat(8)[3]))+
geom_point(alpha=0.5)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
ggsave(gcfile)
ggplot(agg_cov, aes(x=prop_GC, y=motifvar))+
geom_point(alpha=0.5, colour=myPaletteCat(8)[3])+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
ggsave(gcfile)
ggplot(agg_cov, aes(x=prop_GC, y=motifvar))+
geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
ggsave(gcfile)
head(agg_cov)
agg_cov_cor<-agg_cov %>% group_by(Category2) %>% summarise(cor=cor(prop_GC, motifvar))
agg_cov_cor
agg_cov_cor<-agg_cov %>% group_by(Category2) %>% summarise(cor=cor(prop_GC, motifvar, method="spearman"))
agg_cov_cor
ggplot(agg_cov, aes(x=prop_GC, y=nmotifs))+
geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
ggsave(gcfile)
ggplot(agg_cov, aes(x=prop_GC, y=motifvar))+
geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
ggsave(gcfile)
head(agg_5bp_100k)
max(agg_5bp_100k$obs)
if(run_agg){
ptm <- proc.time()
cat("Aggregating data...\n")
source("agg_dat.r")
aggV <- aggData(dat_5bp_100k, 2) #<-modify the adj value for 3bp data
agg_5bp_100k <- aggV$oe
rates1 <- aggV$agg
ratefile <- paste0(parentdir, "/output/", nbp, "bp_", bink, "k_rates.txt")
write.table(rates1, ratefile, col.names=T, row.names=F, quote=F, sep="\t")
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
}
##############################################################################
# Build file of covariates and run PCA
##############################################################################
ptm <- proc.time()
mutcov2file<-paste0(parentdir, "/output/logmod_data/", bink, "kb_mut_cov2.txt")
if(!file.exists(mutcov2file)){
cat("Building covariate data...\n")
source("get_covs.r")
} else {
cat("Reading existing covariate datafile:", mutcov2file, "...\n")
mut_cov<-read.table(mutcov2file, header=F, stringsAsFactors=F)
}
if(pcs==1){
names(mut_cov)<-c("CHR", "BIN", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13")
} else {
names(mut_cov)<-c("CHR", "BIN", "H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3", "H3K27ac", "H3K27me3", "H3K36me3", "CPGI", "EXON", "TIME", "RATE", "prop_GC", "LAMIN")
}
danames<-names(mut_cov)
covnames<-danames[-c(1:2, 14)]
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
##############################################################################
# Run negative binomial regression models and plot
##############################################################################
if(negbin_model){
cat("Initializing negbin regression model...\n")
ptm <- proc.time()
source("negbin_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Finished in", tottime, "seconds \n")
}
ggplot(agg_cov, aes(x=prop_GC, y=med))+
geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
geom_errorbar(ymax=maxn, ymin=minn)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
head(aggcov)
head(agg_5bp_100k)
gg_cov <- merge(agg_5bp_100k, mut_cov, by=c("CHR", "BIN"))
agg_cov$ratio <- agg_cov$exp/agg_cov$obs
agg_cov <- filter(agg_cov, ratio<5)
agg_cov$med <- (agg_cov$maxn, agg_cov$minn)
ggplot(agg_cov, aes(x=prop_GC, y=med))+
geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
geom_errorbar(ymax=maxn, ymin=minn)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
agg_cov <- merge(agg_5bp_100k, mut_cov, by=c("CHR", "BIN"))
agg_cov$ratio <- agg_cov$exp/agg_cov$obs
agg_cov <- filter(agg_cov, ratio<5)
agg_cov$med <- (agg_cov$maxn, agg_cov$minn)
ggplot(agg_cov, aes(x=prop_GC, y=med))+
geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
geom_errorbar(ymax=maxn, ymin=minn)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
agg_cov <- merge(agg_5bp_100k, mut_cov, by=c("CHR", "BIN"))
agg_cov$ratio <- agg_cov$exp/agg_cov$obs
agg_cov <- filter(agg_cov, ratio<5)
agg_cov$med <- mean(agg_cov$maxn, agg_cov$minn)
ggplot(agg_cov, aes(x=prop_GC, y=med))+
geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
geom_errorbar(ymax=maxn, ymin=minn)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
ggsave(gcfile)
head(agg_cov)
agg_cov <- merge(agg_5bp_100k, mut_cov, by=c("CHR", "BIN"))
agg_cov$ratio <- agg_cov$exp/agg_cov$obs
agg_cov <- filter(agg_cov, ratio<5) %>% mutate(med=mean(maxn, minn))
agg_cov <- merge(agg_5bp_100k, mut_cov, by=c("CHR", "BIN"))
agg_cov$ratio <- agg_cov$exp/agg_cov$obs
agg_cov <- filter(agg_cov, ratio<5) %>% mutate(med=(maxn+minn)/2)
ggplot(agg_cov, aes(x=prop_GC, y=med))+
geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
geom_errorbar(ymax=maxn, ymin=minn)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
head(agg_cov)
ggplot(agg_cov, aes(x=prop_GC, y=med))+
geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
geom_errorbar(aes(ymax=maxn, ymin=minn))+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
ggsave(gcfile)
quantile(agg_cov$med)
quantile(agg_cov$med)[2]
ggplot(agg_cov, aes(x=motifvar, y=prop_GC))+
geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
# geom_errorbar(aes(ymax=maxn, ymin=minn))+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
ggsave(gcfile)
ggplot(agg_cov, aes(x=prop_GC, y=motifvar))+
geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
# geom_errorbar(aes(ymax=maxn, ymin=minn))+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
ggsave(gcfile)
head(summagg2)
head(agg_5bp_100k)
agg_cov2 <- merge(summagg2, mut_cov, by=c("CHR", "BIN"))
ggplot(agg_cov2, aes(x=obs, y=prop_GC, colour=Sequence))+
geom_point(alpha=0.3)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcobsfile <- paste0(parentdir, "/images/gc_vs_obs.png")
head(agg_cov2)
head(mut_cov)
head(summagg2)
ggplot(agg_cov, aes(x=prop_GC.x, y=motifvar))+
geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
# geom_errorbar(aes(ymax=maxn, ymin=minn))+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
ggsave(gcfile)
ggplot(agg_cov, aes(x=prop_GC.y, y=motifvar))+
geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
# geom_errorbar(aes(ymax=maxn, ymin=minn))+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
ggsave(gcfile)
head(agg_cov)
head(agg_cov2)
ggplot(agg_cov2, aes(x=obs, y=prop_GC.x, colour=Sequence))+
geom_point(alpha=0.3)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcobsfile <- paste0(parentdir, "/images/gc_vs_obs.png")
ggsave(gcobsfile)
ggplot(agg_cov2, aes(x=prop_GC.x, y=obs, colour=Sequence))+
geom_point(alpha=0.3)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcobsfile <- paste0(parentdir, "/images/gc_vs_obs.png")
ggsave(gcobsfile)
ggplot(agg_cov2, aes(x=prop_GC.x, y=exp, colour=Sequence))+
geom_point(alpha=0.3)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcobsfile <- paste0(parentdir, "/images/gc_vs_obs.png")
ggsave(gcobsfile)
head(summagg2)
agg_cov2 <- summagg2 %>%
group_by(Category2, Sequence) %>%
summarise(cor=cor(obs, prop_GC, use="complete.cases"))
agg_cov2 <- summagg2 %>%
group_by(Category2, Sequence) %>%
summarise(cor=cor(obs, prop_GC, use="complete.obs"))
head(agg_cov2)
max(agg_cov2$cor)
head(summagg2)
agg_cov2 <- summagg2 %>%
group_by(Category2, Sequence) %>%
summarise(cor=cor(obs, prop_GC, use="complete.obs"), nmotifs=sum(Count))
ggplot(agg_cov2, aes(x=nmotifs, y=cor))+
geom_point(alpha=0.3)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
head(agg_cov2)
dim(agg_cov2)
agg_cov2 %>% summarise(cor=cor(cor, nmotifs))
agg_cov2 %>% summarise(cor=cor(cor, nmotifs, method="spearman"))
head(summagg2)
head(agg_cov2)
head(summagg2)
agg_cov2 <- summagg2 %>%
group_by(Category2, Sequence) %>%
summarise(cor=cor(obs, prop_GC, use="complete.obs"), nmotifs=sum(Count), rel_prop=mean(rel_prop, na.rm=T))
ggplot(agg_cov2, aes(x=nmotifs, y=cor, size=rel_prop))+
geom_point(alpha=0.3)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
ggplot(agg_cov2, aes(x=nmotifs, y=cor, size=-log(rel_prop)))+
geom_point(alpha=0.3)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
ggplot(filter(agg_cov2, Category2=="GC_TA"), aes(x=nmotifs, y=cor, size=rel_prop))+
geom_point(alpha=0.3)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
ggplot(filter(agg_cov2, Category2=="AT_GC"), aes(x=nmotifs, y=cor, size=rel_prop))+
geom_point(alpha=0.3)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
ggplot(filter(agg_cov2, Category2=="GC>AT"), aes(x=nmotifs, y=cor, size=rel_prop))+
geom_point(alpha=0.3)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
ggplot(filter(agg_cov2, Category2=="GC_AT"), aes(x=nmotifs, y=cor, size=rel_prop))+
geom_point(alpha=0.3)+
facet_wrap(~Category2, scales="free")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
ggplot(filter(agg_cov2, Category2=="GC_TA"), aes(x=nmotifs, y=cor, size=rel_prop))+
geom_point(alpha=0.3)+
facet_wrap(~Category2, scales="free")+
xlab("# mutable sites")+
ylab("%GC~singleton correlation")+
theme_bw()+
theme(legend.title=element_text("Relative Mutation Rate"))
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)ggplot(filter(agg_cov2, Category2=="GC_TA"), aes(x=nmotifs, y=cor, size=rel_prop))+
geom_point(alpha=0.3)+
scale_colour_brewer(name="Relative Mutation Rate")+
facet_wrap(~Category2, scales="free")+
xlab("# mutable sites")+
ylab("%GC~singleton correlation")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggplot(filter(agg_cov2, Category2=="GC_TA"), aes(x=nmotifs, y=cor, size=rel_prop))+
geom_point(alpha=0.3)+
scale_colour_brewer(name="Relative Mutation Rate")+
facet_wrap(~Category2, scales="free")+
xlab("# mutable sites")+
ylab("%GC~singleton correlation")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
ggplot(filter(agg_cov2, Category2=="GC_TA"), aes(x=nmotifs, y=cor, size=rel_prop))+
geom_point(alpha=0.3)+
scale_size_continuous(name="Relative Mutation Rate")+
facet_wrap(~Category2, scales="free")+
xlab("# mutable sites")+
ylab("%GC~singleton correlation")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
ggplot(filter(agg_cov2, Category2=="AT_CG"), aes(x=nmotifs, y=cor, size=rel_prop))+
geom_point(alpha=0.3)+
scale_size_continuous(name="Relative Mutation Rate")+
facet_wrap(~Category2, scales="free")+
xlab("# mutable sites")+
ylab("%GC~singleton correlation")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
ggplot(filter(agg_cov2, Category2=="GC_TA"), aes(x=nmotifs, y=cor, size=rel_prop))+
geom_point(alpha=0.3)+
scale_size_continuous(name="Relative Mutation Rate")+
facet_wrap(~Category2, scales="free")+
xlab("# mutable sites")+
ylab("%GC~singleton correlation")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
agg_cov2 <- summagg2 %>%
group_by(Category2, Sequence) %>%
summarise(corgc=cor(obs, prop_GC, use="complete.obs"), 
cornm=cor(obs, nmotifs, use="complete.obs"),
nmotifs=sum(Count), rel_prop=mean(rel_prop, na.rm=T))
# agg_cov2 %>% summarise(cor=cor(cor, nmotifs, method="spearman"))
ggplot(aggcov2, aes(x=corgc, y=cornm, size=rel_prop))+
geom_point(alpha=0.3)+
scale_size_continuous(name="Relative Mutation Rate")+
facet_wrap(~Category2, scales="free")+
xlab("# mutable sites~singleton correlation")+
ylab("%GC~singleton correlation")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
ggplot(agg_cov2, aes(x=corgc, y=cornm, size=rel_prop))+
geom_point(alpha=0.3)+
scale_size_continuous(name="Relative Mutation Rate")+
facet_wrap(~Category2, scales="free")+
xlab("# mutable sites~singleton correlation")+
ylab("%GC~singleton correlation")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
head(agg_cov2)
agg_cov2 <- summagg2 %>%
group_by(Category2, Sequence) %>%
summarise(corgc=cor(obs, prop_GC, use="complete.obs"), 
cornm=cor(obs, nmotifs, use="complete.obs"),
nmotifs=sum(Count), rel_prop=mean(rel_prop, na.rm=T))
head(summagg2)
agg_cov2 <- summagg2[complete.cases(summagg2),] %>%
group_by(Category2, Sequence) %>%
summarise(corgc=cor(obs, prop_GC, use="complete.obs"), 
cornm=cor(obs, nmotifs, use="complete.obs"),
nmotifs=sum(Count), rel_prop=mean(rel_prop, na.rm=T))
agg_cov2 <- summagg2 %>%
group_by(Category2, Sequence) %>%
summarise(
cornm=cor(obs, nmotifs, use="complete.obs"),
nmotifs=sum(Count), rel_prop=mean(rel_prop, na.rm=T))
head(summagg2)
agg_cov2 <- summagg2 %>%
group_by(Category2, Sequence) %>%
summarise(corgc=cor(obs, prop_GC, use="complete.obs"), 
cornm=cor(obs, Count, use="complete.obs"),
nmotifs=sum(Count), rel_prop=mean(rel_prop, na.rm=T))
ggplot(agg_cov2, aes(x=corgc, y=cornm, size=rel_prop))+
geom_point(alpha=0.3)+
scale_size_continuous(name="Relative Mutation Rate")+
facet_wrap(~Category2, scales="free")+
xlab("# mutable sites~singleton correlation")+
ylab("%GC~singleton correlation")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
ggplot(agg_cov2, aes(x=corgc, y=cornm, colour=rel_prop))+
geom_point(alpha=0.3)+
scale_colour_continuous(name="Relative Mutation Rate")+
facet_wrap(~Category2, scales="free")+
xlab("# mutable sites~singleton correlation")+
ylab("%GC~singleton correlation")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
agg_cov2$mgc <- countChars("G", agg_cov2$Sequence)+countChars("C", agg_cov2$Sequence)
countChars <- function(char, s) {
    s2 <- gsub(char,"",s)
    return (nchar(s) - nchar(s2))
}
agg_cov2$mgc <- countChars("G", agg_cov2$Sequence)+countChars("C", agg_cov2$Sequence)
agg_cov2$mgc <- countChars(G, agg_cov2$Sequence)+countChars(C, agg_cov2$Sequence)
gcCount2 <-  function(line, st, sp){
  length(gregexpr('[GCgc]', substr(line, st, sp))[[1]])
}
agg_cov2$mgc <- gcCount2(agg_cov2$Sequence, 1, 12)
head(agg_cov2)
tail(agg_cov2)
agg_cov2$mgc<-gregexpr('[GC]', aggcov2$Sequence)
agg_cov2$mgc<-gregexpr('[GC]', agg_cov2$Sequence)
tail(agg_cov2)
agg_cov2$mgc<-length(gregexpr('[GC]', agg_cov2$Sequence))
tail(agg_cov2)
agg_cov2$mgc<-length(gregexpr('[GC]', agg_cov2$Sequence))[[1]]
tail(agg_cov2)
agg_cov2$mgc<-length(gregexpr('[GC]', agg_cov2$Sequence)[[1]])
tail(agg_cov2)
agg_cov2$mgc <- nchar(gsub("[AT]", "", agg_cov2$Sequence))
tail(agg_cov2)
agg_cov2$mgc <- nchar(gsub("[AT]", "", agg_cov2$Sequence))-2
head(agg_cov2)
ggplot(agg_cov2, aes(x=corgc, y=cornm, colour=mgc))+
geom_point(alpha=0.3)+
scale_colour_continuous(name="Relative Mutation Rate")+
facet_wrap(~Category2, scales="free")+
xlab("# mutable sites~singleton correlation")+
ylab("%GC~singleton correlation")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
ggplot(agg_cov2, aes(x=corgc, y=cornm, colour=mgc, size=nmotifs))+
geom_point(alpha=0.3)+
scale_colour_continuous(name="Relative Mutation Rate")+
facet_wrap(~Category2, scales="free")+
xlab("# mutable sites~singleton correlation")+
ylab("%GC~singleton correlation")+
theme_bw()
gccorfile <- paste0(parentdir, "/images/nmotifs_vs_gccor.png")
ggsave(gccorfile)
ggplot(agg_cov, aes(x=prop_GC, y=motifvar))+
geom_point(alpha=0.3, colour=myPaletteCat(8)[3])+
# geom_errorbar(aes(ymax=maxn, ymin=minn))+
facet_wrap(~Category2, scales="free")+
theme_bw()
gcfile <- paste0(parentdir, "/images/gc_vs_var.png")
ggsave(gcfile)
savehistory("/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/cor_plots.r")
