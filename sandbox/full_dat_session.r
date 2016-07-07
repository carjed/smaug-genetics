dim(cad)
cad
head(data.frame(compare.all))
head(data.frame(compare.all)[compare.all$obs>10,])
c.gp<-group_by(data.frame(compare.all)[compare.all$obs>50,], Category2, res)
ca<-summarise(c.gp, mean=mean(err))
cad<-dcast(data.frame(ca), Category2~res)
cad
dim(data.frame(compare.all)[compare.all$obs>50,])
dim(compare.all)
c.gp<-group_by(data.frame(compare.all)[compare.all$obs>10,], Category2, res)
ca<-summarise(c.gp, mean=mean(err))
cad<-dcast(data.frame(ca), Category2~res)
cad
ca.gp <- group_by(compare.all, Category2, res)
ca<-summarise(ca.gp, mean=mean(err), meanobs=mean(obs))
cad<-dcast(data.frame(ca), Category2~res)
cad
cad<-dcast(data.frame(ca), Category2~res, value.var=mean)
cad<-dcast(data.frame(ca), Category2~res, value.var="mean")
cad
c.gp<-group_by(data.frame(compare.all)[compare.all$obs>10,], Category2, res)
ca.gp <- group_by(compare.all, Category2, res)
ca<-summarise(c.gp, mean=mean(err), meanobs=mean(obs))
cad<-dcast(data.frame(ca), Category2~res, value.var=mean)
cad<-dcast(data.frame(ca), Category2~res, value.var="mean")
cad
mean(cad)
mean(cad$features)
mean(cad[c(1:3,7:9),]$features)
mean(cad[c(1:3,7:9),]$5bp)
names(cad)<-c("Category2", "features", "motifs", "both")
mean(cad[c(1:3,7:9),]$motifs)
mean(cad[c(1:3,7:9),]$both)
q()
head(ca)
ggplot(ca, aes(x=Category2, y=mean, group=res, fill=res))+geom_bar(stat=identity, position=dodge)
ggplot(ca, aes(x=Category2, y=mean, group=res, fill=res))+geom_bar(stat="identity", position="dodge")
ggsave("/net/bipolar/jedidiah/mutation/images/mod_comp.png")
limits <- aes(ymax = ca$mean + ca$sd, ymin=ca$mean-ca$sd)
dodge <- position_dodge(width=0.9)
ggplot(ca, aes(x=Category2, y=mean, group=res, fill=res))+
geom_bar(stat="identity", position="dodge")+
geom_errorbar(limits, position=dodge)
ggplot(ca, aes(x=Category2, y=mean, fill=res))+
geom_bar(stat="identity", position="dodge")+
geom_errorbar(limits, position=dodge)
limits
head(ca)
ca<-summarise(c.gp, mean=mean(err), sd=sd(err), meanobs=mean(obs))
cad<-dcast(data.frame(ca), Category2~res, value.var="mean")
limits <- aes(ymax = ca$mean + ca$sd, ymin=ca$mean-ca$sd)
dodge <- position_dodge(width=0.9)
ggplot(ca, aes(x=Category2, y=mean, fill=res))+
geom_bar(stat="identity", position="dodge")+
geom_errorbar(limits, position=dodge)
ggsave("/net/bipolar/jedidiah/mutation/images/mod_comp.png")
cad
head(agg_5bp_100k)
mean(agg_5bp_100k$diff)
caout<-paste0(parentdir, "/output/", nbp, "bp_err.txt")
write.table(ca, caout, col.names=T, row.names=F, quote=F, sep="\t")
ca<-summarise(c.gp, mean=mean(err), sd=std(err), meanobs=mean(obs))
caout<-paste0(parentdir, "/output/", nbp, "bp_err.txt")
write.table(ca, caout, col.names=T, row.names=F, quote=F, sep="\t")
d1<-read.table("/net/bipolar/jedidiah/mutation/output/3bp_err.txt", header=T)
d2<-read.table("/net/bipolar/jedidiah/mutation/output/5bp_err.txt", header=T)
d1a<-d1[d1$res=="3bp",]
d3<-rbind(d2, d1a)
limits <- aes(ymax = d3$mean + d3$sd, ymin=d3$mean-d3$sd)
dodge <- position_dodge(width=0.9)
ggplot(d3, aes(x=Category2, y=mean, fill=res))+
geom_bar(stat="identity", position="dodge")+
geom_errorbar(limits, position=dodge)
ggsave("/net/bipolar/jedidiah/mutation/images/mod_comp.png")
d3
head(agg_cov)
sum(agg_cov$obs)
dim(dat_5bp_100k$summ)
dim(ca)
cad
head(d1)
head(d1a)
source("compare_motif_err.r")
head(aggcat)
if(run_agg==1){
ptm <- proc.time()
cat("Aggregating data...\n")
source("agg_dat.r")
aggV <- aggData(dat_5bp_100k, 2) #<-modify the adj value for 3bp data
agg_5bp_100k <- aggV$oe
rates1 <- aggV$agg
write.table(rates1, "/net/bipolar/jedidiah/mutation/output/5bp_100k_rates.txt", col.names=T, row.names=F, quote=F, sep="\t")
agg_5bp_100k <- agg_5bp_100k[order(agg_5bp_100k$BIN),]
agg_5bp_100k$diff <- agg_5bp_100k$obs-agg_5bp_100k$exp
agg_5bp_100k$Category2 <- as.character(agg_5bp_100k$Category2)
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
}
ref
if(run_agg==1){
ptm <- proc.time()
cat("Aggregating data...\n")
source("agg_dat.r")
aggV <- aggData(dat_5bp_100k, 2) #<-modify the adj value for 3bp data
agg_5bp_100k <- aggV$oe
rates1 <- aggV$agg
write.table(rates1, "/net/bipolar/jedidiah/mutation/output/5bp_100k_rates.txt", col.names=T, row.names=F, quote=F, sep="\t")
agg_5bp_100k <- agg_5bp_100k[order(agg_5bp_100k$BIN),]
agg_5bp_100k$diff <- agg_5bp_100k$obs-agg_5bp_100k$exp
agg_5bp_100k$Category2 <- as.character(agg_5bp_100k$Category2)
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
}
source("agg_dat.r")
ls()
if(run_agg==1){
ptm <- proc.time()
cat("Aggregating data...\n")
source("agg_dat.r")
aggV <- aggData(dat_5bp_100k, 2) #<-modify the adj value for 3bp data
agg_5bp_100k <- aggV$oe
rates1 <- aggV$agg
write.table(rates1, "/net/bipolar/jedidiah/mutation/output/5bp_100k_rates.txt", col.names=T, row.names=F, quote=F, sep="\t")
agg_5bp_100k <- agg_5bp_100k[order(agg_5bp_100k$BIN),]
agg_5bp_100k$diff <- agg_5bp_100k$obs-agg_5bp_100k$exp
agg_5bp_100k$Category2 <- as.character(agg_5bp_100k$Category2)
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
}
if(run_agg==1){
ptm <- proc.time()
cat("Aggregating data...\n")
source("agg_dat.r")
aggV <- aggData(dat_5bp_100k, 2) #<-modify the adj value for 3bp data
agg_5bp_100k <- aggV$oe
rates1 <- aggV$agg
write.table(rates1, "/net/bipolar/jedidiah/mutation/output/5bp_100k_rates.txt", col.names=T, row.names=F, quote=F, sep="\t")
agg_5bp_100k <- agg_5bp_100k[order(agg_5bp_100k$BIN),]
agg_5bp_100k$diff <- agg_5bp_100k$obs-agg_5bp_100k$exp
agg_5bp_100k$Category2 <- as.character(agg_5bp_100k$Category2)
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
}
if(run_agg==1){
ptm <- proc.time()
cat("Aggregating data...\n")
source("agg_dat.r")
aggV <- aggData(dat_5bp_100k, 2) #<-modify the adj value for 3bp data
agg_5bp_100k <- aggV$oe
rates1 <- aggV$agg
write.table(rates1, "/net/bipolar/jedidiah/mutation/output/5bp_100k_rates.txt", col.names=T, row.names=F, quote=F, sep="\t")
agg_5bp_100k <- agg_5bp_100k[order(agg_5bp_100k$BIN),]
agg_5bp_100k$diff <- agg_5bp_100k$obs-agg_5bp_100k$exp
agg_5bp_100k$Category2 <- as.character(agg_5bp_100k$Category2)
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
}
head(summfile)
head(binfile)
binfile <- dat_5bp_100k$bin
aggseq <- count(summfile, Sequence, Category, CAT, COUNT, SEQ)
aggseq$rel_prop <- aggseq$n/aggseq$COUNT
if(run_agg==1){
ptm <- proc.time()
cat("Aggregating data...\n")
source("agg_dat.r")
aggV <- aggData(dat_5bp_100k, 2) #<-modify the adj value for 3bp data
agg_5bp_100k <- aggV$oe
rates1 <- aggV$agg
write.table(rates1, "/net/bipolar/jedidiah/mutation/output/5bp_100k_rates.txt", col.names=T, row.names=F, quote=F, sep="\t")
agg_5bp_100k <- agg_5bp_100k[order(agg_5bp_100k$BIN),]
agg_5bp_100k$diff <- agg_5bp_100k$obs-agg_5bp_100k$exp
agg_5bp_100k$Category2 <- as.character(agg_5bp_100k$Category2)
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
}
if(run_agg==1){
ptm <- proc.time()
cat("Aggregating data...\n")
source("agg_dat.r")
aggV <- aggData(dat_5bp_100k, 2) #<-modify the adj value for 3bp data
agg_5bp_100k <- aggV$oe
rates1 <- aggV$agg
write.table(rates1, "/net/bipolar/jedidiah/mutation/output/5bp_100k_rates.txt", col.names=T, row.names=F, quote=F, sep="\t")
agg_5bp_100k <- agg_5bp_100k[order(agg_5bp_100k$BIN),]
agg_5bp_100k$diff <- agg_5bp_100k$obs-agg_5bp_100k$exp
agg_5bp_100k$Category2 <- as.character(agg_5bp_100k$Category2)
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
}
if(run_agg==1){
ptm <- proc.time()
cat("Aggregating data...\n")
source("agg_dat.r")
aggV <- aggData(dat_5bp_100k, 2) #<-modify the adj value for 3bp data
agg_5bp_100k <- aggV$oe
rates1 <- aggV$agg
write.table(rates1, "/net/bipolar/jedidiah/mutation/output/5bp_100k_rates.txt", col.names=T, row.names=F, quote=F, sep="\t")
agg_5bp_100k <- agg_5bp_100k[order(agg_5bp_100k$BIN),]
agg_5bp_100k$diff <- agg_5bp_100k$obs-agg_5bp_100k$exp
agg_5bp_100k$Category2 <- as.character(agg_5bp_100k$Category2)
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
}
if(run_agg==1){
ptm <- proc.time()
cat("Aggregating data...\n")
source("agg_dat.r")
aggV <- aggData(dat_5bp_100k, 2) #<-modify the adj value for 3bp data
agg_5bp_100k <- aggV$oe
rates1 <- aggV$agg
write.table(rates1, "/net/bipolar/jedidiah/mutation/output/5bp_100k_rates.txt", col.names=T, row.names=F, quote=F, sep="\t")
agg_5bp_100k <- agg_5bp_100k[order(agg_5bp_100k$BIN),]
agg_5bp_100k$diff <- agg_5bp_100k$obs-agg_5bp_100k$exp
agg_5bp_100k$Category2 <- as.character(agg_5bp_100k$Category2)
tottime<-(proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
}
aggseq <- count(summfile, Sequence, Category, CAT, COUNT, SEQ)
aggseq$rel_prop <- aggseq$n/aggseq$COUNT
b<-c("A", "C", "G", "T")
cats<-unique(aggseq$Category)
for(i in 1:6){
aggcat <- aggseq[aggseq$Category==cats[i],]
# aggcat <- aggcat[aggcat$freq>25,]
# test <- chisq.test(aggcat$obs, p=aggcat$exp/sum(aggcat$exp))
# cat <- cats[i]
# msg <- paste0("checking category ",cat)
# print(msg)
for(j in 1:4){
for(k in 1:4){
ref <- unique(substr(aggcat$Sequence,3,3))
motif <- paste0(b[j],ref,b[k])
dat <- aggcat[substr(aggcat$Sequence,2,4)==motif,]
dat$exp <- dat$COUNT*(sum(dat$n)/sum(dat$COUNT))
# print(head(dat))
cat <- cats[i]
test <- chisq.test(dat$n, p=dat$exp/sum(dat$exp))
if(test$p.value>0.05/1536){
print(cat)
print(motif)
print(test$p.value)
# print(dat)
}
}
}
}
test$p.value
test
chisq.test(dat$n, p=dat$exp/sum(dat$exp))$p.value
head(summfile)
head(mut_cov)
log_model<-1
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
head(summfile1)
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}
if(log_model==1){
cat("Initializing logistic regression model...\n")
ptm <- proc.time()
source("log_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Model and predictions finished in", tottime, "seconds \n")
}

rates1$Seq3<-substring(rates1$Sequence,2,4)
rates.gp<-group_by(rates1, Category2, Seq3)
rates2<-summarise(rates.gp, max=max(rel_prop), min=min(rel_prop), fold=max/min)
rates3<-arrange(ungroup(rates2), desc(fold))
rates4<-group_by(rates2, Category2) %>% summarise(max=max(max), min=min(min), fold=max/min)
rates4

ls()
head(aggcat)
head(binfile)
dim(binfile)
head(mut_cov2)
head(mutcov2)
ls()
head(mut_cov)
head(arrange(mut_cov, desc(RATE)))
head(arrange(mut_cov, desc(CHR, BIN)))
head(arrange(mut_cov, asc(CHR, BIN)))
head(arrange(mut_cov, -desc(CHR, BIN)))
head(arrange(mut_cov, !desc(CHR, BIN)))
head(arrange(mut_cov, desc(CHR, BIN)))
head(arrange(mut_cov, CHR, BIN))
min(mut_cov$BIN)
head(mut_cov[mut_cov$CHR==22,])
50*100000
citation(speedglm)
citation(package="speedglm")
ls()
head(aggseq)
savehistory("/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/full_dat_session.r")
