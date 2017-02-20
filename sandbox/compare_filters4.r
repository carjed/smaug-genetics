##############################################################################
# Used to assess quality of full data vs. heavily filtered data
# will be modifying to a proper script soon
##############################################################################

rel<-read.table("/net/bipolar/jedidiah/mutation/output/5bp_100k_full/full.summary", header=F, stringsAsFactors=F, skip=1)
names(rel)<-c("CHR", "POS", "REF", "ALT", "DP", "AN", "SEQ", "ALTSEQ", "GC")

strict<-read.table("/net/bipolar/jedidiah/mutation/output/5bp_100k_strict/full.summary", header=F, stringsAsFactors=F, skip=1)
names(strict)<-c("CHR", "POS", "REF", "ALT", "DP", "AN", "SEQ", "ALTSEQ", "GC")

binw<-100000

adj<-2
chr22<-rel
chr22$BIN <- ceiling(chr22$POS/binw)
chr22$CAT <- paste(chr22$REF, chr22$ALT, sep="")
# Manually remove bins near chr20 centromere
# chr22 <- chr22[ which(chr22$BIN<260 | chr22$BIN>300),]
chr22$Category[chr22$CAT=="AC" | chr22$CAT=="TG"] <- "AT_CG"
chr22$Category[chr22$CAT=="AG" | chr22$CAT=="TC"] <- "AT_GC"
chr22$Category[chr22$CAT=="AT" | chr22$CAT=="TA"] <- "AT_TA"
chr22$Category[chr22$CAT=="GA" | chr22$CAT=="CT"] <- "GC_AT"
chr22$Category[chr22$CAT=="GC" | chr22$CAT=="CG"] <- "GC_CG"
chr22$Category[chr22$CAT=="GT" | chr22$CAT=="CA"] <- "GC_TA"
chr22$Sequence <- ifelse(
	substr(chr22$SEQ,adj+1,adj+1) < substr(chr22$ALTSEQ,adj+1,adj+1),
	paste0(chr22$SEQ,"(",chr22$ALTSEQ,")"),
	paste0(chr22$ALTSEQ,"(",chr22$SEQ,")")
)
rel<-chr22

chr22<-strict
chr22$BIN <- ceiling(chr22$POS/binw)
chr22$CAT <- paste(chr22$REF, chr22$ALT, sep="")
# Manually remove bins near chr20 centromere
# chr22 <- chr22[ which(chr22$BIN<260 | chr22$BIN>300),]
chr22$Category[chr22$CAT=="AC" | chr22$CAT=="TG"] <- "AT_CG"
chr22$Category[chr22$CAT=="AG" | chr22$CAT=="TC"] <- "AT_GC"
chr22$Category[chr22$CAT=="AT" | chr22$CAT=="TA"] <- "AT_TA"
chr22$Category[chr22$CAT=="GA" | chr22$CAT=="CT"] <- "GC_AT"
chr22$Category[chr22$CAT=="GC" | chr22$CAT=="CG"] <- "GC_CG"
chr22$Category[chr22$CAT=="GT" | chr22$CAT=="CA"] <- "GC_TA"
chr22$Sequence <- ifelse(
	substr(chr22$SEQ,adj+1,adj+1) < substr(chr22$ALTSEQ,adj+1,adj+1),
	paste0(chr22$SEQ,"(",chr22$ALTSEQ,")"),
	paste0(chr22$ALTSEQ,"(",chr22$SEQ,")")
)
strict<-chr22

rel$Seq3<-substring(rel$Sequence,2,4)
strict$Seq3<-substring(strict$Sequence,2,4)

relct<-count(rel, Seq3)
strictct<-count(strict, Seq3)

relct$sum<-sum(relct$freq)
strictct$sum<-sum(strictct$freq)

z2<-merge(relct, strictct, by=c("Seq3"))
z2$prop.x<-z2$freq.x/z2$sum.x
z2$prop.y<-z2$freq.y/z2$sum.y

z2<-z2[,c(1,6,7)]
names(z2)<-c("Motif", "None", "Strict")

z2m<-melt(z2, id=c("Motif"))
names(z2m)<-c("Motif", "Filter", "value")
z2m$cum.perc<-(1:nrow(z2m))/nrow(z2m)

require(plyr)
df<-ddply(z2m, .(Filter), transform, cum.perc=Reduce('+', list(value/2,cumsum(c(0,head(value,-1))))))

ggplot(df, aes(Filter, value, fill=Motif, label=Motif))+
	geom_bar(stat="identity")+
	geom_text(aes(x=Filter, y=cum.perc, size=value, ymax=cum.perc))+
	scale_fill_manual(values=myPaletteCat(32), guide=F)+
	theme_bw()+
	theme(legend.position="none")+
	ylab("Relative Proportion")
	
ggsave("/net/bipolar/jedidiah/mutation/images/motif_props.png")

count<-0

for(i in 1:32){
	counts<-as.numeric(z2[i,c(2,4)])
	trials<-as.numeric(z2[i,c(3,5)])
	
	z<-prop.test(counts, trials)
	
	if(z$p.value<0.05/32){
		sequence<-z2[i,1]
		pval<-z$p.value
		cat(sequence, pval, "\n")
		# count<-count+1
	}
}

chr20b<-count(chr20, BIN)
hqb<-count(hq, BIN)
bins<-merge(chr20b, hqb, by="BIN")
bins$prop<-bins$freq.y/bins$freq.x
bins$threshold<-bins$prop<0.45

strb<-count(strict, c("CHR", "BIN"))
relb<-count(rel, c("CHR", "BIN"))

dim(strb[strb$freq<50,])
dim(relb[relb$freq<50,])

relb$dat<-"full"
strb$dat<-"strict"

datm<-merge(relb[,1:3], strb[,1:3], by=c("CHR", "BIN"), all=T)

datm[is.na(datm)]<-0

names(datm)<-c("CHR", "BIN", "full data", "filtered data")
datm2<-melt(datm, id=c("CHR", "BIN"))

ggplot(datm2, aes(value, fill=variable))+
	geom_histogram(binwidth=50, position="identity", alpha=0.5)+
	coord_cartesian(xlim=c(-100, 3000))+
	xlab("# Singletons")+
	ylab("# 100kb windows")
ggsave("/net/bipolar/jedidiah/mutation/images/countdist.png")

strict$TS<-ifelse(((strict$REF=="A" | strict$REF=="G") & (strict$ALT=="A" | strict$ALT=="G")) | ((strict$REF=="C" | strict$REF=="T") & (strict$ALT=="C" | strict$ALT=="T")), "TS", "TV")
rel$TS<-ifelse(((rel$REF=="A" | rel$REF=="G") & (rel$ALT=="A" | rel$ALT=="G")) | ((rel$REF=="C" | rel$REF=="T") & (rel$ALT=="C" | rel$ALT=="T")), "TS", "TV")

tstv.str<-count(strict, c("CHR", "BIN", "TS"))
tstv.str2<-dcast(tstv.str, CHR+BIN~TS, value.var="freq")
tstv.str2$TOT<-tstv.str2$TS+tstv.str2$TV
tstv.str2$TSTV<-tstv.str2$TS/tstv.str2$TV
tstv.str2$BINCAT<-ceiling(tstv.str2$TOT/50)*50
tstv.str2$filter<-"strict"

tstv.rel<-count(rel, c("CHR", "BIN", "TS"))
tstv.rel2<-dcast(tstv.rel, CHR+BIN~TS, value.var="freq")
tstv.rel2$TOT<-tstv.rel2$TS+tstv.rel2$TV
tstv.rel2$TSTV<-tstv.rel2$TS/tstv.rel2$TV
tstv.rel2$BINCAT<-ceiling(tstv.rel2$TOT/50)*50
tstv.rel2$filter<-"none"

tstv.full<-rbind(tstv.str2, tstv.rel2)
# tstv.full<-tstv.full[tstv.full$BINCAT<=60,]
tstv.full<-tstv.full[tstv.full$TSTV<=8,]
# tstv.full<-tstv.full[complete.cases(tstv.full),]

update_geom_defaults("point", list(colour = NULL))
ggplot(tstv.full, aes(x=factor(BINCAT), y=TSTV, colour=filter))+
	geom_boxplot(outlier.colour=NULL)+
	scale_y_continuous(breaks=1:8)+
	geom_hline(aes(yintercept=2), colour="red")+
	# facet_wrap(~filter)+
	xlab("Singletons per bin")+
	theme_bw()+
	theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave("/net/bipolar/jedidiah/mutation/images/strict_tstv.png", width=14, height=7)

tstv.full2<-aggregate(TSTV~filter+CHR+BIN, data=tstv.full, FUN=mean)
	ggplot(tstv.full2, aes(x=factor(BINCAT), y=TSTV, fill=filter))+
	geom_bar(stat="identity", position="dodge")+
	scale_y_continuous(breaks=seq(0, 2.5, by=0.1))+
	geom_hline(aes(yintercept=2), colour="red")+
	# facet_wrap(~filter)+
	xlab("Singletons per bin")+
	theme_bw()+
	theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave("/net/bipolar/jedidiah/mutation/images/strict_tstv_bar.png", width=14, height=7)

tstv.full2<-aggregate(TSTV~filter+CHR+BIN, data=tstv.full, FUN=mean)
ggplot(tstv.full2[tstv.full2$CHR==2,], aes(x=BIN, y=TSTV, group=filter, colour=filter))+
	geom_point(alpha=0.5)+
	# geom_line()+
	scale_y_continuous(breaks=seq(0, 10, by=0.5))+
	geom_hline(aes(yintercept=2), colour="red")+
	# facet_wrap(~CHR)+
	xlab("Bin")+
	theme_bw()+
	theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave("/net/bipolar/jedidiah/mutation/images/strict_tstv_chr.png", width=14, height=7)
 
grouped <- group_by(tstv.full, filter, BINCAT)
gpsum<-ddply(tstv.full, c("filter", "BINCAT"), mean=mean(TSTV), sd=sd(TSTV))

nrow(tstv.full[(tstv.full$filter=="strict" & tstv.full$TSTV<1),])/27001
nrow(tstv.full[(tstv.full$filter=="none" & tstv.full$TSTV<1),])/27001

d1<-tstv.full[(tstv.full$filter=="strict" & tstv.full$TSTV<1),]
d2<-tstv.full[(tstv.full$filter=="none" & tstv.full$TSTV<1),]

prop.test(c(nrow(tstv.full[(tstv.full$filter=="strict" & tstv.full$TSTV<1.5),]), nrow(tstv.full[(tstv.full$filter=="none" & tstv.full$TSTV<1.5),])), c(27001, 27001))
