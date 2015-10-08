##############################################################################
# Looks at singletons per 1kb in each individual
##############################################################################

dat<-read.table(" /net/bipolar/jedidiah/chr20sing.singletons", header=T, stringsAsFactors=F, sep="\t")
dat<-read.table("/net/bipolar/jedidiah/chr20sing.singletons", header=T, stringsAsFactors=F, sep="\t")
head(dat)
require(ggplot2)
binw<-1000
dat$BIN<-ceiling(dat$POS/binw)
head(dat)
dat2<-count(.~BIN+INDV)
require(plyr)
require(reshape2)
dat2<-count(.~BIN+INDV, dat)
dat2<-aggregate(.~BIN+INDV, dat)
dat2<-aggregate(.~BIN+INDV, dat, sum)
dat2<-count(dat, c("INDV", "BIN"))
head(dat2)
dim(dat)
dim(dat2)
tail(dat2)
max(dat2$freq)
which.max(dat2$freq)
dat2[which.max(dat2$freq),]
ggplot(dat2, aes(x=BIN, y=freq, colour=INDV))+geom_line()+theme(legend.position="none")
ggsave("/net/bipolar/jedidiah/images/indv_cluster.png")
ggsave("/net/bipolar/jedidiah/mutation/images/indv_cluster.png")
savehistory("indv_cluster.r")
