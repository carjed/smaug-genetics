chr22<-read.table("chr19.expanded.summary", header=T, stringsAsFactors=F)
bins<-read.table("chr19.bin_out.txt", header=T, stringsAsFactors=F)
head(bins)
require(plyr)
require(reshape2)
bins2<-colsums(bins)
bins2<-colSums(bins)
dim(bins2)
bins2
bins2<-melt(bins, id=names()
dim(bins2)
bins2[1:3,]
dim(bins)
bins[1:4,4:8]
516*3
bins2<-melt(bins, id=names(bins)[5:516])
dim(bins2)
bins2[1:4,1:4]
bins2[1:4,512:514]
colSums[bins2][1]
names(bins)[1:4]
bins2<-melt(bins, id=names(bins)[5:516], variable="BIN")
dim(bins2)
bins2[1:4,1:4]
bins2[1:4,512:514]
sum(bins$AAAAA.TTTTT.)
411637+421009
bins2<-aggregate(bins, by=names(bins)[4:516], sum)
bins2<-aggregate(bins, by=, sum)
class(names(bins))
bins2<-aggregate(bins, by=as.list(names(bins)[5:516]), sum)
bins2<-aggregate(bins, by=colnames(bins)[5:516], sum)
bins2<-bins[, lapply(sum), by=colnames(bins)[5:516]]
bins2<-bins[, lapply(.SD,sum), by=colnames(bins)[5:516]]
bins2<-data.frame(colSums(bins))
dim(bins2)
bins2[1:4,]
bins2<-colSums(bins)
bins2[1]
bins2[5]
class(bins2)
names(bins2)[1]
bins2<-data.frame(as.list(colSums(bins)[5:516]))
dim(bins2)
bins2[1:4,]
dim(bins2)
bins2[,1:4]
bins2<-data.frame(t(as.list(colSums(bins)[5:516])))
dim(bins2)
bins2<-data.frame(as.list(t(colSums(bins)[5:516])))
dim(bins2)
bins2[,1:4]
bins2<-data.frame(t(as.list(colSums(bins)[5:516])))
bins2[,1:4]
bins3<-melt(bins2)
dim(bins3)
bins3<-t(bins2)
dim(bins3)
bins3[1:4,]
head(bins3)
class(bins3)
bins2<-t(data.frame(as.list(colSums(bins)[5:516])))
class(bins2)
bins2[1]
bins2<-data.frame(t(data.frame(as.list(colSums(bins)[5:516]))))
bins2[1]
bins2[1:4,]
names(bins2)
dim(bins2)
bins2<-data.frame(t(data.frame(as.list(colSums(bins)[5:516]))), row.names=F)
savehistory("reshape_bins.R")
