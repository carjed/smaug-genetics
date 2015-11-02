suppressMessages(require(ggplot2))
suppressMessages(require(plyr))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))
binw<-100000
myPaletteCat <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
rb<-c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g<-myPaletteG(6)[1:3]
chr22<-read.table("chr10.expanded.summary", header=T, stringsAsFactors=F)
bins<-read.table("chr10.bin_out.txt", header=T, stringsAsFactors=F, check.names=F)
chr22$BIN<-ceiling(chr22$POS/binw)
chr22$CAT<-paste(chr22$REF, chr22$ALT, sep="")
# chr22<-chr22[ which(chr22$BIN<260 | chr22$BIN>300),]
chr22$Category[chr22$CAT=="AC" | chr22$CAT=="TG"]<-"AT_CG"
chr22$Category[chr22$CAT=="AG" | chr22$CAT=="TC"]<-"AT_GC"
chr22$Category[chr22$CAT=="AT" | chr22$CAT=="TA"]<-"AT_TA"
chr22$Category[chr22$CAT=="GA" | chr22$CAT=="CT"]<-"GC_AT"
chr22$Category[chr22$CAT=="GC" | chr22$CAT=="CG"]<-"GC_CG"
chr22$Category[chr22$CAT=="GT" | chr22$CAT=="CA"]<-"GC_TA"
xmax<-floor(max(chr22$BIN)/100)*100
exons<-read.table("exon_loc.bed", header=F, sep="\t", stringsAsFactors=F)
#rename cols
names(exons)<-c("chr", "start", "end", "gene", "V5", "v6")
#get bins and identify mismatches
exons$bin.a<-ceiling(exons$start/100000)
exons$bin.b<-ceiling(exons$end/100000)
#subset to get unique exon positions
exons<-unique(exons[,c(1,2,3,5,6,7,8)])
#get number exonic sites for bins w/o overlap
exons$length<-exons$end-exons$start+1
subeq<-exons[exons$bin.a==exons$bin.b,]
z<-ddply(subeq, c("bin.a"), summarize, tot=sum(length))
#get number exonic sites for bins w/ overlap
subneq<-exons[exons$bin.a!=exons$bin.b,]
subneq$length.a<-(subneq$bin.a*100000)-subneq$start
subneq$length.b<-subneq$end-(subneq$bin.a*100000)+1
za<-ddply(subneq, c("bin.a"), summarize, tot=sum(length.a))
zb<-ddply(subneq, c("bin.b"), summarize, tot=sum(length.b))
#bind data
names(zb)<-c("bin.a", "tot")
zc<-rbind(za, zb)
#add to bins w/o overlap
zd<-ddply(merge(z, zc, all=T), .(bin.a), summarize, tot=sum(tot))
#rename cols
names(zd)<-c("bin", "tot")
#get %exonic
zd$pct<-zd$tot/100000
zd$qt1<-ifelse(zd$pct<=quantile(zd$pct)[2], 1, 0)
zd$qt3<-ifelse(zd$pct>=quantile(zd$pct)[4], 1, 0)
zdmin<-zd[,c(1,3,4,5)]
chr22z<-merge(chr22, zdmin, by.x="BIN", by.y="bin")
chr22<-chr22z[chr22z$qt1==1,]
head(chr22)
ggplot(chr22, aes(x=POS))+
geom_histogram(binwidth=binw, position="identity", alpha=0.5, aes(colour=pct, fill=pct))+
ggtitle("main_dist_title2")+
scale_fill_gradientn(colours=myPalette(100))+
theme_bw()+
# scale_y_continuous(trans = "log")+
theme(panel.border=element_blank(),
legend.position="none",
axis.text.x = element_text(angle = 90, hjust = 0.5),
axis.text.y = element_text(angle = 90, hjust = 0.5))
savehistory("testrun.R")
