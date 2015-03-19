suppressMessages(require(ggplot2))
suppressMessages(require(plyr))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))
	myPaletteCat <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
rb<-c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g<-myPaletteG(6)[1:3]
binw<-100000
chr22<-read.table("chr22.expanded.summary", header=T, stringsAsFactors=F)
bins<-read.table("chr22.bin_out.txt", header=T, stringsAsFactors=F, check.names=F)
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
chr<-"22"
mac<-"singleton"
imgdir<-"/net/bipolar/jedidiah/images/chr22"
{
main_rel_title<-paste0("Chr",chr," ",mac, " Relative Mutation Rate")
main_rel_out<-paste0(imgdir,"/chr",chr,"_",mac,"_mutation_prop.png")
main_scale_title<-paste0("Chr",chr," ",mac, " Mutations--scaled proportions")
main_scale_out<-paste0(imgdir,"/chr",chr,"_",mac,"_mutation_prop2.png")
main_dist_title<-paste0("Chr",chr," ",mac, " Distribution by Mutation Type")
main_dist_out<-paste0(imgdir,"/chr",chr,"_",mac,"_dist_count.png")
main_dist_title2<-paste0("Distribution of Chr",chr," singletons (low exon density bins)")
main_dist_out2<-paste0(imgdir,"/chr",chr,"_",mac,"_dist_count_all.png")
}
ggplot(chr22, aes(x=POS, colour=Category, fill=Category))+
geom_histogram(binwidth=binw, position="identity", alpha=0.5)+
facet_wrap(~Category)+
scale_fill_manual(values=rb)+
scale_colour_manual(values=rb)+
ggtitle(main_dist_title)+
theme_bw()+
# scale_y_log10()+
theme(panel.border=element_blank(),
legend.position="none",
axis.text.x = element_text(angle = 90, hjust = 0.5),
axis.text.y = element_text(angle = 90, hjust = 0.5))
suppressMessages(ggsave(main_dist_out))
chr22b<-merge(chr22, bins, by="BIN")
# count<-count(chr22b, c("Category", "BIN", "AT", "CG", "prop_GC", "pct"))
count<-count(chr22b, c("Category", "BIN", "AT", "CG", "prop_GC"))
rm(chr22b)
count<-merge(count, aggregate(freq~BIN, data=count, sum), by="BIN", all=TRUE)
count$rel_prop<-count$freq.x/count$freq.y
count<-count[order(count$BIN, count$Category),]
ggplot(count, aes(x=factor(BIN), y=rel_prop, colour=Category, fill=Category))+
geom_bar(position="stack", stat="identity", alpha=0.5)+
# scale_x_discrete(breaks=seq(0,xmax,50))+
scale_fill_manual(values=rb)+
scale_colour_manual(values=rb)+
xlab("Bin")+
ylab("Proportion")+
ggtitle(main_scale_title)+
theme_bw()+
theme(panel.border=element_blank(),
axis.text.x = element_text(angle = 90, hjust = 0.5),
axis.text.y = element_text(angle = 90, hjust = 0.5))
suppressMessages(ggsave(main_scale_out))
countAT<-count[grep("^AT", count$Category),]
countGC<-count[grep("^GC", count$Category),]
countNULL<-count[count$prop_GC==0,]
countAT$prop<-countAT$freq.x/countAT$AT
countGC$prop<-countGC$freq.x/countGC$CG
countNULL$prop<-countNULL$AT/countNULL$freq.x
count2<-rbind(countAT, countGC, countNULL)
count2[is.na(count2)]<-"AT_CG"
count2$prop[count2$prop>0.5]<-0
AT_s<-aggregate(freq.x~Category, data=countAT, sum)
AT_t<-aggregate(AT~Category, data=countAT, sum)
AT_s$prop<-AT_s$freq.x/AT_t$AT
GC_s<-aggregate(freq.x~Category, data=countGC, sum)
GC_t<-aggregate(CG~Category, data=countGC, sum)
GC_s$prop<-GC_s$freq.x/GC_t$CG
cwa<-rbind(AT_s, GC_s)
names(cwa)<-c("variable", "freq.x", "avg")
ggplot(count2, aes(x=factor(BIN), y=prop, colour=Category, fill=Category))+
geom_bar(position="identity", stat="identity", alpha=0.5)+
facet_wrap(~Category)+
# scale_x_discrete(breaks=seq(0,xmax,50))+
scale_fill_manual(values=rb)+
scale_colour_manual(values=rb)+
xlab("Bin")+
ylab("Relative Mutation Rate")+
#ggtitle(main_rel_title)+
theme_bw()+
theme(panel.border=element_blank(),
legend.position="none",
axis.text.x = element_text(angle = 90, hjust = 0.5),
axis.text.y = element_text(angle = 90, hjust = 0.5))
suppressMessages(ggsave(main_rel_out))
chr22$Sequence<-paste0(pmin(chr22$SEQ, chr22$ALTSEQ),"(",pmax(chr22$SEQ, chr22$ALTSEQ),")")

chr22$Sequence<-ifelse(
	substr(chr22$SEQ,2,2)<substr(chr22$ALTSEQ,2,2),
	paste0(chr22$SEQ,"(",chr22$ALTSEQ,")"),
	paste0(chr22$ALTSEQ,"(",chr22$SEQ,")")
)

sort(unique(chr22$Sequence))

revcomp = function(DNAstr) {
	step1 = chartr("ACGT","TGCA",DNAstr)
	step2 = unlist(strsplit(step1, split=""))
	step3 = rev(step2)
	step4 = paste(step3, collapse="")
	return(step4)
}

for(i in 5:68){
	names(bins)[i]<-paste0(names(bins)[i], "(", revcomp(names(bins)[i]), ")" )
}

bins2<-bins[,names(bins)%in%unique(chr22$Sequence)]


)paste0(pmin(chr22$SEQ, chr22$ALTSEQ),"(",pmax(chr22$SEQ, chr22$ALTSEQ),")")

cats<-factor(chr22$Category)
#Output datasets and heatmaps
chr22r<-chr22[,c('BIN', 'Category', 'Sequence')]
#Counts per bin
pc<-dcast(chr22r, BIN~Category+Sequence)
pcm<-merge(bins, pc, by="BIN", all=TRUE)
pcm<-pcm[,names(pc)]
head(pcm)
savehistory("new_trinuc.R")
