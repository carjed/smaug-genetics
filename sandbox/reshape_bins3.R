chr22<-read.table("chr19.expanded.summary", header=T, stringsAsFactors=F) 
bins<-read.table("chr19.bin_out.txt", header=T, stringsAsFactors=F) 
binw<-100000
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
chr22$Sequence<-paste0(pmin(chr22$SEQ, chr22$ALTSEQ),"(",pmax(chr22$SEQ, chr22$ALTSEQ),")")
cats<-factor(chr22$Category)
chr22$Sequence<-ifelse(chr22$REF == pmin(substr(chr22$SEQ,adj+1,adj+1), substr(chr22$ALTSEQ,adj+1,adj+1)), chr22$SEQ,chr22$ALTSEQ)
bins2<-melt(bins[,4:((4^(adj*2+1))/2+4)], id="BIN")
bins2<-aggregate(data=bins2, value ~ variable, sum)
names(bins2)<-c("Sequence", "COUNT")
bins2$Sequence<-sub("[.]", "(", bins2$Sequence)
bins2$Sequence<-sub("[.]", ")", bins2$Sequence)
require(reshape2)
chr22$Sequence<-ifelse(chr22$REF == pmin(substr(chr22$SEQ,adj+1,adj+1), substr(chr22$ALTSEQ,adj+1,adj+1)), chr22$SEQ,chr22$ALTSEQ)
bins2<-melt(bins[,4:((4^(adj*2+1))/2+4)], id="BIN")
bins2<-aggregate(data=bins2, value ~ variable, sum)
names(bins2)<-c("Sequence", "COUNT")
bins2$Sequence<-sub("[.]", "(", bins2$Sequence)
bins2$Sequence<-sub("[.]", ")", bins2$Sequence)
adj<-2
chr22$Sequence<-ifelse(chr22$REF == pmin(substr(chr22$SEQ,adj+1,adj+1), substr(chr22$ALTSEQ,adj+1,adj+1)), chr22$SEQ,chr22$ALTSEQ)
bins2<-melt(bins[,4:((4^(adj*2+1))/2+4)], id="BIN")
bins2<-aggregate(data=bins2, value ~ variable, sum)
names(bins2)<-c("Sequence", "COUNT")
bins2$Sequence<-sub("[.]", "(", bins2$Sequence)
bins2$Sequence<-sub("[.]", ")", bins2$Sequence)
chr22<-merge(chr22, bins2, by="Sequence")
head(chr22)
head(bins2)
chr22<-read.table("chr19.expanded.summary", header=T, stringsAsFactors=F) 
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
head(chr22)
chr22$Sequence<-ifelse(chr22$REF == pmin(substr(chr22$SEQ,adj+1,adj+1), substr(chr22$ALTSEQ,adj+1,adj+1)), chr22$SEQ,chr22$ALTSEQ)
head(chr22)
tail(bins2)
max(chr22$Sequence)
chr22[,chr22$Sequence=="TTTTT"]
chr22[,chr22$Sequence=="AAAAA"]
chr22[,chr22$Sequence="AAAAA"]
chr22[,chr22$Sequence == "AAAAA"]
chr22[,chr22$Sequence == "TTCTT"]
dim(chr22)
head(chr22)
chr22[,chr22$Sequence=="CCAGT"]
chr22[chr22$Sequence=="TTTTT", ]
chr22[chr22$Sequence=="TTCTT", ]
chr22[chr22$Sequence=="AAAAA", ]
dim(chr22)
head(chr22)
chr22[chr22$Sequence=="TTCTT", ]
chr22[419224,]
chr22[chr22$Sequence=="AAGAA", ]
bins3<-read.table("chr19.bin_out2.txt", header=T, stringsAsFactors=F)
chr22m<-merge(chr22, bins2, by="SEQ")
chr22m<-merge(chr22, bins3, by="SEQ")
head(chr22m)
tal(chr22m)
tail(chr22m)
chr22m2<-merge(chr22m, bins3, by.x="ALTSEQ", by.y="SEQ")
max(chr22m$Sequence)
head(chr22m2)
chr22m2$COUNT=chr22m2$COUNT.x+chr22m2$COUNT.y
head(chr22m2)
head(bins2)
tail(bins2)
length(unique(chr22$Sequence))
length(unique(chr22m2$Sequence))
head(chr22)
chr22[chr22$POS==89309,]
chr22m2[chr22m2$POS==89309,]
bins2$SEQ1<-substr(bins2$Sequence, 0, adj*2+1)
head(bins2)
bins2$SEQ2<-substr(bins2$Sequence, 6, adj*2+1)
head(bins2)
bins2$SEQ2<-substr(bins2$Sequence, 6, 6+adj*2+1)
head(bins2)
bins2$SEQ2<-substr(bins2$Sequence, 7, 7+adj*2+1)
head(bins2)
bins2$SEQ2<-substr(bins2$Sequence, 7, 6+adj*2+1)
head(bins2)
bins2$SEQ2<-substr(bins2$Sequence, (adj*2+1)+2, (adj*2+1)+(adj*2+1))
head(bins2)
bins2$SEQ2<-substr(bins2$Sequence, (adj*2+1)+2, (adj*2+2)+(adj*2+1))
head(bins2)
z<-"AAA(TTT)"
substr(z, 3+2, 4+3)
chr22m3<-merge(chr22, bins2, by.x="Sequence", by.y="SEQ")
chr22m3<-merge(chr22, bins2, by.x="Sequence", by.y="SEQ1")
head(chr22m3)
tail(chr22m3)
dim(chr22)
dim(chr22m3)
dim(chr22m2)
tail(chr22m3)
head(bins2)
ls()
bins4<-bins2[,2:4]
dim(chr22)
chr22m3<-merge(chr22, bins2, by.x="Sequence", by.y="SEQ1")
chr22m3<-merge(chr22, bins4, by.x="Sequence", by.y="SEQ1")
dim(chr22m3)
chr22m3<-merge(chr22, bins4, by.x="Sequence", by.y=c("SEQ1", "SEQ2"))
max(bins2$SEQ1)
chr22m3<-merge(chr22, bins4, by.x=c("SEQ", "ALTSEQ"), by.y=c("SEQ1", "SEQ2"))
dim(chr22m3)
dim(chr22)
chr22$COUNT<-ifelse(chr22$Sequence==(bins2$SEQ1 | bins2$SEQ2), bins2$COUNT, 0)
chr22$COUNT<-ifelse(chr22$Sequence==bins2$SEQ1, bins2$COUNT, 0)
head(chr22)
chr22$COUNT<-ifelse((chr22$SEQ | chr22$ALTSEQ)==bins2$SEQ1, bins2$COUNT, 0)
chr22$COUNT<-ifelse((chr22$SEQ==bins2$SEQ1) | (chr22$ALTSEQ==bins2$SEQ1), bins2$COUNT, 0)
chr22m3<-merge(chr22, bins4, by.x="Sequence", by.y="SEQ1", all.x=T)
head(chr22m3)
min(chr22m3$COUNT.y)
min(chr22m3$COUNT.x)
max(chr22m3$COUNT.x)
max(chr22m3$COUNT.y)
tail(chr22m3)
chr22m3<-merge(chr22m3, bins4, by.x="Sequence", by.y="SEQ2", all.x=T)
tail(chr22m3)
head(chr22m3)
head(chr22)
head(bins2)
bins2$SEQ<-pmin(bins2$SEQ1, bins2$SEQ2)
head(bins2)
chr22$SEQMIN<-pmin(chr22$SEQ, chr22$ALTSEQ)
chr22m3<-merge(chr22, bins2, by.x="SEQMIN", by.y="SEQ")
dim(chr22m3)
dim(chr22m2)
savehistory("reshape_bins3.R")
