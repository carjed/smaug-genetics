#! /usr/bin/env Rscript

##############################################################################
# Compare correlation between 100kb windowed distribution of singletons
# and common variants.
# -Uses data now stored in /net/bipolar/jedidiah/mutation/output/supp/
##############################################################################

##############################################################################
# Load packages and define color palettes
##############################################################################
suppressMessages(require(ggplot2))
suppressMessages(require(plyr))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
rb<-c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g<-myPaletteG(6)[1:3]

##############################################################################
# Read data
##############################################################################
# sing<-read.table("expanded.summary", header=T, stringsAsFactors=F)
# spec<-read.table("chr10_spec_expanded.summary", header=T, stringsAsFactors=F)

# sing_bin<-read.table("bin_out_sing.txt", header=T, stringsAsFactors=F, check.names=F)
# spec_bin<-read.table("bin_out_spec.txt", header=T, stringsAsFactors=F, check.names=F)

chr10<-read.table("chr10.expanded.summary", header=T)
sing_bin<-read.table("chr10.bin_out.txt", header=T, check.names=F)
spec_bin<-sing_bin

output<-matrix(ncol=6,nrow=1000)

for(i in 1:1000){
	chr10$gp<-sample(0:1, nrow(chr10), replace=T)

	sing<-chr10[chr10$gp==0,]
	spec<-chr10[chr10$gp==1,]



	##############################################################################
	# Add category and bin columns
	##############################################################################
	binw<-100000
	sing$BIN<-ceiling(sing$POS/binw)

	# chr22<-sing
	sing$CAT<-paste(sing$REF, sing$ALT, sep="")
	sing$Category[sing$CAT=="AC" | sing$CAT=="TG"]<-"AT_CG"
	sing$Category[sing$CAT=="AG" | sing$CAT=="TC"]<-"AT_GC"
	sing$Category[sing$CAT=="AT" | sing$CAT=="TA"]<-"AT_TA"
	sing$Category[sing$CAT=="GA" | sing$CAT=="CT"]<-"GC_AT"
	sing$Category[sing$CAT=="GC" | sing$CAT=="CG"]<-"GC_CG"
	sing$Category[sing$CAT=="GT" | sing$CAT=="CA"]<-"GC_TA"

	xmax<-floor(max(sing$BIN)/100)*100

	spec$BIN<-ceiling(spec$POS/binw)
	spec$CAT<-paste(spec$REF, spec$ALT, sep="")
	spec$Category[spec$CAT=="AC" | spec$CAT=="TG"]<-"AT_CG"
	spec$Category[spec$CAT=="AG" | spec$CAT=="TC"]<-"AT_GC"
	spec$Category[spec$CAT=="AT" | spec$CAT=="TA"]<-"AT_TA"
	spec$Category[spec$CAT=="GA" | spec$CAT=="CT"]<-"GC_AT"
	spec$Category[spec$CAT=="GC" | spec$CAT=="CG"]<-"GC_CG"
	spec$Category[spec$CAT=="GT" | spec$CAT=="CA"]<-"GC_TA"

	##############################################################################
	# Calculate subtype-specific relative rates per bin
	##############################################################################
	chr22b<-merge(sing, sing_bin, by="BIN", all=TRUE)
	countsing<-count(chr22b, c("BIN", "AT", "CG", "prop_GC", "Category"))
	countsing<-merge(countsing, aggregate(freq~BIN, data=countsing, sum), by="BIN", all=TRUE)

	countsing<-countsing[order(countsing$BIN),]


	chr22c<-merge(spec, spec_bin, by="BIN", all=TRUE)
	countspec<-count(chr22c, c("BIN", "AT", "CG", "prop_GC", "Category"))
	countspec<-merge(countspec, aggregate(freq~BIN, data=countspec, sum), by="BIN", all=TRUE)

	countspec<-countspec[order(countspec$BIN),]

	names(countsing)[names(countsing)=="freq.x"]<-"sing_count"
	names(countspec)[names(countspec)=="freq.x"]<-"spec_count"

	##############################################################################
	# Merge and sort combined data
	##############################################################################
	fullcount<-merge(countsing, countspec, by=c("BIN", "Category"))
	fullcount<-fullcount[order(fullcount$BIN),]
	fullcount<-fullcount[complete.cases(fullcount),]

	##############################################################################
	# Calculate correlation in each subtype and create df
	##############################################################################
	cors<-c(
			cor(fullcount[fullcount$Category=="AT_CG",]$sing_count, fullcount[fullcount$Category=="AT_CG",]$spec_count),
			cor(fullcount[fullcount$Category=="AT_GC",]$sing_count, fullcount[fullcount$Category=="AT_GC",]$spec_count),
			cor(fullcount[fullcount$Category=="AT_TA",]$sing_count, fullcount[fullcount$Category=="AT_TA",]$spec_count),
			cor(fullcount[fullcount$Category=="GC_AT",]$sing_count, fullcount[fullcount$Category=="GC_AT",]$spec_count),
			cor(fullcount[fullcount$Category=="GC_CG",]$sing_count, fullcount[fullcount$Category=="GC_CG",]$spec_count),
			cor(fullcount[fullcount$Category=="GC_TA",]$sing_count, fullcount[fullcount$Category=="GC_TA",]$spec_count)
		)

	# cors<-paste0("r=",cors)

	output[i,]<-cors

	# xpos<-rep(median(fullcount$sing_prop), 6)
	# ypos<-0.9*rep(max(fullcount$sing_prop), 6)

	# cats<-sort(unique(fullcount$Category))

	# dat<-data.frame(xpos, ypos, cats, cors)
	# names(dat)<-c("x", "y", "Category", "val")
}



##############################################################################
# Get mean and cis for each category
##############################################################################
cis<-matrix(nrow=6, ncol=4)
for(i in 1:6){

	cis[i,]<-c(names(output)[i], mean(output[,i]), t.test(output[,i])$conf.int[1], t.test(output[,i])$conf.int[2])

}

cis<-data.frame(cis)
names(cis)<-c("Category", "mean", "lcl", "ucl")

##############################################################################
# Scatter plot of relative rates in singleton vs. divergent sites
##############################################################################
# my_grob = grobTree(textGrob("This text stays in place!", x=0.1,  y=0.95, hjust=0,
  # gp=gpar(col="blue", fontsize=12, fontface="italic")))

# ggplot()+
	# geom_point(data=fullcount, aes(x=sing_prop, y=spec_prop, colour=Category, group=Category), alpha=0.4)+
	# scale_colour_manual(values=rb)+
	# xlab("Singleton Relative Rate")+
	# ylab("Divergent Site Relative Rate")+
	# xlab("Group 1 Relative Rate")+
	# ylab("Group 2 Relative Rate")+
	# theme_bw()+
	# geom_text(data=dat, aes(x=-Inf,y=Inf, label=val, hjust=0, vjust=1))+
	# facet_wrap(~Category)
# ggsave("sing_spec_compare_no_sub.png", width=8.4, height=4.4)
