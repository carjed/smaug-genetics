##############################################################################
# Process Options/Args
# Define color palettes
# Define functions
##############################################################################
suppressMessages(require(ggplot2))
suppressMessages(require(plyr))
suppressMessages(require(reshape2))
suppressMessages(require(RColorBrewer))
suppressMessages(require(MASS))
suppressMessages(require(grid))

source("/net/bipolar/jedidiah/mutation/smaug-genetics/get_functions.R")

adj <- 2
nbp <- adj*2+1
binw <- 100000

myPaletteCat <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
rb <- c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g <- myPaletteG(6)[1:3]

summfile <- paste0("/net/bipolar/jedidiah/mutation/output/5bp_100k/chr2.expanded.summary")
binfile <- paste0("/net/bipolar/jedidiah/mutation/output/5bp_100k/chr2.bin_out.txt")
chr22 <- read.table(summfile, header=T, stringsAsFactors=F)
bins <- read.table(binfile, header=T, stringsAsFactors=F)

# chr22 <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.expanded.summary", header=T, stringsAsFactors=F)
# bins <- read.table("/net/bipolar/jedidiah/mutation/output/chr20.bin_out.txt", header=T, stringsAsFactors=F)

# chr22 <- chr22[-grep(",", chr22$ALT),]

##############################################################################
# Add columns to data
##############################################################################
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
	substr(chr22$SEQ,adj+1,adj+1)<substr(chr22$ALTSEQ,adj+1,adj+1),
	paste0(chr22$SEQ,"(",chr22$ALTSEQ,")"),
	paste0(chr22$ALTSEQ,"(",chr22$SEQ,")")
)

# get complement of sequence columns in bin file and remove duplicates
for(i in 6:((4^(adj*2+1))+4)){
	names(bins)[i] <- paste0(names(bins)[i], "(", revcomp(names(bins)[i]), ")" )
}

bins2 <- bins[,names(bins)%in%unique(chr22$Sequence)]
bins <- cbind(bins[,1:5],bins2)
xmax <- floor(max(chr22$BIN)/100)*100

bins2 <- melt(bins[,5:((4^(adj*2+1))/2+4)], id="BIN")
bins2 <- aggregate(data=bins2, value ~ variable, sum)
names(bins2) <- c("Sequence", "COUNT")
bins2$Sequence <- sub("[.]", "(", bins2$Sequence)
bins2$Sequence <- sub("[.]", ")", bins2$Sequence)
bins2$SEQ1 <- substr(bins2$Sequence, 0, adj*2+1)
bins2$SEQ2 <- substr(bins2$Sequence, (adj*2+1)+2, (adj*2+2)+(adj*2+1))
bins2$SEQMIN <- pmin(bins2$SEQ1, bins2$SEQ2)
bins2 <- data.frame(bins2$COUNT, bins2$SEQMIN)
names(bins2) <- c("COUNT", "SEQMIN")

chr22$SEQMIN <- pmin(chr22$SEQ, chr22$ALTSEQ)
chr22 <- merge(chr22, bins2, by="SEQMIN")

aggseq <- count(chr22, c("Sequence", "Category", "COUNT"))
aggseq$rel_prop <- aggseq$freq/aggseq$COUNT

##############################################################################
# Get observed counts for each category/bin
##############################################################################
ct <- count(chr22,c("BIN","Category"))
ct.ord <- ct[order(ct$Category),]

##############################################################################
# Get expected counts for each 
##############################################################################

# Transpose bin data to merge with 
binsT <- setNames(data.frame(t(bins[,-c(1:4)])), paste0("Bin",seq(1,nrow(bins))))
binsT$Sequence <- rownames(binsT)

# merge counts per sequence/category with counts of mutable motifs
aggseq.m <- merge(aggseq, binsT, by="Sequence")

# get expected counts for each row
aggseq.m[,6:ncol(aggseq.m)] <- aggseq.m$rel_prop*aggseq.m[,6:ncol(aggseq.m)]

# sum all sequence combinations for each category/bin
atgc <- colSums(aggseq.m[aggseq.m$Category=="AT_GC",-c(1:5)])
atcg <- colSums(aggseq.m[aggseq.m$Category=="AT_CG",-c(1:5)])
atta <- colSums(aggseq.m[aggseq.m$Category=="AT_TA",-c(1:5)])
gcat <- colSums(aggseq.m[aggseq.m$Category=="GC_AT",-c(1:5)])
gccg <- colSums(aggseq.m[aggseq.m$Category=="GC_CG",-c(1:5)])
gcta <- colSums(aggseq.m[aggseq.m$Category=="GC_TA",-c(1:5)])

# vector of expected counts
exp <- c(atcg, atgc, atta, gcat, gccg, gcta)

# vector of categories
catrep <- c(rep("AT_CG",ncol(aggseq.m)-5), rep("AT_GC",ncol(aggseq.m)-5), rep("AT_TA",ncol(aggseq.m)-5), 
          rep("GC_AT",ncol(aggseq.m)-5), rep("GC_CG",ncol(aggseq.m)-5), rep("GC_TA",ncol(aggseq.m)-5))
		  
# create data frame with observed, expected, and bin
oe2 <- data.frame(exp, Category=catrep, BIN=rep(1:(ncol(aggseq.m)-5),6))
oe2 <- merge(oe2,ct.ord, by=c("Category","BIN"))
names(oe2)[4] <- "obs"

# get odds ratio for each bin/category
oe2$odds <- oe2$obs/oe2$exp

# order by OR and add column of ranks
oe2.ord <- oe2[order(oe2$Category, oe2$odds),]
# oe2.ord$rk <- rep(1:length(unique(oe2.ord$BIN)),6)
oe2.ord <- ddply(oe2.ord, .(Category), transform, rk=seq_along(Category))

# Plot odds profiles of each category
ggplot(oe2.ord, aes(x=rk, y=odds))+
	geom_point()+
	facet_wrap(~Category)
ggsave("/net/bipolar/jedidiah/mutation/images/odds_profiles.png")


oe2.ord$res <- paste0(nbp,"bp")
data_out <- paste0("oe_",nbp,"bp.txt")
# write.table(oe2.ord, data_out, sep="\t", quote=F, row.names=F, col.names=T)

##############################################################################
# Get 1bp predictions
##############################################################################
aggseq6<-ddply(aggseq, .(Category), summarize, COUNT=sum(COUNT), freq=sum(freq))
aggseq6$rel_prop<-aggseq6$freq/aggseq6$COUNT

atcg<-aggseq6[1,4]*bins[,1]
atgc<-aggseq6[2,4]*bins[,1]
atta<-aggseq6[3,4]*bins[,1]
gcat<-aggseq6[4,4]*bins[,2]
gccg<-aggseq6[5,4]*bins[,2]
gcta<-aggseq6[6,4]*bins[,2]

exp<-c(atcg, atgc, atta, gcat, gccg, gcta)

# vector of categories
catrep<-c(rep("AT_CG",length(atcg)), 
		  rep("AT_GC",length(atgc)), 
		  rep("AT_TA",length(atta)), 
          rep("GC_AT",length(gcat)), 
		  rep("GC_CG",length(gccg)), 
		  rep("GC_TA",length(gcta)))
		  
oe1bp<-data.frame(exp, Category=catrep, BIN=rep(1:length(atcg),6))
oe1bp<-merge(oe1bp,ct.ord, by=c("Category","BIN"))
names(oe1bp)[4]<-"obs"

oe1bp$odds<-oe1bp$obs/oe1bp$exp

# order by OR and add column of ranks
oe1bp.ord<-oe1bp[order(oe1bp$Category, oe1bp$odds),]
oe1bp.ord <- ddply(oe1bp.ord, .(Category), transform, rk=seq_along(Category))
oe1bp.ord$res<-paste0(1,"bp")

##############################################################################
# Compare 3bp/5bp resolution
##############################################################################
oe5 <- read.table("/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/oe_5bp.txt", sep="\t", header=T, stringsAsFactors=F)
oe3 <- read.table("/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/oe_3bp.txt", sep="\t", header=T, stringsAsFactors=F)
oe.full <- rbind(oe3, oe5)

oe.full2 <- oe.full[(oe.full$odds>0.5 & oe.full$odds<2),]
# oe.full2 <- oe.full[(oe.full$obs>5 & oe.full$odds<5),]
ggplot(oe.full2, aes(x=rk, y=odds, group=res, colour=res))+
	geom_point()+
	facet_wrap(~Category, scales="free")
ggsave("/net/bipolar/jedidiah/mutation/images/3bp_vs_5bp.png")

ggplot(oe.full, aes(factor(res), exp))+
	geom_boxplot()+
	facet_wrap(~Category, scales="free")+
	theme_bw()+
	theme(strip.text.x=element_text(size=14))
ggsave("/net/bipolar/jedidiah/mutation/images/3bp_vs_5bp_box.png", width=13.5, height=8.5)

# include 1bp for comparison
oe.full <- rbind(oe1bp.ord, oe.full)


ggplot(oe.full, aes(x=rk, y=odds, group=res, colour=res))+
	geom_point()+
	facet_wrap(~Category, scales="free")
ggsave("/net/bipolar/jedidiah/mutation/images/1_bp_vs_3bp_vs_5bp.png")

oe.full$logodds <- log(oe.full$odds)
ggplot(oe.full, aes(sample = logodds, colour=res))+
	geom_point(stat = "qq")+
	facet_wrap(~Category, scales="free")
ggsave("/net/bipolar/jedidiah/mutation/images/res_qq.png")

cats <- unique(aggseq$Category)

# test for significant difference in odds between 3bp/5bp across all bins for each category
for(i in 1:6){
	oe3.exp <- oe.full[(oe.full$Category==cats[i] & oe.full$res=="3bp"),8]
	oe5.exp <- oe.full[(oe.full$Category==cats[i] & oe.full$res=="5bp"),8]
	oe1.exp <- oe.full[(oe.full$Category==cats[i] & oe.full$res=="1bp"),8]
	
	test <- t.test(oe3.exp, oe5.exp)
	test2 <- t.test(oe1.exp, oe3.exp)
	
	cat <- cats[i]
	print(cat)
	print(test$p.value)
	print(test2$p.value)
}

# [1] "AT_CG"
# [1] 0.8278162
# [1] "AT_GC"
# [1] 0.9289108
# [1] "AT_TA"
# [1] 0.9737023
# [1] "GC_AT"
# [1] 0.909962
# [1] "GC_CG"
# [1] 0.9254659
# [1] "GC_TA"
# [1] 0.8593596

ggplot(oe.full, aes(factor(res), exp))+
	geom_boxplot()+
	facet_wrap(~Category, scales="free")
ggsave("/net/bipolar/jedidiah/mutation/images/1bp_vs_3bp_vs_5bp_box.png")

# Test for uniformity 
# order by bins
oe2.ord2 <- oe2[order(oe2$BIN),]

for(i in 1:6){
	# aggcat <- oe5[oe5$Category==cats[i],]
	aggcat <- oe2.ord2[oe2.ord2$Category==cats[i],]
	aggcat <- aggcat[aggcat$obs>25,]
	test <- chisq.test(aggcat$obs, p=aggcat$exp/sum(aggcat$exp))
	cat <- cats[i]
	print(cat)
	print(test$p.value)
}

# Test for uniformity of 5bp motifs that share a 3bp motif

b <- c("A", "C", "G", "T")

for(i in 1:6){
	
	aggcat <- aggseq[aggseq$Category==cats[i],]
	aggcat <- aggcat[aggcat$freq>25,]
	# test <- chisq.test(aggcat$obs, p=aggcat$exp/sum(aggcat$exp))
	cat <- cats[i]
	msg <- paste0("checking category ",cat)
	print(msg)
	# print(cat)
	# print(test$p.value)

	# print(head(aggcat))

	for(j in 1:4){
		for(k in 1:4){
			ref <- unique(substr(aggcat$Sequence,3,3))
			motif <- paste0(b[j],ref,b[k])
			dat <- aggcat[substr(aggcat$Sequence,2,4)==motif,]
			
			dat$exp <- dat$COUNT*(sum(dat$freq)/sum(dat$COUNT))
			# print(head(dat))
			cat <- cats[i]

			test <- chisq.test(dat$freq, p=dat$exp/sum(dat$exp))
			if(test$p.value>0.05/96){
				print(cat)
				print(motif)
				print(test$p.value)
				# print(dat)
			}
		}
	}
}

oe2.ord2$diff <- oe2.ord2$obs-oe2.ord2$exp
oe2.ord2$pred <- "5bp"

# Plot residual spatial variation after 5bp estimate
ggplot(oe2.ord2, aes(x=BIN, y=diff))+
	geom_bar(stat="identity", position="identity")+
	facet_wrap(~Category)
ggsave("/net/bipolar/jedidiah/mutation/images/diffs.png")

# [1] "AT_CG"
# [1] 1.060231e-05
# [1] "AT_GC"
# [1] 9.193585e-60
# [1] "AT_TA"
# [1] 8.038046e-54
# [1] "GC_AT"
# [1] 0
# [1] "GC_CG"
# [1] 5.570329e-33
# [1] "GC_TA"
# [1] 4.950446e-271
# [1] "checking category AT_CG"
# [1] "checking category AT_GC"
# [1] "checking category AT_TA"
# [1] "AT_TA"
# [1] "TAC"
# [1] 0.1405975
# [1] "checking category GC_AT"
# [1] "checking category GC_CG"
# [1] "checking category GC_TA"
# [1] "GC_TA"
# [1] "CCG"
# [1] 0.03218077


##############################################################################
# Compare 5bp predictions vs 
# Replication timing, recombination rate, gc content model
##############################################################################
reptime <- read.table("/net/bipolar/jedidiah/mutation/reference_data/lymph_rep_time.txt", header=F, stringsAsFactors=F, sep="\t")
names(reptime) <- c("CHR", "POS", "TIME")
reptime$BIN <- ceiling(reptime$POS/binw)
rtagg <- aggregate(TIME~CHR+BIN, reptime, mean)

# Plot replication timing by chromosome
# ggplot(rtagg, aes(x=BIN, y=TIME))+
	# geom_point()+
	# facet_wrap(~CHR)
# ggsave("/net/bipolar/jedidiah/mutation/images/rep_time.png")

rcrate <- read.table("/net/bipolar/jedidiah/mutation/reference_data/recomb_rate.bed", header=T, stringsAsFactors=F, sep="\t")
rcrate$CHR <- as.numeric(substring(rcrate$CHR, 4))
rcrate$POS <- (rcrate$START+rcrate$END)/2
rcrate$BIN <- ceiling(rcrate$POS/binw)
rcagg <- aggregate(RATE~CHR+BIN, rcrate, mean)

rtagg20 <- rtagg[rtagg$CHR==2,]
rcagg20 <- rcagg[rcagg$CHR==2,]
gc20 <- bins[,4:5]

mut.cov <- merge(ct.ord, rcagg20, by="BIN")
mut.cov <- merge(mut.cov, rtagg20, by="BIN")
mut.cov <- merge(mut.cov, gc20, by="BIN")
mut.cov$prop_GC <- mut.cov$prop_GC-0.5



d <- data.frame()
for(i in 1:6){
	cat <- cats[i]
	aggcat <- mut.cov[mut.cov$Category==cats[i],]
	# mut.lm <- lm(freq~TIME+RATE+prop_GC, data=aggcat)
	# mut.lm <- glm(freq~TIME+RATE+prop_GC, family="poisson", data=aggcat)
	mut.lm <- glm.nb(freq~TIME+RATE+prop_GC, data=aggcat)
	# z <- mut.lm$coefficients
	
	fits <- mut.lm$fitted.values
	names(fits) <- aggcat$BIN
	# names(fits) <- as.character(ceiling(as.numeric(names(fits))/6))
	
	df <- data.frame(Category=cat, 
				   BIN=as.numeric(names(fits)), 
				   exp=fits, 
				   obs=aggcat$freq,
				   stringsAsFactors=F)
				   
	df$diff <- df$obs-df$exp
	df$pred <- "model"
	
	d <- rbind(d,df)
	
	z <- summary(mut.lm)
	
	print(cat)
	print(z)
}

mut.diff<-merge(mut.cov, oe2.ord2, by=c("Category", "BIN"))

d.diff <- data.frame()
for(i in 1:6){
	cat <- cats[i]
	aggcat <- mut.diff[mut.diff$Category==cats[i],]

	mut.lm <- lm(diff~TIME+RATE+prop_GC, data=aggcat)
	
	fits <- mut.lm$fitted.values
	names(fits) <- aggcat$BIN
	# names(fits) <- as.character(ceiling(as.numeric(names(fits))/6))
	
	df <- data.frame(Category=cat, 
				   BIN=as.numeric(names(fits)), 
				   exp=fits, 
				   obs=aggcat$diff,
				   stringsAsFactors=F)
				   
	df$diff <- df$obs-df$exp
	df$pred <- "model"
	
	d.diff <- rbind(d.diff,df)
	
	z <- summary(mut.lm)
	
	print(cat)
	print(z)
}	

diffm<-melt(d.diff[,1:4], id.vars=c("Category", "BIN"), value.var=c("exp", "obs"))
# diffm$value<-as.numeric(diffm$value)

# Plot residual spatial variation after 5bp estimate
ggplot(diffm, aes(x=factor(variable), y=value))+
	# geom_bar(stat="identity", position="identity")+
	geom_boxplot()+
	facet_wrap(~Category, scales="free")+
	scale_x_discrete(labels=c("exp"="5bp+features", "obs"="5bp only"))+
	ylab("Observed-Expected")+
	xlab(NULL)+
	theme_bw()+
	theme(strip.text.x=element_text(size=14))
ggsave("/net/bipolar/jedidiah/mutation/images/hier_diffs_box.png", width=13.5, height=8.5)

ggplot(diffm, aes(x=BIN, y=value, colour=variable))+
	# geom_bar(stat="identity", position="identity")+
	geom_point(alpha=0.4)+
	scale_colour_manual("Model", values=myPaletteCat(4)[3:4], labels=c("exp"="5bp+features", "obs"="5bp only"))+
	facet_wrap(~Category, scales="free")+
	# scale_x_discrete(labels=c("exp"="5bp+features", "obs"="5bp only"))+
	ylab("Observed-Expected")+
	xlab(NULL)+
	theme_bw()+
	theme(strip.text.x=element_text(size=14))
ggsave("/net/bipolar/jedidiah/mutation/images/hier_diffs_pt.png", width=13.5, height=8.5)


# Points with ideogram track



diffm$BIN2<-diffm$BIN*100000

p <- plotIdeogram(hg19IdeogramCyto, "chr20", xlabel=TRUE, alpha=0, zoom.region=c(min(diffm[diffm$Category=="AT_GC",]$BIN2),max(diffm[diffm$Category=="AT_GC",]$BIN2)))

p2<-ggplot(diffm[diffm$Category=="AT_GC",], aes(x=BIN2, y=value, colour=variable))+
	# geom_bar(stat="identity", position="identity")+
	geom_point(alpha=0.5, size=5)+
	scale_colour_manual("Model", values=myPaletteCat(8)[6:7], labels=c("exp"="5bp+features", "obs"="5bp only"))+
	# facet_wrap(~Category, scales="free", ncol=1)+
	# scale_x_discrete(labels=c("exp"="5bp+features", "obs"="5bp only"))+
	ylab("Error")+
	xlab(NULL)+
	theme_bw()+
	theme(axis.title.y=element_text(size=16), 
	      axis.text.y=element_text(size=14),
		  legend.title=element_text(size=16),
		  legend.text=element_text(size=14),
		  axis.text.x=element_blank(), 
		  axis.ticks.x=element_blank())
	
tracks(p,p2, heights=c(1.5,7))
ggsave("/net/bipolar/jedidiah/mutation/images/hier_diffs_pt2.png", width=13.5, height=9)


oe2.ord3 <- oe2.ord2[,-5]

compare.mod <- rbind(oe2.ord3, d)
compare.mod$odds <- compare.mod$obs/compare.mod$exp

# order by OR and add column of ranks
compare.mod.ord <- compare.mod[order(compare.mod$Category, compare.mod$odds),]
compare.ord <- ddply(compare.mod.ord, .(Category, pred), transform, rk=seq_along(Category))
compare.ord <- compare.ord[compare.ord$obs>10,]

# Plot odds profile comparing 5bp predictions vs model predictions
# ggplot(compare.ord, aes(x=BIN, y=exp, colour=pred))+
	# geom_point(alpha=0.4)+
	# facet_wrap(~Category, scales="free")
# ggsave("/net/bipolar/jedidiah/mutation/images/residual_spatial_variation.png")

compare.ord$logodds <- log(compare.ord$odds)
compare.ord2<-compare.ord[compare.ord$logodds>-1,]
ggplot(compare.ord2, aes(sample = logodds, colour=pred))+
	geom_point(stat = "qq")+
	facet_wrap(~Category, scales="free")
ggsave("/net/bipolar/jedidiah/mutation/images/model_qq.png")

# T-test for difference in odds (or expected counts--use col 3) 
# between 5bp and model predictions
# Mean should be close to 1 for both--need a different test
# for(i in 1:6){
	# mod.exp <- compare.ord[(compare.ord$Category==cats[i] & compare.ord$pred=="model"),3]
	# oe5.exp <- compare.ord[(compare.ord$Category==cats[i] & compare.ord$pred=="5bp"),3]
	
	# test <- t.test(mod.exp, oe5.exp)
	
	# cat <- cats[i]
	# print(cat)
	# print(test$p.value)
# }

compare.mod2 <- melt(compare.mod[,c(1,2,3,4,6)], id.vars=c("Category", "BIN", "pred"), value.var=c("exp", "obs"))
compare.mod3 <- compare.mod2[compare.mod2$variable!="diff",]
compare.mod3$model <- paste0(compare.mod3$pred, compare.mod3$variable)

# compare.mod$log_abs_diff<-log(abs(compare.mod$diff))
# ggplot(compare.mod, aes(x=BIN, y=diff, colour=pred))+
	# geom_point(alpha=0.4)+
	# facet_wrap(~Category, scales="free")
# ggsave("/net/bipolar/jedidiah/mutation/images/residual_spatial_variation.png")
 
names(oe.full)[7]<-"pred"
oe.full<-oe.full[,c(1,2,3,4,7)]
compare.all<-rbind(compare.mod[,c(1,2,3,4,6)], oe.full)
mod.corr<-ddply(compare.all, .(Category, pred), summarize, cor=cor(exp, obs))

ggplot(mod.corr, aes(x=Category, y=cor, group=pred, fill=pred, colour=pred))+
	geom_bar(stat="identity", position="dodge")+
	scale_colour_brewer("Predictor", palette="Dark2")+
	scale_fill_brewer("Predictor", palette="Dark2")+
	ylab("Correlation with observed count")+
	theme_bw()
ggsave("/net/bipolar/jedidiah/mutation/images/5bp_vs_mod.png", width=14, height=7)

##############################################################################
# Compare 10k and 100k window size
##############################################################################
oe5.10 <- read.table("/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/oe_5bp_10k.txt", sep="\t", header=T, stringsAsFactors=F)
oe5$res <- "100k"
oe5.10$res <- "10k"
oe5.10$rk <- oe5.10$rk/10
oe.full10 <- rbind(oe5, oe5.10)
oe.full10 <- oe.full10[(oe.full10$odds>0.5 & oe.full10$odds<3),]

ggplot(oe.full10, aes(x=rk, y=odds, group=res, colour=res))+
	geom_point()+
	scale_colour_brewer("Bin Size", palette="Dark2")+
	facet_wrap(~Category, scales="free")+
	theme_bw()+
	theme(strip.text.x=element_text(size=14))
ggsave("/net/bipolar/jedidiah/mutation/images/100k_vs_10k.png", width=15, height=6)

oe.full10$logodds <- log(oe.full10$odds)
ggplot(oe.full10, aes(sample = logodds, colour=res))+
	geom_point(stat = "qq")+
	facet_wrap(~Category, scales="free")
ggsave("/net/bipolar/jedidiah/mutation/images/100k_vs_10k_qq.png")
