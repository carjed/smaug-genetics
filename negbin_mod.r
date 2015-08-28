##############################################################################
# Negative binomial regression models and plots
##############################################################################
agg_cov <- merge(agg_5bp_100k[,c(1:6)], mut_cov, by=c("CHR", "BIN"))

# Motifs + genomic features
compare.all <- data.frame()
mut_cats <- unique(agg_5bp_100k$Category2)
for(i in 1:length(mut_cats)){
	cat1 <- mut_cats[i]
	aggcat <- agg_cov[agg_cov$Category2==mut_cats[i],]
	
	
	# pc.dat <- prcomp(aggcat[,6:17], center=T, scale=T)
	# scores <- as.data.frame(pc.dat$x)
	# aggcat <- cbind(aggcat, scores)

	feat_mod_formula<-as.formula(paste("obs~", paste(covnames[1:12], collapse="+")))
	# motif_mod_formula<-as.formula(paste("obs~", paste(covnames, collapse="+")))
	full_mod_formula<-as.formula(paste("obs~exp+", paste(covnames[1:12], collapse="+")))
	mut_lm_feat <- glm.nb(feat_mod_formula, data=aggcat)
	mut_lm_motif <- glm.nb(obs~exp, data=aggcat)
	mut_lm_full <- glm.nb(full_mod_formula, data=aggcat)
	
	fits_feat <- mut_lm_feat$fitted.values
	names(fits_feat) <- paste0(aggcat$CHR,".",aggcat$BIN)
	
	fits_motif <- mut_lm_motif$fitted.values
	names(fits_motif) <- paste0(aggcat$CHR,".",aggcat$BIN)
	
	fits_full <- mut_lm_full$fitted.values
	names(fits_full) <- paste0(aggcat$CHR,".",aggcat$BIN)
	
	BIN <- as.integer(gsub(".*\\.", "", names(fits_feat)))
	CHR <- as.integer(gsub("\\..*", "", names(fits_feat)))
	
	model_dat_feat <- data.frame(CHR,Category2=cat1, BIN, 
				   exp=fits_feat,
				   obs=aggcat$obs,
				   stringsAsFactors=F)			   
	model_dat_feat$res <- "features"
	
	model_dat_motif <- data.frame(CHR,Category2=cat1, BIN, 
				   exp=fits_motif, 
				   # exp=aggcat$exp,
				   obs=aggcat$obs,
				   stringsAsFactors=F)		   
	model_dat_motif$res <- "5bp"
	
	model_dat_full <- data.frame(CHR,Category2=cat1, BIN, 
				   exp=fits_full, 
				   obs=aggcat$obs,
				   stringsAsFactors=F)			   
	model_dat_full$res <- "5bp+features"
	
	compare.all <- rbind(compare.all,model_dat_feat,model_dat_motif,model_dat_full)
	
	# z <- summary(mut_lm_full)
	# print(cat1)
	# print(z)
}

# Update comparison data frame
# compare.all<-rbind(agg_5bp_100k[,c(1,2,3,4,5,8)], d, d1) #<- this version uses direct estimates from motif rates; not fair comparison
compare.all$res <- factor(compare.all$res, levels = c("features", "5bp", "5bp+features"))
compare.all$diff <- compare.all$obs-compare.all$exp
compare.all$Category2 <- factor(compare.all$Category2)
levels(compare.all$Category2) <- c("AT>CG", "AT>GC", "AT>TA", "(CpG)GC>AT", "(CpG)GC>CG", "(CpG)GC>TA", "GC>AT", "GC>CG", "GC>TA")

ca.gp <- group_by(compare.all, Category2, res)
compare.all<-mutate(ca.gp, dm=mean(diff), dsd=sd(diff), zscore=(diff-dm)/dsd)
compare.all$err<-abs(compare.all$diff)/compare.all$obs

# get mean absolute error for each category/model
c.gp<-group_by(data.frame(compare.all)[compare.all$obs>10,], Category2, res)
ca.gp <- group_by(compare.all, Category2, res)
ca<-summarise(c.gp, mean=mean(err), sd=std(err), meanobs=mean(obs))
caout<-paste0(parentdir, "/output/", nbp, "bp_err.txt")
write.table(ca, caout, col.names=T, row.names=F, quote=F, sep="\t")
cad<-dcast(data.frame(ca), Category2~res, value.var="mean")



	
# compare.all<-merge(compare.all, ca.sum, by=c("Category2", "res"))
# compare.all$zscore<-(compare.all$diff-compare.all$mean)/compare.all$sd

# Function to compute standard error for correlations
corSE<-function(corval, ct){
	sqrt((1-corval^2)/(ct-2))
}

p2 <- ggplot(agg_5bp_100k[agg_5bp_100k$res=="5bp",], aes(x=obs, y=nmotifs))+
	# geom_bar(stat="identity", position="identity")+
	geom_point(alpha=0.2, size=3)+
	# geom_abline(intercept=0, slope=1)+
	scale_colour_manual("Model", values=myPaletteCat(8)[6:7])+
	facet_wrap(~Category2, scales="free", ncol=3)+
	# scale_x_discrete(labels=c("exp"="5bp+features", "obs"="5bp only"))+
	ylab("Mutable motifs")+
	xlab("Observed count")+
	theme_bw()+
	theme(axis.title.y=element_text(size=16), 
	      axis.text.y=element_text(size=14),
		  axis.title.x=element_text(size=16), 
	      axis.text.x=element_text(size=14),
		  legend.title=element_text(size=16),
		  legend.text=element_text(size=14))

hierfile8<-paste0(parentdir, "/images/hier_diffs_pt8.png")
ggsave(hierfile8, width=18, height=18)

a2<-mutate(agg_5bp_100k, rel=obs/nmotifs)
ggplot(a2, aes(x=rel, y=nmotifs))+geom_point()+facet_wrap(~Category2, scales="free")
ggsave("/net/bipolar/jedidiah/mutation/images/rel_ct_cor.png")

a3<-a2 %>% group_by(Category2) %>% summarise(cor=cor(rel, nmotifs))

plotdat<-compare.all[abs(compare.all$diff)<300 & compare.all$res=="5bp",]
p2 <- ggplot(plotdat, aes(x=obs, y=exp, colour=res))+
	# geom_bar(stat="identity", position="identity")+
	geom_point(alpha=0.2, size=3)+
	# geom_point(alpha=0.2, size=3, data=filter(plotdat, res=="features"))+
	# geom_point(alpha=0.2, size=3, data=filter(plotdat, res=="5bp+features"))+
	geom_abline(intercept=0, slope=1)+
	scale_colour_manual("Model", values=myPaletteCat(8)[7])+
	facet_wrap(~Category2, ncol=3)+
	# scale_x_discrete(labels=c("exp"="5bp+features", "obs"="5bp only"))+
	ylab("Predicted count")+
	xlab("Observed count")+
	theme_bw()+
	theme(axis.title.y=element_text(size=16), 
	      axis.text.y=element_text(size=14),
		  axis.title.x=element_text(size=16), 
	      axis.text.x=element_text(size=14),
		  legend.title=element_text(size=16),
		  legend.text=element_text(size=14))

hierfile7<-paste0(parentdir, "/images/hier_diffs_pt7.png")
ggsave(hierfile7, width=18, height=18)

plotdat<-compare.all[compare.all$CHR==2 & compare.all$diff>-150,]
# Plot scatterplot across chr of model errors (using negbin motif predictions)
p2 <- ggplot(plotdat, aes(x=BIN, y=zscore, colour=res))+
	# geom_bar(stat="identity", position="identity")+
	geom_point(alpha=0.4, size=4)+
	geom_point(alpha=0.4, size=4, data=filter(plotdat, res=="features"))+
	geom_point(alpha=0.4, size=4, data=filter(plotdat, res=="5bp+features"))+
	scale_colour_manual("Model", values=myPaletteCat(8)[6:8])+
	facet_wrap(~Category2, scales="free", ncol=3)+
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

hierfile5<-paste0(parentdir, "/images/hier_diffs_pt5.png")
ggsave(hierfile5, width=18, height=18)

# Summarize predictions by average pct error
# compare.all$pct<-abs(compare.all$diff)/compare.all$obs
# mod.corr2 <- ddply(compare.all, .(Category2, res), summarize, meanerr=mean(pct))

# QQ plot comparison of 3 models
qqplot.data <- function(dat){
	vec<-dat$obs
	# following four lines from base R's qqline()
	y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
	x <- qnorm(c(0.25, 0.75))
	slope <- diff(y)/diff(x)
	int <- y[1L] - slope * x[1L]

	d <- data.frame(resids = vec)

	ggplot(dat, aes(sample = exp, colour=res))+
		stat_qq(alpha=0.5)+
		geom_abline(slope = slope, intercept = int)+
		coord_cartesian(ylim=c(0,500))+
		theme_bw()
}

# qqplot.data(caqq)
# ggsave("/net/bipolar/jedidiah/mutation/images/model_qq_full.png")

# Plot barcharts comparing obs/exp correlation for different models
mod.gp <- group_by(compare.all[abs(compare.all$diff)<300,], Category2, res)
mod.corr <- summarise(mod.gp, num=length(exp), cor=cor(exp, obs, method="pearson"))
mod.corr$SE <- corSE(mod.corr$cor, mod.corr$num)
limits <- aes(ymax = mod.corr$cor + mod.corr$SE, ymin=mod.corr$cor - mod.corr$SE)
dodge <- position_dodge(width=0.9)

ggplot(mod.corr, aes(x=Category2, y=cor, fill=res))+
	geom_bar(stat="identity", position=dodge)+
	scale_colour_manual("Predictor", values=myPaletteCat(8)[6:8])+
	scale_fill_manual("Predictor", values=myPaletteCat(8)[6:8])+
	xlab("Category")+
	ylab("Correlation with observed count")+
	geom_errorbar(limits, position=dodge, width=0.25)+
	theme_bw()+
	theme(legend.title = element_text(size=18),
		  legend.text = element_text(size=16),
		  axis.title.x = element_text(size=20),
		  axis.title.y = element_text(size=20),
		  axis.text.y = element_text(size=16), 
		  axis.text.x = element_text(size=16, angle = 45,  vjust=1, hjust=1.01))

modelbar<-paste0(parentdir, "/images/gw_5bp_vs_mod.png")		  
ggsave(modelbar, width=7, height=7)

# old code--model motif error as function of features
mut.diff<-merge(agg_5bp_100k, mut_cov, by=c("CHR", "BIN"))
d.diff <- data.frame()
for(i in 1:length(mut_cats)){
	cat1 <- mut_cats[i]
	aggcat <- mut.diff[mut.diff$Category2==cat1,]
	
	negbin_mod_formula<-as.formula(paste("diff~", paste(covnames[1:12], collapse="+")))
	mut.lm <- lm(negbin_mod_formula, data=aggcat)
	
	fits <- mut.lm$fitted.values
	names(fits) <- paste0(aggcat$CHR,".",aggcat$BIN)
	
	BIN <- as.integer(gsub(".*\\.", "", names(fits)))
	CHR <- as.integer(gsub("\\..*", "", names(fits)))
	
	df <- data.frame(CHR, Category2=cat1, BIN,
				   exp=fits, 
				   obs=aggcat$diff,
				   stringsAsFactors=F)
				   
	df$diff <- df$obs-df$exp
	df$pred <- "model"
	
	d.diff <- rbind(d.diff,df)
	
	z <- summary(mut.lm)
	
	print(cat1)
	print(z)
}	

##############################################################################
# Plot error from 5bp motif model predictions vs 5bp+features model predictions
# with ideogram track
##############################################################################
diffm <- melt(d.diff[,1:5], id.vars=c("CHR", "Category2", "BIN"), value.var=c("exp", "obs"))
diffm <- diffm[diffm$CHR==2,]
diffm$BIN2 <- diffm$BIN*100000

p <- plotIdeogram(hg19IdeogramCyto, "chr2", xlabel=TRUE, alpha=0, 
				  zoom.region=c(min(diffm[diffm$Category=="AT_GC",]$BIN2),max(diffm[diffm$Category=="AT_GC",]$BIN2)))

# p2 <- ggplot(diffm[diffm$Category=="GC_AT",], aes(x=BIN2, y=value, colour=variable))+
# p2 <- ggplot(diffm, aes(x=BIN, y=value, colour=variable))+
	# geom_point(alpha=0.2, size=3)+
	# scale_colour_manual("Model", values=myPaletteCat(8)[6:7], labels=c("exp"="5bp+features", "obs"="5bp only"))+
	# facet_wrap(~Category2, scales="free", ncol=1)+
	# scale_x_continuous(breaks=pretty_breaks(n=25))+
	# ylab("Error")+
	# xlab(NULL)+
	# theme_bw()+
	# theme(axis.title.y=element_text(size=16), 
	      # axis.text.y=element_text(size=14),
		  # legend.title=element_text(size=16),
		  # legend.text=element_text(size=14))
		  
		  # ,
		  # axis.text.x=element_blank(), 
		  # axis.ticks.x=element_blank()

c2b<-diffm[diffm$Category=="AT_GC" & diffm$variable=="obs",]
fixed(p)<-FALSE
p2 <- ggplot(c2b, aes(x=BIN, y=value, colour=variable))+
	geom_point(alpha=0.4, size=4)+
	scale_colour_manual("Model", values=myPaletteCat(8)[7], labels=c("obs"="5bp"))+
	# facet_wrap(~Category2, scales="free", ncol=1)+
	# scale_x_discrete(labels=c("exp"="5bp+features", "obs"="5bp only"))+
	ylab("Error")+
	xlab(NULL)+
	theme_bw()+
	theme(axis.title.y=element_text(size=16), 
	      axis.text.y=element_text(size=14),
		  legend.title=element_text(size=16),
		  legend.text=element_text(size=14))

labeled(p2)<-FALSE
labeled(p)<-FALSE
tracks(p2,p, heights=c(6,2))		  
hierfile6<-paste0(parentdir, "/images/hier_diffs_pt6.png")
ggsave(hierfile6, width=16, height=8)		  

c2b<-diffm[diffm$Category2=="AT_GC",]

c2b.gp <- group_by(c2b, variable)
c2b2 <- mutate(c2b.gp, mean=mean(value), sd=sd(value))

# c2b2<-merge(c2b, c2b1, by="variable")
c2b2$zscore<-(c2b2$value-c2b2$mean)/c2b2$sd
# c2b2$zscore<-(c2b2$value)/24.3

p2 <- ggplot(c2b2, aes(x=BIN, y=zscore, colour=variable))+
	geom_point(alpha=0.4, size=4)+
	scale_colour_manual("Model", 
						values=myPaletteCat(8)[6:7], 
						labels=c("obs"="5bp", "exp"="5bp+features"))+
	scale_y_continuous(breaks=seq(-10,10,by=1))+
	# facet_wrap(~Category2, scales="free", ncol=1)+
	# scale_x_discrete(labels=c("exp"="5bp+features", "obs"="5bp only"))+
	ylab("z score")+
	xlab(NULL)+
	geom_hline(yintercept=2)+
	geom_hline(yintercept=-2)+
	theme_bw()+
	theme(axis.title.y=element_text(size=16), 
	      axis.text.y=element_text(size=14),
		  legend.title=element_text(size=16),
		  legend.text=element_text(size=14))

# tracks(p2,p, heights=c(6,2))	
		  
hierfile3<-paste0(parentdir, "/images/hier_diffs_pt3.png")
ggsave(hierfile3, width=16, height=8)