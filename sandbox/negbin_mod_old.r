# QQ plot function (currently unused)
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