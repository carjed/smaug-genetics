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

	feat_mod_formula<-as.formula(paste("obs~", paste(covnames[1:12], collapse="+")))
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

# Function to compute standard error for correlations
corSE<-function(corval, ct){
	sqrt((1-corval^2)/(ct-2))
}

# plot raw 5bp motif predictions vs observed
p2 <- ggplot(agg_5bp_100k[agg_5bp_100k$res=="5bp",], aes(x=obs, y=nmotifs))+
	geom_point(alpha=0.2, size=3)+
	scale_colour_manual("Model", values=myPaletteCat(8)[6:7])+
	facet_wrap(~Category2, scales="free", ncol=3)+
	ylab("Mutable motifs")+
	xlab("Observed count")+
	theme_bw()+
	theme(axis.title.y=element_text(size=16), 
	      axis.text.y=element_text(size=14),
		  axis.title.x=element_text(size=16), 
	      axis.text.x=element_text(size=14),
		  legend.title=element_text(size=16),
		  legend.text=element_text(size=14))

hierfile8<-paste0(parentdir, "/images/raw_motif_pred_vs_obs.png")
ggsave(hierfile8, width=18, height=18)

a2<-mutate(agg_5bp_100k, rel=obs/nmotifs)
ggplot(a2, aes(x=rel, y=nmotifs))+
	geom_point()+
	facet_wrap(~Category2, scales="free")
ggsave("/net/bipolar/jedidiah/mutation/images/rel_ct_cor.png")

a3<-a2 %>% group_by(Category2) %>% summarise(cor=cor(rel, nmotifs))

# plot negbin model 5bp motif predictions vs observed
plotdat<-compare.all[compare.all$res=="5bp",]
p2 <- ggplot(plotdat, aes(x=obs, y=exp, colour=res))+
	geom_point(alpha=0.2, size=3)+
	scale_colour_manual("Model", values=myPaletteCat(8)[7])+
	facet_wrap(~Category2, ncol=3, scales="free")+
	ylab("Predicted count")+
	xlab("Observed count")+
	theme_bw()+
	theme(axis.title.y=element_text(size=16), 
	      axis.text.y=element_text(size=14),
		  axis.title.x=element_text(size=16), 
	      axis.text.x=element_text(size=14),
		  legend.title=element_text(size=16),
		  legend.text=element_text(size=14))

hierfile7<-paste0(parentdir, "/images/negbin_motif_pred_vs_obs.png")
ggsave(hierfile7, width=18, height=18)

# Plot scatterplot across chr of model errors (using negbin motif predictions)
plotdat<-compare.all[compare.all$CHR==2 & compare.all$diff>-150,]
p2 <- ggplot(plotdat, aes(x=BIN, y=zscore, colour=res))+
	geom_point(alpha=0.4, size=4)+
	geom_point(alpha=0.4, size=4, data=filter(plotdat, res=="features"))+
	geom_point(alpha=0.4, size=4, data=filter(plotdat, res=="5bp+features"))+
	scale_colour_manual("Model", values=myPaletteCat(8)[6:8])+
	facet_wrap(~Category2, scales="free", ncol=3)+
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

# Plot barcharts comparing obs/exp correlation for different models
mod.corr <- compare.all %>% 
			group_by(Category2, res) %>%
			summarise(num=length(exp), cor=cor(exp, obs, method="pearson"))
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