# plot negbin model 5bp motif predictions vs observed
p2 <- ggplot(compare.all, aes(x=obs, y=exp, colour=res))+
	geom_point(alpha=0.2, size=3)+
	geom_point(alpha=0.2, size=3, data=filter(plotdat, res=="features"))+
	geom_point(alpha=0.2, size=3, data=filter(plotdat, res=="motifs_ext+features"))+
	scale_colour_manual("Model", values=myPaletteCat(8)[4:8])+
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

hierfile4 <- paste0(parentdir, "/images/negbin_all_pred_vs_obs.png")
ggsave(hierfile4, width=18, height=18)

# Plot scatterplot across chr of model errors (using negbin motif predictions)
plotdat <- compare.all[compare.all$CHR==2 & compare.all$diff>-150,]
p2 <- ggplot(plotdat, aes(x=BIN, y=zscore, colour=res))+
	geom_point(alpha=0.4, size=4)+
	geom_point(alpha=0.4, size=4, data=filter(plotdat, res=="features"))+
	geom_point(alpha=0.4, size=4, data=filter(plotdat, res=="motifs_ext+features"))+
	scale_colour_manual("Model", values=myPaletteCat(8)[4:8])+
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

hierfile5 <- paste0(parentdir, "/images/hier_diffs_pt5.png")
ggsave(hierfile5, width=18, height=18)
