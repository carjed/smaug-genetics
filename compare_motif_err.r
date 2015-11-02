##############################################################################
# Simple script used to plot barchart comparing absolute error from the 
# following negative binomial models:
# -5bp only
# -3bp only
# -genomic features only
# -5bp + features
##############################################################################

d1<-read.table("/net/bipolar/jedidiah/mutation/output/3bp_err.txt", header=T)
d2<-read.table("/net/bipolar/jedidiah/mutation/output/5bp_err.txt", header=T)

d1a<-d1[d1$res=="5bp",]
d1a$res<-"3bp"

d3<-rbind(d2, d1a)

limits <- aes(ymax = d3$mean + d3$sd, ymin=d3$mean-d3$sd)
dodge <- position_dodge(width=0.9)

ggplot(d3, aes(x=Category2, y=mean, fill=res))+
	geom_bar(stat="identity", position=dodge)+
	geom_errorbar(limits, position=dodge, width=0.25)+
	scale_colour_brewer("Model", palette="Dark2")+
	xlab("Category")+
	ylab("Mean absolute error")+
	theme_bw()+
	theme(legend.title = element_text(size=18),
		  legend.text = element_text(size=16),
		  axis.title.x = element_text(size=20),
		  axis.title.y = element_text(size=20),
		  axis.text.y = element_text(size=16), 
		  axis.text.x = element_text(size=16, angle = 45,  vjust=1, hjust=1.01))

ggsave("/net/bipolar/jedidiah/mutation/images/mod_comp.png")