#DEPRECATED
##############################################################################
# Plot correlation comparing per-bin counts
##############################################################################
names(agg_5bp_100k_common)[4]<-"common_obs"
bincts<-merge(agg_5bp_100k, agg_5bp_100k_common, by=c("CHR", "BIN", "Category2"))

bincts %>% group_by(Category2) %>% summarise(cor=cor(obs, common_obs))

ggplot(data=bincts,
		aes(x=obs, y=common_obs, colour=Category2, group=Category2))+
	geom_point(alpha=0.4)+
	geom_smooth(method=lm, se=FALSE, colour="black")+
	scale_colour_manual(values=rbg)+
	xlab("Singletons")+
	ylab("Polymorphisms (MAC>10)")+
	theme_bw()+
	theme(legend.position="none")+
	facet_wrap(~Category2, scales="free")
ggsave(paste0(parentdir, "/images/sing_com_corr.png"), width=6.2, height=4.2)

##############################################################################
# Simulate correlations
##############################################################################
output<-data.frame()
for(i in 1:50){
	cat("running simulation", i, "of 50...\n")
	simdat1<-dat_5bp_100k$summ[sample(nrow(dat_5bp_100k$summ), 24176074),]
	simdat1$gp<-sample(0:1, nrow(simdat1), replace=T)

	simdat<-simdat1 %>%
		group_by(CHR, BIN, Category2, gp) %>%
		summarise(count=n())

	simdat$gp<-ifelse(simdat$gp==0,"gp1", "gp2")
	simdat2<-simdat %>% spread(gp, count)
	simdat2<-simdat2[complete.cases(simdat2),]

	simcor<-simdat2 %>% group_by(Category2) %>% summarise(cor=cor(gp1, gp2))
	output<-rbind(output, simcor)
}

simcor<-output %>%
	group_by(Category2) %>%
	summarise(corsim=mean(cor)) %>%
	mutate(gp="sim")

obscor<-bincts %>%
	group_by(Category2) %>%
	summarise(corobs=cor(obs, common_obs)) %>%
	mutate(gp="obs")

tmp<-r.test(2897, simcor$corsim, obscor$corobs)$z

corcomb<-cbind(merge(simcor, obscor, by="Category2"), tmp) %>%
	mutate(cor.p=-2*pnorm(-abs(tmp), log.p=T))

names(simcor)[2]<-"cor"
names(obscor)[2]<-"cor"
corplot<-rbind(simcor, obscor)


corplot$Category2<-factor(corplot$Category2, levels=orderedcats)
corplot$gp<-as.factor(corplot$gp)
corplot$gp<-relevel(corplot$gp, "sim")
corplot<-corplot %>%
	group_by(Category2) %>%
	mutate(cormax=max(cor)) %>%
	arrange(Category2)

corplot$xst<-seq(0.75,9.25,0.5)
corplot$xend<-rep(1:9,each=2)

corplot2<-merge(corplot, corcomb)
corplot2$cor.p<-round(corplot2$cor.p, 1)

ggplot(corplot2, aes(x=Category2, y=cor, fill=gp))+
	geom_bar(position="dodge", stat="identity")+
	geom_segment(data=corplot2, aes(x=xst, xend=xst, yend=cormax+0.1))+
	geom_segment(data=corplot2, aes(x=xst, xend=xend, y=cormax+0.1, yend=cormax+0.1))+
	geom_text(data=corplot2, aes(y=cormax+0.2, label=cor.p))+
	scale_fill_manual("Correlation", values=myPaletteCat(4)[1:2], labels=c("Simulated", "Observed"))+
	scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1))+
	# geom_errorbar(limits, position="dodge", width=1)+
	theme_classic()+
	theme(legend.title=element_blank(),
		  legend.text=element_text(size=16),
		  axis.title.x=element_text(size=20),
		  axis.title.y=element_text(size=20),
		axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
		  axis.text.y=element_text(size=16))+
	xlab("Category")+
	ylab("Correlation")
ggsave(paste0(parentdir, "/images/sing_com_cor_bars.png"), width=9.6, height=4.8)


##############################################################################
# Test for heteroscedasticity
##############################################################################
orderedcats<-c("AT_CG", "AT_GC", "AT_TA",
	"GC_AT", "GC_CG", "GC_TA",
	"cpg_GC_AT", "cpg_GC_CG", "cpg_GC_TA")

dat<-data.frame()
for(i in orderedcats2){
  r5m1 <- r5m %>% filter(Category2==i)
  r5m2 <- r5m1 %>% filter(num.x>=100 & num.y>=100)
  p1 <- lm(rel_prop~common_rel_prop+num.x+num.y, data=r5m1)
  p2 <- lm(rel_prop~common_rel_prop+num.x+num.y, data=r5m2)
  t1 <- bptest(p1)
  t2 <- bptest(p2)
  fullcor <- cor(r5m1$rel_prop, r5m1$common_rel_prop, method="spearman")
  subcor <- cor(r5m2$rel_prop, r5m2$common_rel_prop, method="spearman")

  fullcor <- cor(r5m1$rel_prop, r5m1$common_rel_prop)
  subcor <- cor(r5m2$rel_prop, r5m2$common_rel_prop)

  newrow1 <- data.frame(Category=i, Data="full", Cor=fullcor, Dim=nrow(r5m1), pval=t1$p.value)
  newrow2 <- data.frame(Category=i, Data="sub", Cor=subcor, Dim=nrow(r5m2), pval=t2$p.value)
  dat<-rbind(dat, newrow1)
  dat<-rbind(dat, newrow2)
}

r5m$Category2<-factor(r5m$Category2, levels=orderedcats)

nsing<-sum(r5m$num.x)
ncommon<-sum(r5m$num.y)
pvals<-rep(0, nrow(r5m))
for(i in 1:nrow(r5m)){
	row<-r5m[i,]

	test<-prop.test(c(ceiling(row$num.x/3.073), row$num.y), c(ceiling(nsing/3.073), ncommon))
	pvals[i]<-test$p.value
}
r5m$pval<-pvals
