##############################################################################
# Create figures
##############################################################################
nrsq <- expression(
	paste("Fraction of variance explained (Nagelkerke ", R^2, ")", sep=""))

# Define common theme for Nagelkerke R^2 plots
theme_nrsq <- function(base_size = 12, base_family = ""){
  theme_bw(base_size = base_size) %+replace%
  theme(legend.text=element_text(size=14),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=16, angle=90),
    axis.text.x=element_blank(),
    axis.text.y=element_text(size=14))
}

rsq_full_dat <- combineddat %>%
  # filter(mod %in% mod[1:5])
	filter(mod %in% mod[1:6])

rsq_full_dat <- combineddat[1:6,]

Lv7v5v3_full <- ggplot(rsq_full_dat, aes(x=category, y=rsq, fill=mod))+
  geom_bar(stat="identity", position="dodge")+
	# scale_fill_manual("Model", values=c(iwhPalette[c(1,3,4,5,9)]))+
	scale_fill_manual("Model", values=c(iwhPalette[c(1,3,4,5,9,10)]))+
	theme_nrsq()+
	ylab(nrsq)
# ggsave(paste0(parentdir, "/images/Lv7v5v3_rsq_full.png"), width=4, height=6)


fd2 <- fulldat %>%
	filter(mod %in% paste0("logmod", c(3,5,7,"L")))

fd3 <- fd2 %>%
	dplyr::select(-aic) %>%
	spread(mod, rsq)

fd3$fold<-fd3[,5]/fd3[,4]
fd3$fold2<-fd3[,5]/fd3[,3]

buildSegmentData <- function(){
	ngroups <- 9
	nmodels <- 4
	nbars <- ngroups*nmodels

	st_once <- c(seq(1, nbars, nmodels), seq(nmodels, nbars, nmodels))
	st_twice <- setdiff(1:nbars, st_once)
	st_selection <- sort(c(st_once, rep(st_twice, each=2)))

	end_none <- seq(1, nbars, nmodels)
	end_twice <- setdiff(1:nbars, end_none)
	end_selection <- sort(c(rep(end_twice, each=2)))

	vertst <- sort(c(rep(c(1:ngroups-1/9), each=2),
		rep(c(1:ngroups+1/9), each=2),
		c(1:ngroups+1/3),
		c(1:ngroups-1/3)))
	vertend <- vertst

	horizst <- sort(c(
		c(1:ngroups-1/3),
		c(1:ngroups-1/9),
		c(1:ngroups+1/9)))
	horizend <- sort(c(
		c(1:ngroups-1/9),
		c(1:ngroups+1/9),
		c(1:ngroups+1/3)))

	# Build data for geom_segment call in plot
	corplot<-data.frame(
		xst=vertst,
	  xend=vertend,
	  yst=fd2$rsq[st_selection]+.001,
	  yend=fd2$rsq[end_selection]+.01)

	corplot2 <- data.frame(
		xst=horizst,
	  xend=horizend,
	  yst=fd2$rsq[unique(end_selection)]+.01,
	  yend=fd2$rsq[unique(end_selection)]+.01,
	  pval=signif(lrtestdat$pvals, 2))

	corplot2$pval <- ifelse(corplot2$pval<0.001, "***",
		ifelse(corplot2$pval<0.01, "**", ifelse(corplot2$pval<0.05, "*", " ")))

	out <- list()
	out$corplot <- corplot
	out$corplot2 <- corplot2
  return(out)
}

segment_data <- buildSegmentData()
corplot <- segment_data$corplot
corplot2 <- segment_data$corplot2

# Plot pseudo-r^2 for 7-mers vs 5-mers vs 3-mers
rsqdat <- fd2 %>%
  # fulldat %>%
  # filter(mod %in% c("7-mers", "5-mers", "3-mers", "7-mers+features")) %>%
  filter(group=="FULL") %>%
  mutate(Category =
      factor(plyr::mapvalues(category, orderedcats2, orderedcats2),
				levels=orderedcats2))

Lv7v5v3 <- ggplot(rsqdat)+
  geom_bar(aes(x=Category, y=rsq, fill=mod), stat="identity", position="dodge")+
  # scale_fill_manual("Model", values=cbbPalette[c(4,6,7,8)])+
	scale_fill_manual("Model", values=c(iwhPalette[c(3,4,5,9)]))+
  geom_segment(data=corplot, aes(x=xst, xend=xend, y=yst, yend=yend))+
  geom_segment(data=corplot2, aes(x=xst, xend=xend, y=yst, yend=yend))+
  geom_text(data=corplot2,
		aes(x=xst+1/9, y=yst+.005, label=pval, angle=90), size=6)+
  scale_y_continuous(limits=c(0,0.1))+
  theme_bw()+
  theme(legend.text=element_text(size=14),
    axis.title.x=element_text(size=16),
    axis.title.y=element_text(size=16),
    axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
    axis.text.y=element_text(size=14))+
  xlab("Mutation Type")+
  ylab(nrsq)
# ggsave(paste0(parentdir, "/images/Lv7v5v3_rsq.png"), width=12, height=6)

Lv7v5v3_full <- ggplotGrob(Lv7v5v3_full)
Lv7v5v3 <- ggplotGrob(Lv7v5v3)

g <- arrangeGrob(Lv7v5v3_full, Lv7v5v3, nrow=1, widths=c(1,2))
ggsave(paste0(parentdir, "/images/Lv7v5v3_rsq_combined.svg"),
	width=12, height=6, g)

# Plot pseudo-r^2 for 7-mers vs logit
rsqdatL <- fulldat %>%
  filter(mod=="7-mers" | mod=="7-mers+features") %>%
  filter(group=="FULL") %>%
  mutate(Category =
      factor(plyr::mapvalues(category, orderedcats2, orderedcats2),
				levels=orderedcats2))

ggplot(rsqdatL)+
  geom_bar(aes(x=Category, y=rsq, fill=mod),
		stat="identity", position="dodge")+
  scale_fill_manual("Model", values=cbbPalette[c(7,8)])+
	theme_nrsq()+
	ylab(nrsq)+
  xlab("Mutation Type")
ggsave(paste0(parentdir, "/images/7v5v3_rsqL.png"), width=12, height=6)

# Plot pseudo-r^2 for ERVs vs Common
rsqdatEC <- fulldat %>%
	filter(mod=="AV" | mod=="Common" | mod=="ERVs") %>%
  filter(group=="FULL") %>%
  mutate(Category =
      factor(plyr::mapvalues(category, orderedcats2, orderedcats2),
				levels=orderedcats2))
ggplot(rsqdatEC)+
	geom_bar(aes(x=Category, y=rsq, fill=mod),
		stat="identity", position="dodge")+
	scale_fill_manual("Model", values=c("grey30", cbbPalette[c(3,7)]))+
	theme_bw()+
	theme(legend.text=element_text(size=14),
		axis.title.x=element_text(size=16),
		axis.title.y=element_text(size=16),
		axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
		axis.text.y=element_text(size=14))+
	xlab("Mutation Type")+
  ylab(nrsq)
ggsave(paste0(parentdir, "/images/AvCvS_rsq.png"), width=12, height=6)

oldmodnames <- unique(fulldat$mod)

newmodnames <- c("3-mers", "5-mers", "7-mers", "7-mers+features",
	"7-mers (downsampled BRIDGES ERVs)",
	"7-mers (BRIDGES MAC10+ variants)",
	"7-mers (1KG EUR intergenic variants)")

newmodnamesord <- c("3-mers", "5-mers", "7-mers",
	"7-mers (downsampled BRIDGES ERVs)",
	"7-mers (BRIDGES MAC10+ variants)",
	"7-mers (1KG EUR intergenic variants)",
	"7-mers+features")

rsqdatFULL <- fulldat %>%
  filter(group=="FULL") %>%
  mutate(Category =
      factor(plyr::mapvalues(category, orderedcats2, orderedcats2),
				levels=orderedcats2)) %>%
	mutate(mod=
		factor(plyr::mapvalues(mod, oldmodnames, newmodnames),
			levels=newmodnamesord))


alldat <- combineddat[2:8,] %>%
	mutate(mod=factor(plyr::mapvalues(mod, oldmodnames, newmodnames),
		levels=newmodnamesord))

all_full <- ggplot(alldat)+
  geom_bar(aes(x=category, y=rsq, fill=mod),
		stat="identity", position="dodge")+
	scale_fill_manual("Model", values=c(iwhPalette[c(3:9)]))+
  theme_nrsq(legend.position="none")+
  ylab(nrsq)

all <- ggplot(rsqdatFULL)+
  geom_bar(aes(x=Category, y=rsq, fill=mod),
		stat="identity", position="dodge")+
	scale_fill_manual("Model", values=c(iwhPalette[c(3:9)]))+
  theme_bw()+
  theme(
      legend.text=element_text(size=14),
			legend.position=c(0.4,0.7),
			legend.background = element_rect(colour = "black"),
      axis.title.x=element_text(size=16),
      axis.title.y=element_text(size=16),
    axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
      axis.text.y=element_text(size=14))+
  xlab("Mutation Type")+
  ylab(nrsq)

all_full <- ggplotGrob(all_full)
all <- ggplotGrob(all)
g <- arrangeGrob(all_full, all, nrow=1, widths=c(1,2))
ggsave(paste0(parentdir, "/images/all_rsq_combined.png"), width=12, height=6, g)

evcdat <- rsqdatFULL %>%
	filter(grepl("BRIDGES", mod))
EvC <- ggplot(evcdat)+
  geom_bar(aes(x=Category, y=rsq, fill=mod),
		stat="identity", position="dodge")+
  # scale_fill_manual("Model", values=brewer.pal(8, "Set3")[5:6])+
	scale_fill_manual("Model",
		values=c(iwhPalette[c(6,7)]),
		guide = guide_legend(nrow=2))+
  theme_bw()+
  theme(
      legend.text=element_text(size=14),
			legend.position="bottom",
      axis.title.x=element_text(size=16),
      axis.title.y=element_text(size=16),
    axis.text.x=element_text(size=14, angle=45, hjust=1, vjust=1),
      axis.text.y=element_text(size=14))+
  xlab("Mutation Type")+
  ylab(nrsq)
# ggsave(paste0(parentdir, "/images/EvC_rsq.png"), width=12, height=6)

EvC_full_dat <- combineddat %>%
  # filter(mod=="AV" | mod=="Common" | mod=="ERVs") %>%
	filter(mod=="Common" | mod=="ERVs") %>%
	mutate(mod=plyr::mapvalues(mod, c("Common", "ERVs"),
		c("7-mers (BRIDGES MAC10+ variants)", "7-mers (downsampled BRIDGES ERVs)"))) %>%
# combineddat %>%
  # filter(mod=="AV" | mod=="Common" | mod=="ERVs") %>%
	# filter(mod=="Common" | mod=="ERVs") %>%
	mutate(mod=factor(plyr::mapvalues(mod, oldmodnames, newmodnames),
		levels=newmodnamesord)) %>%
	filter(grepl("BRIDGES", mod))

EvC_full <- ggplot(EvC_full_dat)+
  geom_bar(aes(x=category, y=rsq, fill=mod), stat="identity", position="dodge")+
  # scale_fill_manual("Model", values=brewer.pal(8, "Set3")[5:6])+
	scale_fill_manual("Model", values=c(iwhPalette[c(6,7)]))+
	theme_nrsq()+
  ylab(nrsq)
# ggsave(paste0(parentdir, "/images/EvC_rsq_full.png"), width=8, height=6)

EvC_full <- ggplotGrob(EvC_full)
EvC <- ggplotGrob(EvC)
g <- arrangeGrob(EvC_full, EvC, nrow=1, widths=c(1,2))
ggsave(paste0(parentdir, "/images/EvC_rsq_combined.png"), width=12, height=8, g)
