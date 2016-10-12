##############################################################################
# Plot K-mer heatmaps
##############################################################################
rrheat <- function(dat, f, levels, facetvar, nbp){
	p <- ggplot()+
	# log(v4*10000+1,2)
	# limits=c(min(dat$v4), max(dat$v4))
	geom_tile(data=dat, aes(x=v2a, y=v3, fill=v4))+
	# geom_text(data=dat, aes(x=v2a, y=v3, label=v4a, family="Courier", size=0.1))+
	geom_rect(data=f, size=0.8, colour="grey30",
		aes(xmin=xlo, xmax=xhi, ymin=ylo, ymax=yhi), fill=NA)+
	scale_fill_gradientn("Relative Rate",
		# colours=myPalette((nbp-1)^4),
    colours=myPalette(11),
		trans="log10",
		breaks=10^(seq(-3.65,-.84,0.281)),
		labels=round(seq(-3.65,-.84,0.281), 2),
		limits=c(0.0002, 0.2))+
	xlab("5' flank")+
	ylab("3' flank")+
  theme_classic()+
	theme(
		# legend.position="none",
		legend.title = element_text(size=18),
	  legend.text = element_text(size=16),
	  strip.text.x = element_text(size=20),
		axis.text.y = element_blank(),
		axis.text.x = element_blank(),
    axis.title.y = element_blank(),
		axis.title.x = element_blank(),
    legend.key.height = unit(1.5, "cm"))+
	facet_wrap(as.formula(paste("~", facetvar)), ncol=6, scales="free_x")

  # if(nbp==7){
  #   p <- p+
  #     scale_x_discrete(labels=function(x) substr(levels,1,1))+
  #     # theme(axis.text.x=element_text(size=8, vjust=c(-1,1)))
  #     theme(axis.text.x=element_text(size=8, vjust = grid::unit(c(0, 0.75), "points")))
  # }
	return(p)
}

for(i in 1:3){
  adj <- i
  nbp <- adj*2+1
  rates <- read.table(paste0("/net/bipolar/jedidiah/mutation/output/", nbp, "bp_1000k_rates.txt"), header=T, stringsAsFactors=F)
  rates$v2 <- substr(rates$Sequence,1,adj)
  rates$v2a <- as.character(lapply(as.vector(rates$v2), reverse_chars))
  rates$v2a <- factor(rates$v2a)
  rates$v3 <- substr(rates$Sequence, adj+2, adj*2+1)
  rates$v4 <- rates$rel_prop
  rates$Category <- gsub("cpg_", "", rates$Category2)
  rates$Category <- gsub("_", ">", rates$Category)
  rates$v5 <- factor(rates$Category)

  nbox <- length(unique(rates$v2a))
  nint <- nbox/(4^(adj-1))
  xhi <- rep(1:(4^(adj-1)),4^(adj-1))*nint+0.5
  xlo <- xhi-nint
  yhi <- rep(1:(4^(adj-1)),each=4^(adj-1))*nint+0.5
  ylo <- yhi-nint
  f <- data.frame(xlo,xhi,ylo,yhi)

  levs_a <- as.character(lapply(as.vector(levels(rates$v2a)), reverse_chars))
  rrheat(rates, f, levs_a, "v5", nbp)
  ggsave(paste0("/net/bipolar/jedidiah/mutation/images/",nbp,"bp_heatmap.png"), height=4, width=24)
}

rrheat1 <- function(dat, facetvar){
    p <- ggplot()+
    # log(v4*10000+1,2)
    # limits=c(min(dat$v4), max(dat$v4))
    geom_tile(data=dat, aes(x=v5, y=1, fill=v4))+
    # geom_text(data=dat, aes(x=v2a, y=v3, label=v4a, family="Courier", size=0.1))+
    scale_fill_gradientn("Relative Rate\n",
  		colours=myPalette((nbp-1)^4),
  		trans="log",
  		breaks=c(0.0002, 0.008, 0.15),
  		labels=c(0.0002, 0.008, 0.15),
  		limits=c(0.0002, 0.15))+
    ylab(" ")+
    theme_classic()+
    theme(
        # legend.position="none",
        legend.title = element_text(size=18),
        legend.key.size = unit(0.2, "in"),
      legend.text = element_text(size=16),
      strip.text.x = element_text(size=20),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=20),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x= element_blank(),
      axis.ticks.y= element_blank())+
    facet_wrap(as.formula(paste("~", facetvar)), ncol=6, scales="free_x")

    return(p)
}

#map_t1<-map_t %>% group_by(v5) %>% summarise(v4=mean(v4))

rates3 <- read.table("/net/bipolar/jedidiah/mutation/output/3bp_1000k_rates.txt", header=T, stringsAsFactors=F)
rates3$v5 <- factor(gsub("cpg_", "", rates3$Category2))
rates3$v5 <- gsub("_", ">", rates3$v5)
rates1 <- rates3 %>%
  group_by(v5) %>%
  summarise(num=sum(num), COUNT=sum(COUNT), v4=num/COUNT)
rrheat1(rates1, "v5")
ggsave("/net/bipolar/jedidiah/mutation/images/1bp_heatmap.png", height=2, width=24)
