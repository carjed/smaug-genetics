##############################################################################
# Plot K-mer heatmaps
##############################################################################
rrheat2 <- function(dat, f, levels, facetvar, nbp){
	p <- ggplot()+
	# log(v4*10000+1,2)
	# limits=c(min(dat$v4), max(dat$v4))
	geom_tile(data=dat, aes(x=v3, y=v2a, fill=v4))+
	# geom_text(data=dat, aes(x=v2a, y=v3, label=v4a, family="Courier", size=0.1))+
	geom_rect(data=f, size=0.6, colour="grey70",
		aes(xmin=xlo, xmax=xhi, ymin=ylo, ymax=yhi), fill=NA)+
	scale_fill_gradientn("Relative Rate",
		# colours=myPalette((nbp-1)^4),
    colours=myPaletteO(11),
		trans="log10",
		breaks=10^(seq(-3.65,-.84,0.281)),
		labels=signif(10^(seq(-3.65,-.84,0.281)), 2),
		limits=c(0.0002, 0.2))+
	xlab("3' flank")+
	ylab("5' flank")+
  theme_classic()+
	theme(
		# legend.position="none",
		legend.title = element_text(size=18),
		legend.text = element_text(size=16),
		legend.key.height = unit(1.5, "cm"),
		axis.ticks.x = element_blank(),
		axis.ticks.y = element_blank(),
		axis.text.y = element_blank(),
		axis.text.x = element_blank(),
    axis.title.y = element_blank(),
		axis.title.x = element_blank())

	return(p)
}

for(i in 1:3){
  adj <- i
  nbp <- adj*2+1
  rates <- read.table(paste0(parentdir, "/output/", nbp, "bp_1000k_rates.txt"),
		header=T, stringsAsFactors=F)
  rates$v2 <- substr(rates$Sequence,1,adj)
  rates$v2a <- as.character(lapply(as.vector(rates$v2), reverse_chars))
  rates$v2a <- factor(rates$v2a)
  rates$v3 <- substr(rates$Sequence, adj+2, adj*2+1)
  rates$v4 <- rates$rel_prop
  rates$Category <- gsub("cpg_", "", rates$Category2)
  rates$v5 <- gsub("_", ">", rates$Category)
  rates$v5 <- factor(rates$v5)

  nbox <- length(unique(rates$v2a))
  nint <- nbox/(4^(adj-1))
  xhi <- rep(1:(4^(adj-1)),4^(adj-1))*nint+0.5
  xlo <- xhi-nint
  yhi <- rep(1:(4^(adj-1)),each=4^(adj-1))*nint+0.5
  ylo <- yhi-nint
  f <- data.frame(xlo,xhi,ylo,yhi)

  levs_a <- as.character(lapply(as.vector(levels(rates$v2a)), reverse_chars))
  # p1<-rrheat(rates, f, levs_a, "v5", nbp)
	for(j in 1:6){
		categ <- orderedcats[j]
		p1 <- rrheat2(rates[rates$Category==categ,], f, levs_a, "v5", nbp)
		p1a <- p1+theme(legend.position="none")
		 png(paste0(parentdir, "/images/", categ, "_", nbp, "bp_heatmap.png"),
		 	height=5, width=5, units="in", res=300)
		 pushViewport(viewport(width=unit(5, "in"), height=unit(5, "in")))
		 grid.draw(ggplotGrob(p1a))
		 dev.off()
	}
}

# trim whitespace on panels with imagemagick mogrify
trimcmd <- paste0("mogrify -trim ", parentdir, "/images/*_*bp_heatmap.png")
system(trimcmd)

# extract legend
legend <- get_legend(p1)
png(paste0(parentdir, "/images/heatmap_legend.png"),
	height=8, width=3, units="in", res=300)
grid.draw(legend)
dev.off()

##############################################################################
# 1-mer heatmap
##############################################################################
rrheat1 <- function(dat, facetvar){
    p <- ggplot()+
    geom_tile(data=dat, aes(x=v5, y=1, fill=v4))+
    scale_fill_gradientn("Relative Rate\n",
  		colours=myPaletteO(11),
  		trans="log",
  		breaks=c(0.0002, 0.008, 0.15),
  		labels=c(0.0002, 0.008, 0.15),
  		limits=c(0.0002, 0.15))+
    ylab(" ")+
    theme_classic()+
    theme(
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

rates3 <- read.table(paste0(parentdir, "/output/3bp_1000k_rates.txt"), header=T, stringsAsFactors=F)
rates3$v5 <- factor(gsub("cpg_", "", rates3$Category2))
rates3$v5 <- gsub("_", ">", rates3$v5)
rates1 <- rates3 %>%
  group_by(v5) %>%
  summarise(num=sum(num), COUNT=sum(COUNT), v4=num/COUNT)
rrheat1(rates1, "v5")
ggsave(paste0(parentdir, "/images/1bp_heatmap.png"), height=2, width=24)

rates7out <- read.table(paste0(parentdir, "/output/7bp_1000k_rates.txt"),
	header=T, stringsAsFactors=F)
rates7out <- mutate(rates7out, Category=gsub("cpg_", "", Category2))

rates7out <- merge(rates7out, r5m, by=c("Category", "Sequence"), all.x=T)
rates7out <- rates7out %>% dplyr::select(Type=Category, Motif=Sequence,
	nERVs=num, nMotifs=COUNT, ERV_rel_rate=rel_prop.x, nERVs_DS=num.x,
	ERV_DS_rel_rate=rel_prop.y, nMAC10=num.y, MAC10_rel_rate=common_rel_prop)



rates5 <- read.table(paste0(parentdir, "/output/5bp_1000k_rates.txt"),
	header=T, stringsAsFactors=F)
rates5out <- rates5 %>%
	mutate(Category=gsub("cpg_", "", Category2)) %>%
	dplyr::select(Type=Category, Motif=Sequence,
	nERVs=num, nMotifs=COUNT, ERV_rel_rate=rel_prop)

rates3 <- read.table(paste0(parentdir, "/output/3bp_1000k_rates.txt"),
	header=T, stringsAsFactors=F)
rates3out <- rates3 %>%
	mutate(Category=gsub("cpg_", "", Category2)) %>%
	dplyr::select(Type=Category, Motif=Sequence,
	nERVs=num, nMotifs=COUNT, ERV_rel_rate=rel_prop)

rates1out <- rates3out %>%
	group_by(Type) %>%
	summarise(nERVs=sum(nERVs), nMotifs=sum(nMotifs), ERV_rel_rate=nERVs/nMotifs)

write.table(rates7out, paste0(parentdir, "/output/7bp_final_rates.txt"),
	col.names=T, row.names=F, quote=F, sep="\t")

write.table(rates5out, paste0(parentdir, "/output/5bp_final_rates.txt"),
	col.names=T, row.names=F, quote=F, sep="\t")

write.table(rates3out, paste0(parentdir, "/output/3bp_final_rates.txt"),
	col.names=T, row.names=F, quote=F, sep="\t")

write.table(rates1out, paste0(parentdir, "/output/1bp_final_rates.txt"),
	col.names=T, row.names=F, quote=F, sep="\t")
