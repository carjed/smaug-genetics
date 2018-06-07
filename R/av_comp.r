rates7 <- ratelist[[4]]
# r5m <- merge(rates7, rates7_common, by=c("Type", "Motif"))

avrates <- read.xlsx(paste0(parentdir, "/reference_data/AV_rates.xlsx"),
  sheet=10)
names(avrates) <- c("Motif", "alt", "afr", "asn", "eur", "refrev", "altrev")

avrates$CAT <- paste0(substr(avrates$Motif,4,4), substr(avrates$alt, 4,4))
avrates$Category[avrates$CAT=="AC" | avrates$CAT=="TG"] <- "AT_CG"
avrates$Category[avrates$CAT=="AG" | avrates$CAT=="TC"] <- "AT_GC"
avrates$Category[avrates$CAT=="AT" | avrates$CAT=="TA"] <- "AT_TA"
avrates$Category[avrates$CAT=="GA" | avrates$CAT=="CT"] <- "GC_AT"
avrates$Category[avrates$CAT=="GC" | avrates$CAT=="CG"] <- "GC_CG"
avrates$Category[avrates$CAT=="GT" | avrates$CAT=="CA"] <- "GC_TA"

avrates$Motif <- paste0(avrates$Motif, "(", avrates$refrev, ")")

avrates <- avrates %>%
  dplyr::select(Type=Category, Motif, eur)

r5m <- merge(rates7, avrates, by=c("Type", "Motif"))

r5m$Category2 <- ifelse(substr(r5m$Motif,4,5)=="CG",
  paste0("cpg_",r5m$Type),
  r5m$Type)

r5m <- r5m %>%
  mutate(Category2 = plyr::mapvalues(Category2, orderedcats1, orderedcats2))

r5m$Category2 <- factor(r5m$Category2, levels=orderedcats2)

r5m$prop_diff <- (r5m$eur/(mean(r5m$eur)/mean(r5m$ERV_rel_rate)))/r5m$ERV_rel_rate
r5m$prop_diff4 <- r5m$prop_diff
r5m$prop_diff4[r5m$prop_diff< 0.5] <- 0.5
r5m$prop_diff4[r5m$prop_diff>2] <- 2

r5m$v2 <- substr(r5m$Motif,1,3)
r5m$v2a <- as.character(lapply(as.vector(r5m$v2), reverse_chars))
r5m$v2a <- factor(r5m$v2a)
r5m$v3 <- substr(r5m$Motif,3+2,3*2+1)

r5m2 <- r5m %>%
  mutate(eur2=eur/(mean(r5m$eur)/mean(r5m$ERV_rel_rate)))

r5m3 <- r5m %>%
  dplyr::select(Type, Motif, ERV_rel_rate, eur) %>%
  gather(gp, val, ERV_rel_rate, eur)

##############################################################################
# Get BRIDGES MAC10+ rates
##############################################################################
commonfile <- paste0(parentdir, "/summaries/common.full.summary")
bindir <- paste0(parentdir, "/motif_counts/", nbp, "-mers/full")
common_data <- getData(summfile=commonfile, bindir=bindir)
common_data$aggseq <- get_aggseq(common_data$sites, common_data$mct)

i <- 3
gpdat <- common_data$aggseq %>%
  mutate(Type=gsub("cpg_", "", Category2),
    SEQA=substr(Motif, cbp-i, cbp+i),
    SEQB=substr(Motif, cbp*3-i, cbp*3+i),
    Motif=paste0(SEQA, "(", SEQB, ")")) %>%
    dplyr::select(Type, Motif, nERVs) %>%
    group_by(Type, Motif) %>%
    summarise(nMAC10=sum(nERVs))

r5m <- merge(r5m, gpdat, by=c("Type", "Motif")) %>%
  mutate(MAC10_rel_rate=nMAC10/nMotifs)

set.seed(sum(r5m$nMAC10))
ervs_down <- full_data$sites[sample(nrow(full_data$sites), sum(r5m$nMAC10)),]
ervs_down_aggseq <- get_aggseq(ervs_down, common_data$mct)

i <- 3
gpdat <- ervs_down_aggseq %>%
  mutate(Type=gsub("cpg_", "", Category2),
    SEQA=substr(Motif, cbp-i, cbp+i),
    SEQB=substr(Motif, cbp*3-i, cbp*3+i),
    Motif=paste0(SEQA, "(", SEQB, ")")) %>%
    dplyr::select(Type, Motif, nERVs) %>%
    group_by(Type, Motif) %>%
    summarise(nERVs_down=sum(nERVs))

r5m <- merge(r5m, gpdat, by=c("Type", "Motif")) %>%
  mutate(ERV_down_rel_rate=nERVs_down/nMotifs)

rm2 <- r5m %>% dplyr::select(-Category2)

kspvals <- r5m %>%
  mutate(eur_s=eur*mean(ERV_rel_rate)/mean(eur)) %>%
	group_by(Category2) %>%
	summarise(pval=ks.test(ERV_rel_rate, eur_s)$p.value)

##############################################################################
# Theme & styling for 7-mer rate comparison plots
##############################################################################
theme_mac_comp <- function(base_size = 12, base_family = ""){
  theme_bw(base_size = base_size) %+replace%
  theme(legend.position="bottom",
    # strip.text.x=element_text(size=12),
    # strip.background = element_rect(size=3),
    legend.title=element_text(size=12),
    axis.title.x=element_text(size=14),
    axis.title.y=element_text(size=14, angle=90),
		axis.text.x=element_text(size=12, angle=45, hjust=1, vjust=1),
		axis.text.y=element_text(size=12))
}

format_mac_comp <- list(
  # coord_fixed(),
  scale_colour_manual(values=gp_cols),
  guides(colour = guide_legend(title=NULL,
      nrow=3,
      override.aes = list(alpha=1, shape=16))),
  scale_x_log10(#expand=c(0,0),
    labels=c(0.0001, 0.001, 0.01, 0.1),
    breaks=c(0.0001, 0.001, 0.01, 0.1),
    limits=c(0.00005, 0.22)),
  scale_y_log10(#expand=c(0,0),
    labels=c(0.0001, 0.001, 0.01, 0.1),
		breaks=c(0.0001, 0.001, 0.01, 0.1),
		limits=c(0.00005, 0.22)))

format_mac_comp_facet <- function(corlabs){
  list(
    geom_point(data=rm2, alpha=0.3, size=2, colour="grey70"),
    geom_hex(bins=50),
    geom_abline(intercept=0, linetype="dashed"),
    geom_text(data=corlabs, aes(x=x, y=y, label=label),
              hjust=0, vjust=1, size=4),
    scale_fill_gradientn(colours = myPalette(6)),
    facet_wrap(~Category2, dir="v"),
    theme(legend.position="none")
  )
}

format_mac_comp_nf <- list(
  geom_point(aes(group=Category2, colour=Category2),
    alpha=0.3, size=2, shape=1),
  geom_abline(intercept=0, linetype="dashed"))

ervlab <- "Relative mutation rate \n (BRIDGES ERVs)"
ervmasklab <- "Relative mutation rate \n (BRIDGES ERVs passing 1000G strict mask)"
mac10lab <- "Relative mutation rate \n (BRIDGES MAC10+ variants)"
avlab <- "Relative mutation rate \n (1000G EUR intergenic variants)"

##############################################################################
# ERVs vs AV
##############################################################################
p <- ggplot(data=r5m,
    aes(x=ERV_rel_rate, y=eur*mean(ERV_rel_rate)/mean(eur)))+
  xlab(ervlab)+
  ylab(avlab)+
  format_mac_comp+
  theme_mac_comp()

p + format_mac_comp_nf
ggsave(paste0(parentdir, "/images/ERV_vs_AV_corr.pdf"),
	width=4, height=5)

corlabs <- r5m %>%
	group_by(Category2) %>%
	summarise(cor=cor(ERV_rel_rate, eur, method="spearman")) %>%
	mutate(label=paste0("r=", round(cor,2))) %>%
  mutate(x=0.0001, y=0.1)

p + format_mac_comp_facet(corlabs)
ggsave(paste0(parentdir, "/images/ERV_vs_AV_corr_facet.pdf"),
  width=6, height=6)

##############################################################################
# ERVs vs MAC10
##############################################################################
p <- ggplot(data=r5m,
    aes(x=ERV_rel_rate, y=MAC10_rel_rate*mean(ERV_rel_rate)/mean(MAC10_rel_rate)))+
  xlab(ervlab)+
  ylab(mac10lab)+
  format_mac_comp+
  theme_mac_comp()

p + format_mac_comp_nf
ggsave(paste0(parentdir, "/images/ERV_vs_MAC10_corr.pdf"),
	width=4, height=5)

corlabs <- r5m %>%
	group_by(Category2) %>%
	summarise(cor=cor(ERV_rel_rate, MAC10_rel_rate, method="spearman")) %>%
  mutate(label=paste0("r=", round(cor,2))) %>%
  mutate(x=0.0001, y=0.1)

p + format_mac_comp_facet(corlabs)
ggsave(paste0(parentdir, "/images/ERV_vs_MAC10_corr_facet.pdf"),
  width=6, height=6)

##############################################################################
# ERVs vs masked ERVs
##############################################################################

maskgpdat$Category2 <- ifelse(substr(maskgpdat$Motif,4,5)=="CG",
  paste0("cpg_",maskgpdat$Type),
  maskgpdat$Type)

maskgpdat <- maskgpdat %>%
  mutate(Category2 = plyr::mapvalues(Category2, orderedcats1, orderedcats2))

maskgpdat$Category2 <- factor(maskgpdat$Category2, levels=orderedcats2)

rm2 <- maskgpdat %>%
  filter(nERVs_mask>20) %>%
  dplyr::select(-Category2)

p <- ggplot(data=maskgpdat[maskgpdat$nERVs_mask>20,],
    aes(x=ERV_rel_rate, y=ERV_rel_rate_mask*mean(ERV_rel_rate)/mean(ERV_rel_rate_mask)))+
  xlab(ervlab)+
  ylab(ervmasklab)+
  format_mac_comp+
  theme_mac_comp()

p + format_mac_comp_nf
ggsave(paste0(parentdir, "/images/ERV_vs_mask_corr.pdf"),
	width=4, height=5)

corlabs <- maskgpdat %>%
  filter(nERVs_mask>20) %>%
	group_by(Category2) %>%
	summarise(cor=cor(ERV_rel_rate, ERV_rel_rate_mask, method="spearman")) %>%
  mutate(label=paste0("r=", round(cor,2))) %>%
  mutate(x=0.0001, y=0.1)

p + format_mac_comp_facet(corlabs)
ggsave(paste0(parentdir, "/images/ERV_vs_mask_corr_facet.pdf"),
  width=6, height=6)

maskgpdat %>%
  filter(nERVs_mask>20) %>%
  mutate(r1=ERV_rel_rate_mask*mean(ERV_rel_rate)/mean(ERV_rel_rate_mask)) %>%
  mutate(prop=r1/ERV_rel_rate) %>%
  ggplot(aes(x=nERVs_mask, y=log(prop)))+
    geom_point(alpha=0.4)+
    xlab("#ERVs in masked subtype")+
    ylab("log(masked rate/unmasked rate)")

ggsave(paste0(parentdir, "/images/nervs_vs_prop.pdf"),
  width=6, height=6)


##############################################################################
# MAC10 vs AV
##############################################################################
p <- ggplot(data=r5m,
    aes(x=MAC10_rel_rate, y=eur*mean(MAC10_rel_rate)/mean(eur)))+
  xlab(mac10lab)+
  ylab(avlab)+
  format_mac_comp+
  theme_mac_comp()

p + format_mac_comp_nf
ggsave(paste0(parentdir, "/images/MAC10_vs_AV_corr.pdf"),
	width=4, height=5)

corlabs <- r5m %>%
	group_by(Category2) %>%
	summarise(cor=cor(MAC10_rel_rate, eur, method="spearman")) %>%
  mutate(label=paste0("r=", round(cor,2))) %>%
  mutate(x=0.0001, y=0.1)

p + format_mac_comp_facet(corlabs)
ggsave(paste0(parentdir, "/images/MAC10_vs_AV_corr_facet.pdf"),
  width=6, height=6)

oc3 <- c("ALL", gsub(" ", "\n", orderedcats2))

r5m %>%
  mutate(Category2 = gsub(" ", "\n", Category2)) %>%
	mutate(Category2=factor(Category2, levels=oc3[-1])) %>%
  mutate(av_rel_rate=eur*mean(ERV_rel_rate)/mean(eur)) %>%
  # filter(Type %in% c("AT_CG", "AT_GC")) %>%
  mutate(num_GC_flank = str_count(paste(substr(Motif,1,3), substr(Motif,5,7)), "C|G")) %>%
  # mutate(highGC = ifelse(str_count(substr(Motif,1,7), "A|T")<=3,"High","Low")) %>%
  mutate(maxc = ifelse(prop_diff4>=2, ">=2", "<2")) %>%
ggplot(
		aes(x=prop_diff, group=factor(num_GC_flank),
			colour=factor(num_GC_flank)))+
  geom_vline(xintercept=1, linetype="dashed")+
  geom_density()+
  # geom_point(alpha=0.7, size=2)+
	scale_colour_brewer(palette="YlGn")+
	scale_x_continuous(trans = "log", breaks=c(0.25,0.5,1,2,4))+
	# scale_y_log10()+
  facet_wrap(~Category2, ncol=1, scales="free_y", strip.position="left")+
	# scale_shape(solid = FALSE)+
  xlab("1000G:ERV rate ratio")+
  # ylab(avlab)+
	guides(
		colour = guide_legend("# G:C base pairs \n in flanking sequence", override.aes = list(alpha=1)))+
  theme_bw()+
  theme(legend.position = "bottom",
    legend.title=element_text(size=12),
    axis.title.x=element_text(size=14),
    axis.title.y=element_text(size=14))
ggsave(paste0(parentdir, "/images/AT_GC_scatter.pdf"), width=8, height=8)

##############################################################################
# Plot heatmap of change in relative rates
##############################################################################
nbox <- length(unique(r5m$v2a))
nint <- nbox/4
xhi <- rep(1:4,4)*nint+0.5
xlo <- xhi-nint
yhi <- rep(1:4,each=4)*nint+0.5
ylo <- yhi-nint
f <- data.frame(xlo,xhi,ylo,yhi)

levs <- as.character(lapply(as.vector(levels(r5m$v2a)), reverse_chars))

for(j in 1:6){
  categ<-orderedcats[j]
  p1a<-rrheat3(r5m[r5m$Type==categ,])

  png(paste0(parentdir, "/images/rare_common_diff2", categ, "_panel.png"),
    height=5, width=5, units="in", res=300)
  pushViewport(viewport(width=unit(5, "in"), height=unit(5, "in")))
  grid.draw(ggplotGrob(p1a))
  dev.off()
}

# trim whitespace on panels with imagemagick mogrify
trimcmd <- paste0("mogrify -trim ",
  parentdir, "/images/rare_common_diff*_panel.png")
system(trimcmd)

plegend<-ggplot()+
  geom_tile(data=r5m[r5m$Category==categ,],
    aes(x=v3, y=v2a, fill=prop_diff4))+
  scale_fill_gradientn("Rp/Rs\n",
    colours=myPaletteBrBG(nbp),
    trans="log",
    breaks=c(0.5, 1, 2),
    labels=c("<0.5", "1", ">2"),
    limits=c(0.5, 2.2))+
  xlab("3' flank")+
  ylab("5' flank")

legend <- get_legend(plegend)
png(paste0(parentdir, "/images/rare_common_diff2_legend.png"),
  height=8, width=3, units="in", res=300)
grid.draw(legend)
dev.off()
