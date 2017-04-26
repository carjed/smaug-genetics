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

##############################################################################
# Get BRIDGES MAC10+ rates
##############################################################################
summfile <- paste0(parentdir, "/summaries/full.summary")
singfile <- paste0(parentdir, "/singletons/full.singletons")
bindir <- paste0(parentdir, "/motif_counts/", nbp, "-mers/full")

full_data <- getData(commonfile, singfile, bindir)

##############################################################################
# Plot heatmap of change in relative rates
##############################################################################
nbox<-length(unique(r5m$v2a))
nint<-nbox/4
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
