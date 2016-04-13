##############################################################################
# Recombination Hotspot Analysis
# flag is for distance to nearest hotspot heatmap
# also including code for depth heatmap--needs work
##############################################################################

recomb_title <- paste0("Chr",chr," ",mac, " Relative Mutation Rates in Recombination Hotspots")
recomb_out <- paste0(imgdir,"/chr",chr,"_",mac,"_hotspot_rel_mut_rate.png")

sites <- read.table("new_sites.txt", header=T, stringsAsFactors=F)
aggsites <- aggregate(Rate~Start, data=sites, mean)

hot <- read.table("hotspot_counts.txt", header=T, stringsAsFactors=F)
hotm <- melt(hot, id=names(hot)[1:10])
AT <- hotm[grep("^AT", hotm$variable),]
GC <- hotm[grep("^GC", hotm$variable),]
AT$prop <- AT$value/AT$AT
GC$prop <- GC$value/GC$GC
hotm2 <- rbind(AT,GC)

hotm2 <- merge(hotm2, aggsites, by="Start")

#average relative rates across recombination hotspots
AT_s <- aggregate(value~variable, data=AT, sum)
AT_t <- aggregate(AT~variable, data=AT, sum)
AT_s$prop <- AT_s$value/AT_t$AT

GC_s <- aggregate(value~variable, data=GC, sum)
GC_t <- aggregate(GC~variable, data=GC, sum)
GC_s$prop <- GC_s$value/GC_t$GC

cwa2 <- rbind(AT_s, GC_s)
#cwa2 <- aggregate(prop~variable, data=hotm2, mean)

############# PLOT RELATIVE MUTATION RATES ACROSS HOTSPOTS
ggplot(hotm2, aes(x=Centre, y=prop, colour=Rate, width=End-Start))+
	geom_bar(position="identity", stat="identity")+
	facet_wrap(~variable)+
	scale_colour_gradientn("Recombination Rate (cM/Mb)", colours=myPalette(4))+
	geom_hline(data=cwa, aes(yintercept=avg))+
	geom_hline(data=cwa2, aes(yintercept=prop), linetype="dashed")+
	theme_bw()+
	theme(panel.border=element_blank())+
	ggtitle(recomb_title)+
	xlab("Position")+
	ylab("Relative Mutation Rate")
suppressMessages(ggsave(recomb_out, width=12, height=7))


############# PLOT RELATIVE MUTATION RATES FOR 96 SUBTYPES--RECOMBINATION HOTSPOTS ONLY
bins3 <- read.table("bin_out3.txt", header=T, stringsAsFactors=F)
chr22s <- chr22[chr22$DIST==1,]
chr22sm <- merge(chr22s, bins3, by="SEQ")
chr22sm2 <- merge(chr22sm, bins3, by.x="ALTSEQ", by.y="SEQ")
chr22sm2$COUNT=chr22sm2$COUNT.x+chr22sm2$COUNT.y

aggseq <- count(chr22sm2, c("Sequence", "Category", "COUNT"))
aggseq$rel_prop <- aggseq$freq/aggseq$COUNT

aggseq_a <- aggseq[grep("^A", aggseq$Category),]
aggseq_g <- aggseq[grep("^G", aggseq$Category),]

at_rel_rec_prop_title <- paste0("Chr",chr,": AT Relative Mutation Rate by Local Sequence\n(recombination hotspots only)")
gc_rel_rec_prop_title <- paste0("Chr",chr,": GC Relative Mutation Rate by Local Sequence\n(recombination hotspots only)")

at_rel_rec_out <- paste0(imgdir,"/chr",chr,"_",mac,"_AT_rel_rec.png")
gc_rel_rec_out <- paste0(imgdir,"/chr",chr,"_",mac,"_GC_rel_rec.png")

ggplot(aggseq_a, aes(x=Category, y=rel_prop, fill=Sequence))+
	geom_bar(position="dodge", stat="identity")+
	#geom_errorbar(limits_prop, position="dodge")+
	scale_fill_manual(values=myPalette(16))+
	ggtitle(at_rel_rec_prop_title)+
	ylab("Relative Mutation Rate")+
	theme_bw()+
	theme(panel.border=element_blank(),
	axis.ticks.x=element_blank(),
	panel.grid.major.y=element_blank())
ggsave(at_rel_rec_out)

ggplot(aggseq_g, aes(x=Category, y=rel_prop, fill=Sequence))+
	geom_bar(position="dodge", stat="identity")+
	#geom_errorbar(limits_prop, position="dodge")+
	scale_fill_manual(values=myPalette(16))+
	ggtitle(gc_rel_rec_prop_title)+
	ylab("Relative Mutation Rate")+
	theme_bw()+
	theme(panel.border=element_blank(),
	axis.ticks.x=element_blank(),
	panel.grid.major.y=element_blank())
ggsave(gc_rel_rec_out)


############# PLOT HEATMAP OF RECOMBINATION HOTSPOTS (INDICATOR METHOD)
# hotspot_heat_out <- paste0(imgdir,"/chr",chr,"_",mac,"_mutation_vs_hotspot_heatmap.png")

# chr22s <- chr22[chr22$DIST==1,]
# hotspot_agg <- count(chr22s, c("Category", "BIN"))

# ggplot(hotspot_agg, aes(x=BIN, y=Category, fill=freq))+
	# geom_raster()+
	# scale_fill_gradientn(colours=myPalette(4))+
	# theme_bw()+
	# theme(panel.border=element_blank(),
		# axis.text.x = element_text(angle = 90, hjust = 0.5),
		# axis.text.y = element_text(angle = 90, hjust = 0.5))
# suppressMessages(ggsave(hotspot_heat_out))

#RECOMBINATION HOTSPOTS--OLD METHOD
# hotspot_agg <- aggregate(DIST ~ BIN+Category, data=chr22, mean)

# ggplot(hotspot_agg, aes(x=BIN, y=Category, fill=DIST))+
	# geom_raster()+
	# scale_fill_gradientn(colours=myPalette(4))+
	# scale_x_continuous(breaks=seq(0,xmax,50))+
	# theme_bw()
# suppressMessages(ggsave(hotspot_heat_out))

#DEPTH	

# chr22$DP <- as.numeric(chr22$DP)
# chr22$DP[is.na(chr22$DP)]=mean(chr22$DP, na.rm=T)
# chr22$AVGDP <- chr22$DP/(chr22$AN/2)

# depth_heat_out <- paste0(imgdir,"/chr",chr,"_",mac,"_mutation_vs_depth_heatmap.png")

# aggdata <- aggregate(AVGDP ~ BIN+Category, data=chr22, mean)

#subset to exclude bins with very high avg. depth
# agg2 <- aggdata[aggdata$BIN>250,]

#Graphics--heatmap of avg. depth per bin for the 6 mutation categories
# ggplot(agg2, aes(x=BIN, y=Category, fill=AVGDP))+
	# geom_raster()+
	# scale_fill_gradientn(colours=myPalette(4))+
	# scale_x_continuous(breaks=seq(0,xmax,50))+
	# theme_bw()+
	# theme(panel.border=element_blank(),
		# axis.text.x = element_text(angle = 90, hjust = 0.5),
		# axis.text.y = element_text(angle = 90, hjust = 0.5))
# suppressMessages(ggsave(depth_heat_out))
