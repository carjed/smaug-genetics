##############################################################################
# Miscellaneous snippets to possibly include in prop.R
##############################################################################

##############################################################################
# process data for secondary file
##############################################################################
# spec$BIN <- ceiling(spec$POS/binw)
# spec$CAT <- paste(spec$REF, spec$ALT, sep="")

# spec$Category[spec$CAT=="AC" | spec$CAT=="TG"] <- "AT_CG"
# spec$Category[spec$CAT=="AG" | spec$CAT=="TC"] <- "AT_GC"
# spec$Category[spec$CAT=="AT" | spec$CAT=="TA"] <- "AT_TA"
# spec$Category[spec$CAT=="GA" | spec$CAT=="CT"] <- "GC_AT"
# spec$Category[spec$CAT=="GC" | spec$CAT=="CG"] <- "GC_CG"
# spec$Category[spec$CAT=="GT" | spec$CAT=="CA"] <- "GC_TA"

##############################################################################
# Kataegis plot (experimental--not informative yet)
##############################################################################
# chr22$PREV <- c(chr22$POS[1], chr22$POS[-length(chr22$POS)])
# chr22$DIST_NEXT <- chr22$POS-chr22$PREV
# ggplot(chr22, aes(x=POS, y=DIST_NEXT, colour=Category))+geom_point()+scale_y_log10()
# ggsave(imgdir,"/chr22_kataegis.png", width=10, height=4)

##############################################################################
# GC CONTENT HEATMAP	
##############################################################################
# aggdata <- aggregate(GC ~ BIN+Category, data=chr22, mean)

# ggplot(aggdata, aes(x=BIN, y=Category, fill=GC))+
	# geom_raster()+
	# scale_fill_gradientn(colours=myPalette(4))+
	# scale_x_continuous(breaks=seq(0,xmax,50))+
	# scale_y_discrete(breaks=NULL)+
	# ylab(NULL)+
	# ggtitle(gc_heat_title)+
	# theme_bw()+
	# theme(panel.border=element_blank(),
		# axis.text.x = element_text(angle = 90, hjust = 0.5),
		# axis.text.y = element_text(angle = 90, hjust = 0.5))
# suppressMessages(ggsave(gc_heat_out))
