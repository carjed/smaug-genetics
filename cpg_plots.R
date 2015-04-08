#names(chr22a) <- c("POS","PAIR","CPGI")
cpg_scale_title <- paste0("Chr",chr," ",mac, " CpG Mutations--scaled proportions")
cpg_scale_out <- paste0(imgdir,"/chr",chr,"_",mac,"_cpg_mutation_prop2.png")

cpg_dist_title <- paste0("Chr",chr," ",mac, " CpG Distribution by Mutation Type")
cpg_ident_out <- paste0(imgdir,"/chr",chr,"_",mac,"_cpg_dist_ident.png")
cpg_dodge_out <- paste0(imgdir,"/chr",chr,"_",mac,"_cpg_dist_dodge.png")

chr22b$CpG[chr22b$PAIR=="CG" | chr22b$PAIR=="GC"] <- 1
chr22b$CpG[chr22b$PAIR!="CG" & chr22b$PAIR!="GC"] <- 0

chr22cpg <- chr22b[(chr22b$CpG==1 & chr22b$CPGI==0),]
chr22cpg$CAT <- paste(chr22cpg$REF, chr22cpg$ALT, sep="")

chr22cpg$Category[chr22cpg$CAT=="GA" | chr22cpg$CAT=="CT"] <- "GC to AT"
chr22cpg$Category[chr22cpg$CAT=="GC" | chr22cpg$CAT=="CG"] <- "GC to CG"
chr22cpg$Category[chr22cpg$CAT=="GT" | chr22cpg$CAT=="CA"] <- "GC to TA"

############# PLOT DODGED AND STACKED DISTRIBUTIONS OF 3 CPG CATEGORIES
ggplot(chr22cpg, aes(x=POS, colour=Category, fill=Category, group=Category))+
	geom_histogram(binwidth=binw, position="identity", alpha=0.5)+
	ggtitle(cpg_dist_title)+
	scale_colour_manual(values=g)+
	theme_bw()+
	theme(panel.border=element_blank(),
		axis.text.x = element_text(angle = 90, hjust = 0.5),
		axis.text.y = element_text(angle = 90, hjust = 0.5))
ggsave(cpg_ident_out)

ggplot(chr22cpg, aes(x=POS, colour=Category, fill=Category, group=Category))+
	geom_histogram(binwidth=binw, position="dodge", alpha=0.5)+
	ggtitle(cpg_dist_title)+
	scale_colour_manual(values=g)+
	theme_bw()+
	theme(panel.border=element_blank(),
		axis.text.x = element_text(angle = 90, hjust = 0.5),
		axis.text.y = element_text(angle = 90, hjust = 0.5))
ggsave(cpg_dodge_out)

#CpG--Merge bins + summary files and process
countcpg <- count(chr22cpg, c("Category", "BIN"))
countcpg <- merge(countcpg, aggregate(freq~BIN, data=countcpg, sum), by="BIN")
countcpg$rel_prop <- countcpg$freq.x/countcpg$freq.y
countcpg <- countcpg[order(countcpg$BIN, countcpg$Category),]

############# PLOT RELATIVE MUTATION RATES FOR CPG CATEGORIES
ggplot(countcpg, aes(x=factor(BIN), y=rel_prop, colour=Category, fill=Category))+
	geom_bar(position="stack", stat="identity", alpha=0.5)+
	scale_x_discrete(breaks=seq(0,xmax,50))+
	xlab("Bin")+
	ylab("Proportion")+
	scale_colour_manual(values=g)+
	scale_fill_manual(values=g)+
	ggtitle(cpg_scale_title)+
	theme_bw()+
	theme(panel.border=element_blank(),
		axis.text.x = element_text(angle = 90, hjust = 0.5),
		axis.text.y = element_text(angle = 90, hjust = 0.5))
suppressMessages(ggsave(cpg_scale_out))