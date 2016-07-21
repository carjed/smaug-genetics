require(ggplot2)
require(dplyr)
require(tidyr)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##############################################################################
# Read and process data
##############################################################################


# Read data
cat("Reading data...\n")
maxc <- read.table("/net/bipolar/jedidiah/mutation/maxc_7bp.txt", header=T, stringsAsFactors=F)
maxc$Category <- gsub("cpg_", "", maxc$Category2)
chrpf <- read.table("/net/bipolar/jedidiah/mutation/output/rocdat.7bp.txt", header=F)
names(chrpf) <- c("CHR", "POS", "BIN", "MU", "OBS", "Category", "SEQ", "MU_C", "MU_S")

# Remove CpGs and sites with mu=0
# chrpf<-chrpf[substr(chrpf$SEQ, 2, 3)!="CG" & chrpf$MU>0,]
chrpf <- chrpf[chrpf$MU>0,]

# Read DNMs
cat("Reading DNMs...\n")
dnms_full <-  read.table("/net/bipolar/jedidiah/mutation/reference_data/DNMs/GoNL_DNMs.txt", header=T, stringsAsFactors=F)
dnms_full <- dnms_full[,1:5]
names(dnms_full) <- c("ID", "CHR", "POS", "REF", "ALT")

gdnms <- read.table("/net/bipolar/jedidiah/mutation/reference_data/DNMs/goldmann_2016_dnms.txt", header=T, stringsAsFactors=F)
gdnms$ID <- "goldmann"
gdnms <- gdnms[,c(7,1,2,4,5)]
names(gdnms) <- c("ID", "CHR", "POS", "REF", "ALT")
gdnms$CHR <- gsub("chr", "", gdnms$CHR)

dnms_full <- rbind(dnms_full, gdnms)
dnms_full$CAT <- paste(dnms_full$REF, dnms_full$ALT, sep="")

# Manually remove bins near chr20 centromere
# chr22 <- chr22[ which(chr22$BIN<260 | chr22$BIN>300),]
dnms_full$Category[dnms_full$CAT=="AC" | dnms_full$CAT=="TG"] <- "AT_CG"
dnms_full$Category[dnms_full$CAT=="AG" | dnms_full$CAT=="TC"] <- "AT_GC"
dnms_full$Category[dnms_full$CAT=="AT" | dnms_full$CAT=="TA"] <- "AT_TA"
dnms_full$Category[dnms_full$CAT=="GA" | dnms_full$CAT=="CT"] <- "GC_AT"
dnms_full$Category[dnms_full$CAT=="GC" | dnms_full$CAT=="CG"] <- "GC_CG"
dnms_full$Category[dnms_full$CAT=="GT" | dnms_full$CAT=="CA"] <- "GC_TA"

dnms_full <- dnms_full %>%
  dplyr::select(ID, CHR, POS, Category)

orderedcats <- sort(unique(dnms_full$Category))

# for(i in 1:6){
#   catind <- orderedcats[i]
#   dnmsub <- dnms_full %>% filter(Category==catind)
#   outfile <- paste0(parentdir, "/reference_data/DNMs/GoNL_", catind, ".txt")
#   write.table(dnmsub, outfile, col.names=F, row.names=F, quote=F, sep="\t")
# }

# Duplicate data, merge with DNMs to get ID
# chrpf<-chrp
cat("Annotating with ID...\n")
chrpf <- merge(chrpf, dnms_full, by=c("CHR", "POS"), all.x=T)
chrpf$ID[is.na(chrpf$ID)] <- "all"

# Subset to non-DNMs and DNMs
cat("Splitting by DNM status...\n")
chrpfa <- chrpf[chrpf$ID=="all",]
chrpfa <- chrpfa[sample(nrow(chrpfa), 1000000),]
chrpfdnm <- chrpf[chrpf$ID!="all",]

# Combine data
cat("Creating combined data...\n")
chrp <- rbind(chrpfdnm, chrpfa) %>%
  group_by(Category.x) %>%
  mutate(prop=cumsum(OBS)/sum(OBS)) %>%
  arrange(MU, prop)
# chrp$prop <- cumsum(chrp$OBS)/sum(chrp$OBS)

chrp <- chrp %>% filter((SEQ %in% maxc$Sequence) & (Category.x %in% maxc$Category))

# chrpsub <- rbind(chrpfdnm[sample(nrow(chrpfdnm), ndnms),], chrpfa) %>%
#   arrange(MU)
# chrpsub$prop <- cumsum(chrpsub$OBS)/sum(chrpsub$OBS)
# nsamp <- 100000
chrp1 <- chrp %>%
  group_by(Category.x) %>%
  arrange(MU) %>%
  mutate(prop=cumsum(OBS)/sum(OBS)) %>%
  arrange(MU, prop) %>%
  mutate(ntile=ntile(MU, 1000), group="Logit")

chrp2 <- chrp %>%
  group_by(Category.x) %>%
  arrange(MU_S) %>%
  mutate(prop=cumsum(OBS)/sum(OBS)) %>%
  arrange(MU_S, prop) %>%
  mutate(ntile=ntile(MU_S, 1000), group="ERVs")

chrp3 <- chrp %>%
  group_by(Category.x) %>%
  arrange(MU_C) %>%
  mutate(prop=cumsum(OBS)/sum(OBS)) %>%
  arrange(MU_C, prop) %>%
  mutate(ntile=ntile(MU_C, 1000), group="Common")

auctmp1 <- chrp1 %>% summarise(AUC=1-sum(prop)/nrow(chrp))
auctmp2 <- chrp2 %>% summarise(AUC=1-sum(prop)/nrow(chrp))
auctmp3 <- chrp3 %>% summarise(AUC=1-sum(prop)/nrow(chrp))

full_auc_dat <- rbind_all(list(chrp1, chrp2, chrp3))

full_auc <- full_auc_dat %>%
  group_by(group, Category.x) %>%
  summarise(AUC=1-sum(prop)/n()) %>%
  spread(Category.x, AUC)

ggplot()+
  geom_line(data=full_auc_dat,
		aes(x=(1000-ntile)/1000, y=1-prop, group=group, colour=group), size=1.2)+
  geom_abline(intercept=0, slope=1)+
  facet_wrap(~Category.x)+
  # geom_ribbon(data=full_auc_bounds,
	# 	aes(x=(1000-ntile)/1000, y=1-val, ymin=1-ub, ymax=1-lb), alpha=0.2)+
	scale_colour_manual(values=cbbPalette[2:4])+
  coord_cartesian(xlim=c(0,1))+
	xlab("False Positive Rate")+
	ylab("True Positive Rate")+
  theme_bw()+
	theme(
		# legend.position="none",
		axis.title.x=element_text(size=20),
		axis.title.y=element_text(size=20))
  # scale_y_continuous(limits=c(0,1))+
  # scale_x_continuous(limits=c(0,1000), labels=c())
ggsave("/net/bipolar/jedidiah/mutation/images/pseudo_roc_curves.png", height=9, width=12)
