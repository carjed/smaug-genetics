nbp2 <- 7
data2 <- "mask"

# maskdir <- paste0(parentdir, "/motif_counts/", nbp2, "bp_", bink, "k_singletons_", data2)
# maskbinfile <- paste0(datadir, "/full_motif_counts_", bin_scheme, ".txt")
maskdir <- paste0(parentdir, "/motif_counts/", nbp2, "-mers/mask")
p1 <- "motifs_mask.txt"
maskbins <- get_bins(maskdir, p1)

# maskbins <- read.table(maskbinfile, header=T, stringsAsFactors=F, check.names=F)
maskmct <- get_mct(maskbins)

cbp <- 5
i <- 3
masktmp <- full_data$sites %>%
	filter(MASK==0) %>%
	mutate(SEQA=substr(Motif, cbp-i, cbp+i),
		SEQB=substr(Motif, cbp*3-i, cbp*3+i),
		Motif=paste0(SEQA, "(", SEQB, ")"))

maskaggseq <- get_aggseq(masktmp, maskmct)
rm(masktmp)

maskgpdat <- maskaggseq %>%
  mutate(Type=gsub("cpg_", "", Category2)) %>%
  dplyr::select(Type, Motif,
		nERVs_mask=nERVs, nMotifs_mask=nMotifs,
		ERV_rel_rate_mask=rel_prop)

test_mask <- merge(ratelist[[4]], maskgpdat, by=c("Type", "Motif")) %>%
	filter(nERVs >= 10 & nERVs_mask >= 10) %>%
  # mutate(pval=prop.test())
  group_by(Type, Motif) %>%
  do(tidy(prop.test(c(.$nERVs, .$nERVs_mask), c(.$nMotifs, .$nMotifs_mask)))) %>%
	dplyr::select(Type, Motif, estimate1, estimate2, statistic, p.value)

sig_mask <- test_mask %>%
  filter(p.value<0.05/nrow(qctestdat)) %>%
  mutate(prop=estimate1/estimate2,
    gt10=ifelse(abs(log(prop))>0.223, TRUE, FALSE)) %>%
	filter(gt10==TRUE) %>%
	arrange(desc(prop))

# sig_mask %>%
# 	group_by(gt10) %>%
#   summarise(n=n())
