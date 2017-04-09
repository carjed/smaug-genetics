nbp2 <- 7
data2 <- "mask"

maskdir <- paste0(parentdir,
	"/output/", nbp2, "bp_", bink, "k_singletons_", data2)
maskbinfile <- paste0(datadir, "/full_motif_counts_", bin_scheme, ".txt")

maskbins <- read.table(maskbinfile, header=T, stringsAsFactors=F, check.names=F)
maskmct <- get_mct(maskbins)
maskaggseq <- get_aggseq(full_data$sites[full_data$sites$MASK==0,], maskmct)

cbp<-5
i<-3
maskgpdat <- maskaggseq %>%
  mutate(Type=gsub("cpg_", "", Category2),
    SEQA=substr(Sequence, cbp-i, cbp+i),
    SEQB=substr(Sequence, cbp*3-i, cbp*3+i),
    Motif=paste0(SEQA, "(", SEQB, ")")) %>%
  dplyr::select(Type, Motif, nERVs, nMotifs) %>%
  group_by(Type, Motif) %>%
  summarise(nERVs_mask=sum(nERVs), nMotifs_mask=sum(nMotifs)) %>%
  mutate(ERV_rel_rate_mask=nERVs_mask/nMotifs_mask)

qctestdat <- merge(ratelist[[4]], maskgpdat, by=c("Type", "Motif")) %>%
  # mutate(pval=prop.test())
  group_by(Type, Motif) %>%
  do(tidy(prop.test(c(.$nERVs, .$nERVs_mask), c(.$nMotifs, .$nMotifs_mask))))

qctestdat %>%
  filter(p.value<0.05/24576) %>%
  mutate(prop=estimate1/estimate2,
    gt10=ifelse(abs(log(prop))>0.223, TRUE, FALSE)) %>%
  group_by(gt10) %>%
  summarise(n=n())
