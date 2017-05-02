cbp <- 5
i <- 3
nbp2 <- 7

# # mispolarized (HC)
# sites_mp_hc <- full_data$sites %>%
#   filter(ALT==AA)
#
# # mispolarized (HC+LC)
# sites_mp_all <- full_data$sites %>%
#   filter(tolower(ALT)==tolower(AA))
#
# # correct (HC+LC)
# sites_c_all <- full_data$sites %>%
#   filter(tolower(REF)==tolower(AA))
#
# #
# sites_c_lc <- full_data$sites %>%
#   filter(tolower(REF)==AA | AA %in% c("-", "N", "."))

# maskaggseq <- get_aggseq(masktmp, maskmct)

# correct (HC)
sites_c_hc <- full_data$sites %>%
  # filter(REF==AA) %>%
  # filter(MASK==0 & REF==AA) %>%
  filter(tolower(REF)==tolower(AA))

prop <- nrow(sites_c_hc)/nrow(full_data$sites)
# prop <- 0.946

sites_c_hc <- sites_c_hc %>%
  filter(MASK==0) %>%
  mutate(SEQA=substr(Motif, cbp-i, cbp+i),
    SEQB=substr(Motif, cbp*3-i, cbp*3+i),
    Motif=paste0(SEQA, "(", SEQB, ")"))

# bindir <- paste0(parentdir, "/motif_counts/", nbp2, "-mers/full")
# p1 <- "motifs_full.txt"
# bins <- get_bins(bindir, p1)
# mct <- get_mct(bins)
# ancaggseq <- get_aggseq(sites_c_hc, mct)

maskdir <- paste0(parentdir, "/motif_counts/", nbp2, "-mers/mask")
p1 <- "motifs_mask.txt"
maskbins <- get_bins(maskdir, p1)
maskmct <- get_mct(maskbins)
# grepl("^[[:upper:]]+$", s)
ancaggseq <- get_aggseq(sites_c_hc, maskmct)
# rm(masktmp)

ancgpdat <- ancaggseq %>%
  mutate(Type=gsub("cpg_", "", Category2)) %>%
  dplyr::select(Type, Motif,
		nERVs_anc=nERVs, nMotifs_anc=nMotifs,
		ERV_rel_rate_anc=rel_prop)

ancgpdat <- merge(ratelist[[4]], ancgpdat, by=c("Type", "Motif"))

test_anc <- ancgpdat %>%
	# filter(nERVs >= 10 & nERVs_anc >= 10) %>%
  # mutate(pval=prop.test())
  group_by(Type, Motif) %>%
  do(tidy(prop.test(c(.$nERVs, .$nERVs_anc), c(.$nMotifs, .$nMotifs_anc*prop)))) %>%
  # do(tidy(prop.test(c(.$nERVs, .$nERVs_anc), c(.$nMotifs, .$nMotifs_anc*0.946)))) %>%
	dplyr::select(Type, Motif, estimate1, estimate2, statistic, p.value)

test_anc <- merge(test_anc, ancgpdat, by=c("Type", "Motif"))

sig_anc <- test_anc %>%
  filter(p.value<0.05/nrow(test_anc)) %>%
  mutate(prop=estimate1/estimate2,
    gt10=ifelse(abs(log(prop))>log(1.2), TRUE, FALSE)) %>%
	filter(gt10==TRUE) %>%
	arrange(desc(prop))

# anc_motifs <- test_anc %>%
#   ungroup() %>%
#   filter(p.value<0.05/nrow(test_anc)) %>%
#   dplyr::select(Motif)
#
# mask_motifs <- test_mask %>%
#   ungroup() %>%
#   filter(p.value<0.05/nrow(test_mask)) %>%
#   dplyr::select(Motif)
#
#
# sum(anc_motifs$Motif %in% mask_motifs$Motif)
