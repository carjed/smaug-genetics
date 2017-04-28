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

# correct (HC)
sites_c_hc <- full_data$sites %>%
  filter(REF==AA) %>%
  mutate(SEQA=substr(Motif, cbp-i, cbp+i),
    SEQB=substr(Motif, cbp*3-i, cbp*3+i),
    Motif=paste0(SEQA, "(", SEQB, ")"))

prop <- nrow(sites_c_hc)/nrow(full_data$sites)

bindir <- paste0(parentdir, "/motif_counts/", nbp2, "-mers/full")
p1 <- "motifs_full.txt"
bins <- get_bins(bindir, p1)
mct <- get_mct(bins)


# grepl("^[[:upper:]]+$", s)
ancaggseq <- get_aggseq(sites_c_hc, mct)
# rm(masktmp)

ancgpdat <- ancaggseq %>%
  mutate(Type=gsub("cpg_", "", Category2)) %>%
  dplyr::select(Type, Motif,
		nERVs_anc=nERVs, nMotifs_anc=nMotifs,
		ERV_rel_rate_anc=rel_prop)

test_anc <- merge(ratelist[[4]], ancgpdat, by=c("Type", "Motif")) %>%
	filter(nERVs >= 10 & nERVs_anc >= 10) %>%
  # mutate(pval=prop.test())
  group_by(Type, Motif) %>%
  do(tidy(prop.test(c(.$nERVs*prop, .$nERVs_anc), c(.$nMotifs, .$nMotifs_anc)))) %>%
	dplyr::select(Type, Motif, estimate1, estimate2, statistic, p.value)

sig_anc <- test_anc %>%
  filter(p.value<0.05/nrow(qctestdat)) %>%
  mutate(prop=estimate1/estimate2,
    gt10=ifelse(abs(log(prop))>0.223, TRUE, FALSE)) %>%
	filter(gt10==TRUE) %>%
	arrange(desc(prop))
