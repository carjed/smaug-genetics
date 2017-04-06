nbp2 <- 7
data2 <- "mask"

maskdir <- paste0(parentdir,
	"/output/", nbp2, "bp_", bink, "k_singletons_", data2)

maskdata <- getData(maskdir, bin_scheme)

cbp<-4
i<-3
maskgpdat <- maskdata$aggseq %>%
  mutate(Type=gsub("cpg_", "", Category2),
    SEQA=substr(Sequence, cbp-i, cbp+i),
    SEQB=substr(Sequence, cbp*3-i, cbp*3+i),
    Motif=paste0(SEQA, "(", SEQB, ")"))

maskgpdat <- maskgpdat %>%
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

rates_mask <- maskgpdat %>%
  # mutate(SEQ5=substr(Motif, 1, 5)) %>%
  dplyr::select(Category.x=Type, SEQ=Motif, MU_7M=ERV_rel_rate_mask)

# chrp$SEQ5 <- substr(chrp$SEQ, 2, 6)
chrp <- merge(chrp, rates_mask, by=c("Category.x", "SEQ"), all.x=T)

rates9 <- ratelist[[5]] %>%
  # mutate(SEQ5=substr(Motif, 1, 5)) %>%
  dplyr::select(Category.x=Type, SEQ=Motif, MU_9=ERV_rel_rate)

# chrp$SEQ5 <- substr(chrp$SEQ, 2, 6)
chrp <- merge(chrp, rates9, by=c("Category.x", "SEQ"), all.x=T)
