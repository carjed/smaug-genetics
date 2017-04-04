nbp2 <- 7
data2 <- "mask"

datadir2 <- paste0(parentdir,
	"/output/", nbp2, "bp_", bink, "k_singletons_", data2)

# summfile <- paste0(datadir2, "/full.summary")
# binfile <- paste0(datadir2, "/full_motif_counts_", bin_scheme, ".txt")

getData <- function(datadir, bin_scheme){

  summfile <- paste0(datadir, "/full.summary")
  binfile <- paste0(datadir, "/full_motif_counts_", bin_scheme, ".txt")

  if(!file.exists(summfile)){
  	cat("Merged summary file does not exist---Merging now...\n")
  	combine_sites <- paste0(
  		"awk 'FNR==1 && NR!=1{while(/^CHR/) getline; } 1 {print} ' ",
  		datadir, "/chr*.expanded.summary > ", summfile)
  	system(combine_sites)
  }

  if(!file.exists(binfile)){
  	cat("Merged bin file does not exist---Merging now...\n")
  	combine_motifs <- paste0(
  		"awk 'FNR==1 && NR!=1{while(/^CHR/) getline; } 1 {print} ' ",
  		# datadir, "/chr*.motif_counts_", bin_scheme, ".txt > ", binfile)
      datadir, "/chr*.motif_counts.txt > ", binfile)
  	system(combine_motifs)
  }

  # Read in per-site summary data
  cat("Reading summary file:", summfile, "...\n")
  sites <- read.table(summfile, header=T, stringsAsFactors=F)
  sites$BIN <- ceiling(sites$POS/binw)

  # Read in motif counts per chromosome
  cat("Reading bin file:", binfile, "...\n")
  bins <- read.table(binfile, header=T, stringsAsFactors=F, check.names=F)

  # summarize motif counts genome-wide
  mct <- bins %>%
  	# dplyr::select(CHR, BIN, Sequence=MOTIF, nMotifs=COUNT) %>%
  	dplyr::select(CHR, Sequence=MOTIF, nMotifs=COUNT) %>%
  	group_by(Sequence) %>%
  	summarise(nMotifs=sum(nMotifs))

  aggseq <- sites %>%
  	group_by(Sequence, Category2, BIN) %>%
  	summarise(n=n()) %>%
  	summarise(nERVs=sum(n))
  aggseq <- merge(aggseq, mct, by="Sequence")
  aggseq$rel_prop <- aggseq$nERVs/aggseq$nMotifs

  out <- list()
  out$sites <- sites
  out$bins <- bins
  out$mct <- mct
  out$aggseq <- aggseq
  return(out)
}

maskdata <- getData(datadir2, bin_scheme)

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
