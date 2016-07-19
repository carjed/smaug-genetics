library(devtools)
# install_github('carjed/bedr')
library(bedr)
library(dplyr)

sites<-read.table(pipe("cut -f1-5,18,19 /net/bipolar/jedidiah/mutation/output/logmod_data/motifs/AT_CG/dp/AT_CG_AAAAACG_dp.txt"), header=F)

# Get recombination rate at each site
rcrCol <- function(sites, file){
  feat_ranges <- bed_to_granges(file, header=T)
  site_ranges <- GRanges(seqnames=paste0("chr",sites$V2),
                         ranges=IRanges(start=sites$V1, end=sites$V1))

  indices <- findOverlaps(site_ranges, feat_ranges, type="within", select="first")
  indices[is.na(indices)]<-0
  ind_df<-data.frame(V1=sites$V1, V2=sites$V2, indices)

  feat_df<-as.data.frame(feat_ranges)
  feat_df$indices<-seq_along(1:nrow(feat_df))
  rate_table <- merge(ind_df, feat_df, by="indices", all.x=T, incomparables=0) %>%
    arrange(V2, V1)

  rates<-rate_table$id
  rates[is.na(rates)]<-0
  return(as.numeric(rates))
}

binaryCol <- function(sites, bedfile){
  feat_ranges <- bed_to_granges(bedfile, header=F)
  site_ranges <- GRanges(seqnames=paste0("chr",sites$V2),
                         ranges=IRanges(start=sites$V1, end=sites$V1))
  return(as.integer(site_ranges %within% feat_ranges))
}

repCol <- function(sites, repfile, binwidth){
  reptime <- read.table(repfile, header=F, stringsAsFactors=F, sep="\t")
  names(reptime) <- c("CHR", "END", "TIME")
  reptime2<-reptime %>%
    group_by(CHR) %>%
    arrange(CHR, END) %>%
    mutate(START=imputeStart(END))

  feat_ranges <- GRanges(seqnames=paste0("chr",reptime2$CHR),
                         ranges=IRanges(start=reptime2$START, end=reptime2$END), id=reptime2$TIME)
  site_ranges <- GRanges(seqnames=paste0("chr",sites$CHR),
                         ranges=IRanges(start=sites$POS, end=sites$POS))

  indices <- findOverlaps(site_ranges, feat_ranges, type="within", select="first")
  indices[is.na(indices)]<-0
  ind_df <- data.frame(POS=sites$POS, CHR=sites$CHR, indices)

  feat_df <- as.data.frame(feat_ranges)
  feat_df$indices <- seq_along(1:nrow(feat_df))
  rate_table <- merge(ind_df, feat_df, by="indices", all.x=T, incomparables=0) %>%
    arrange(CHR, POS)

  rates <- rate_table$id
  rates[is.na(rates)] <- 0
  return(as.numeric(rates))
  # return(as.integer(site_ranges %within% feat_ranges))
}

gcfile<-"/net/bipolar/jedidiah/mutation/output/3bp_10k/full_bin.txt"
gcCol<-function(sites, gcfile){

  incmd <- paste0("cut -f1,4,5 ", gcfile)
  gcbins <- read.table(pipe(incmd), header=T, stringsAsFactors=F)
  gcbins$start <- gcbins$BIN*10000-10000+1
  gcbins$end <- gcbins$start+10000-1

  gcbins<-gcbins %>% arrange(CHR, start)

  feat_ranges <- GRanges(seqnames=gcbins$CHR,
                         ranges=IRanges(start=gcbins$start, end=gcbins$end), id=gcbins$prop_GC)
  site_ranges <- GRanges(seqnames=paste0("chr",sites$CHR),
                         ranges=IRanges(start=sites$POS, end=sites$POS))

  indices <- findOverlaps(site_ranges, feat_ranges, type="within", select="first")
  indices[is.na(indices)]<-0
  ind_df <- data.frame(POS=sites$POS, CHR=sites$CHR, indices)

  feat_df <- as.data.frame(feat_ranges)
  feat_df$indices <- seq_along(1:nrow(feat_df))
  rate_table <- merge(ind_df, feat_df, by="indices", all.x=T, incomparables=0) %>%
    arrange(CHR, POS)

  rates <- rate_table$id
  rates[is.na(rates)] <- 0
  return(as.numeric(rates))
  # return(as.integer(site_ranges %within% feat_ranges))
}

gctest<-gcCol(sites, "/net/bipolar/jedidiah/mutation/output/3bp_10k/full_bin.txt")

# Loop to add histone marks to site data
hists<-c("H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3", "H3K27ac", "H3K27me3", "H3K36me3")

dflist <- list()
for(i in 1:length(hists)){
  mark <- hists[i]
  file <- paste0("/net/bipolar/jedidiah/mutation/reference_data/histone_marks/broad/sort.E062-", mark, ".bed")
  hist <- binaryCol(sites, file)
  dflist[[i]] <- hist
}

df <- as.data.frame(do.call(cbind, dflist))
names(df) <- hists
sites<-cbind(sites, df)

parentdir<-"/net/bipolar/jedidiah/mutation"
repfile <- paste0(parentdir, "/reference_data/lymph_rep_time.txt")
inside<-repCol(sites, repfile)

reptime$START <- c(0, reptime$POS[-(length(reptime$POS))])

# Function determines start of interval from value in previous row
imputeStart<-function(ends){
  starts<-c(0, ends[-(length(ends))])
  return(starts)
}



approxData <- data.frame(
  with(reptime,
       approx(POS, TIME, xout = seq(1, n, by = 10), method = "linear")
  ),
  method = "approx()"
)


sites$CpGI <- binaryCol(sites, "/net/bipolar/jedidiah/mutation/reference_data/cpg_islands_sorted.bed")
sites$RR <- rcrCol(sites, "/net/bipolar/jedidiah/mutation/reference_data/recomb_rate.bed")
sites$LAMIN <- binaryCol(sites, "/net/bipolar/jedidiah/mutation/reference_data/lamin_B1_LADS2.bed")
sites$DHS <- binaryCol(sites, "/net/bipolar/jedidiah/mutation/reference_data/DHS.bed")
