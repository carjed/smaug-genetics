source("./R/get_functions.r")
source("./R/roc_functions.r")

suppressMessages(usePackage(ggplot2))
suppressMessages(usePackage(dplyr))
suppressMessages(usePackage(tidyr))
suppressMessages(usePackage(RColorBrewer))
suppressMessages(usePackage(MASS))
suppressMessages(usePackage(speedglm))
suppressMessages(usePackage(boot))
suppressMessages(usePackage(devtools))
suppressMessages(usePackage(psych))
suppressMessages(usePackage(lmtest))
suppressMessages(usePackage(fmsb))
suppressMessages(usePackage(stringr))
suppressMessages(usePackage(hexbin))
suppressMessages(usePackage(cowplot))
suppressMessages(usePackage(grid))
suppressMessages(usePackage(gridExtra))
suppressMessages(usePackage(grid))
suppressMessages(usePackage(gtable))

# Set custom library path for bedr package
libpath <- "~/R/x86_64-pc-linux-gnu-library/2.13"
suppressMessages(require(bedr, lib.loc=libpath, quietly=T))

parentdir <- dirname(getwd())

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#e7baea")

iwhPalette <- c("#cd5431", "#a14ad9", "#67b03f", "#604dad", "#c79931",
  "#cc4498", "#4d9f83", "#b54f50", "#5d8cb6", "#7f7d48", "#a67abe", "#965571")
myPaletteCat <- colorRampPalette(brewer.pal(12, "Paired"))
myPaletteCatN <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
myPaletteO <- colorRampPalette(rev(brewer.pal(11, "RdBu")), space="Lab")
myPaletteBrBG <- colorRampPalette(rev(brewer.pal(11, "BrBG")), space="Lab")

rb <- c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
g <- myPaletteG(6)[1:3]
rbg<-c(myPaletteB(6)[1:3],myPaletteR(6)[1:3], myPaletteG(6)[1:3])

orderedcats <- c("AT_CG", "AT_GC", "AT_TA",
  "GC_AT", "GC_CG", "GC_TA",
  "cpg_GC_AT", "cpg_GC_CG", "cpg_GC_TA")
# cols <- myPaletteCat(12)[c(8,10,12,2,4,6,1,3,5)] #<- colors if using this ordering

# reordered to group ts and tv
orderedcats1 <- c("AT_GC", "GC_AT", "cpg_GC_AT",
  "AT_CG", "GC_CG", "cpg_GC_CG",
  "AT_TA", "GC_TA", "cpg_GC_TA")
orderedcats2 <- c("A>G", "C>T", "CpG>TpG",
  "A>C", "C>G", "CpG>GpG",
  "A>T", "C>A", "CpG>ApG")
cols <- myPaletteCat(12)[
  c(10,2,1,
    8,4,3,
    12,6,5)] #<- colors if using this ordering

orderedcats1 <- c("AT_GC", "AT_CG", "AT_TA",
  "GC_AT", "GC_TA", "GC_CG",
  "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG")
orderedcats2 <- c("A>G", "A>C", "A>T",
  "C>T (non-CpG)", "C>A (non-CpG)", "C>G (non-CpG)",
  "CpG>TpG", "CpG>ApG", "CpG>GpG")
cols <- myPaletteCat(12)[
  c(10,8,12,
    2,4,6,
    1,3,5)] #<- colors if using this ordering

ptm <- proc.time()
cat("Plotting K-mer heatmaps...\n")
source("./R/kmer_heatmaps.r")
tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

ptm <- proc.time()
cat("Comparing ERVs and polymorphisms...\n")
source("./R/sing_vs_com.r")
tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

ptm <- proc.time()
cat("Validating models on de novo mutations...\n")
source("./R/roc_max.r")
tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")

ptm <- proc.time()
cat("Analyzing genomic features...\n")
source("./R/coef_summary.r")
tottime <- (proc.time()-ptm)[3]
cat("Done (", tottime, "s)\n")
