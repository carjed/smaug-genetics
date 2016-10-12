suppressMessages(usePackage(ggplot2))
suppressMessages(usePackage(dplyr))
suppressMessages(usePackage(tidyr))
suppressMessages(usePackage(reshape2))
suppressMessages(usePackage(RColorBrewer))
suppressMessages(usePackage(MASS))
suppressMessages(usePackage(speedglm))
suppressMessages(usePackage(boot))
suppressMessages(usePackage(devtools))
suppressMessages(usePackage(psych))

source("./R/get_functions.r")
source("./R/roc_functions.r")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#e7baea")
myPaletteCat <- colorRampPalette(brewer.pal(12, "Paired"))

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
cols <- myPaletteCat(12)[c(10,2,1,8,4,3,12,6,5)] #<- colors if using this ordering

ptm <- proc.time()
cat("Plotting K-mer heatmaps...\n")
source("./R/kmer_heatmap.r")
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
