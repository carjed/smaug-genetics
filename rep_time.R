binw <- 100000
adj <- 2
nbp <- adj*2+1

summfile<-paste0("/net/bipolar/jedidiah/mutation/output/chr20.",nbp,"bp.expanded.summary")
binfile<-paste0("/net/bipolar/jedidiah/mutation/output/chr20.bin_out_",nbp,"bp.txt")
chr22 <- read.table(summfile, header=T, stringsAsFactors=F)
bins <- read.table(binfile, header=T, stringsAsFactors=F)




