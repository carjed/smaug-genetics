# seqs <- unlist(c(covdir$Sequence, lapply(covdir$Sequence, revcomp)))
# write.table(seqs, paste0(parentdir, "/seqs.txt"),
#   col.names=F, row.names=F, quote=F, sep="\t")
#
# # Create fasta file from feature file
# # if(sum(dnmstmp$inside)>0){
# if(!file.exists(paste0(covbase, ".fa"))){
#
#   if(cov=="RR"){
#     modifybedcmd <- paste0("cut -f1-4 ", covbed, " | sed 's/chr//g' | ",
#       awkstr, " | sed /^X/d | sed /^Y/d > ", covbase, ".new.bed")
#   } else if (cov=="TIME"){
#     modifybedcmd <- paste0("cut -f1-3 ", covbed, " | sed 's/chr//g' | ",
#       awkstr, " | sed /^X/d | sed /^Y/d > ", covbase, ".new.bed")
#   } else if (cov=="GC"){
#     modifybedcmd <- paste0("cut -f1-5 ", covbed, " | sed 's/chr//g' | ",
#       awkstr, " | sed /^X/d | sed /^Y/d > ", covbase, ".new.bed")
#   } else {
#     modifybedcmd <- paste0("sed 's/chr//g' ", covbed,
#       " | sed /^X/d | sed /^Y/d > ", covbase, ".new.bed")
#   }
#
#   getfastacmd <- paste0("bedtools getfasta -fi ", parentdir,
#     "/reference_data/human_g1k_v37.fasta -bed ", covbase,
#     ".new.bed -fo ", covbase, ".fa")
#
#   system(modifybedcmd)
#   system(getfastacmd)
# }
#
# # Search feature fasta for motifs in list
# grepcmd <- paste0("grep -o -Ff ",
#   parentdir, "/seqs.txt ", covbase, ".fa | sort | uniq -c > ",
#   parentdir, "/testcounts.txt")
# system(grepcmd)
#
# motifcts <- read.table(paste0(parentdir, "/testcounts.txt"),
#   header=F, stringsAsFactors=F)
#
# names(motifcts) <- c("Count", "SEQ")
# motifcts$REVSEQ <- unlist(lapply(motifcts$SEQ, revcomp))
# motifcts$Sequence <- ifelse(substr(motifcts$SEQ,4,4) %in% c("A", "C"),
#   motifcts$SEQ, motifcts$REVSEQ)
#
# motifcts2 <- motifcts %>%
#   group_by(Sequence) %>%
#   summarise(Count=sum(Count))
# covdir2 <- merge(covdir, motifcts2, by=c("Sequence"))
#
# covdir3 <- merge(covdir2, motifdat, by=c("Category", "Sequence"))
# covdir3$exp <- covdir3$Count*covdir3$rel_prop*1.67e-06*1074
# covdir3$expa <- covdir3$Count*covdir3$rel_prop
# # covdir3$exp <- covdir3$Count*covdir3$rel_prop*(46802/35574420)
#
# exp <- sum(covdir3$exp)
# obs <- sum(dnmstmp$inside, na.rm=T)
#
# r2 <- ratelist[[4]] %>%
#   filter(Motif %in% dnmstmp$SEQ & Type %in% dnmstmp$Category) %>%
#   dplyr::select(Motif, nMotifs) %>%
#   distinct(nMotifs)
#
# t2 <- dnmstmp %>% group_by(inside) %>% tally
# obs_in <- t2[2,]$n
# obs_out <- t2[1,]$n
# total <- nrow(dnmstmp)
#
#
# motifs_in <- sum(covdir3$Count)
# motifs_out <- sum(r2$nMotifs)-motifs_in
#
# if(dir=="Up"){
#   test <- prop.test(c(obs, exp), c(total, total), alternative="greater")
#   # test <- prop.test(c(obs_in, obs_out), c(motifs_in, motifs_out), alternative="greater")
# } else {
#   test <- prop.test(c(obs, exp), c(total, total), alternative="less")
#   # test <- prop.test(c(obs_in, obs_out), c(motifs_in, motifs_out), alternative="less")
# }
# newrow <- data.frame(cov, dir=dir, n=total,
#   obs=obs, exp=exp,
#   propobs=test$estimate[1], propexp=test$estimate[2], pval=test$p.value)
# newrow
