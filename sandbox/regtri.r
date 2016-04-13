# tris<-apply(expand.grid(bases, nts, bases), 1, paste, collapse="")

bases<-c("A", "C", "G", "T")
nts<-c("A", "C")
# Append columns to windowed count data for all motif lengths
for(j in ((nbp-1)/2-1):1){

  # Specify iteration motif length
  mlength <- (j+1)*2+1

  # Define rule for substring evaluation
  griddef <- paste(c(rep("bases", j), "nts", rep("bases", j)), collapse=",")

  # Evaluate substring rule and get vector of submotifs
  tris <- apply(eval(parse(text=paste("expand.grid(",griddef,")"))),
    1, paste, collapse="")

  # Loop through each substring and append column of
  # aggregated counts
  for(i in tris){
    # Generate regex string; j is fixed per iteration
    # (e.g., looking for internal 3-mers or 5-mers)
    # so we search for all 3-mers or 5-mers by allowing
    # any base preceding or following the internal motif
    # regtri <- paste0("^", "[A-Z]{", j, "}", i, "[A-Z]{", j, "}")
    regtri <- paste0("^[A-Z]", i, "[A-Z]")

    # Extract sequences matching this submotif
    z <- names(a3a1)[grepl(regtri, names(a3a1))]

    # Ensure motif match vector includes only sequences
    # corresponding to the appropriate motif length
    z <- z[nchar(head(gsub("_[A-Z]*", "", z)))==mlength]

    # Create column and append to df
    tripct <- a3a1 %>%
      mutate_(.dots=setNames(paste(z, collapse="+"), i)) %>%
      select_(.dots=i)
    a3a1 <- cbind(a3a1, tripct)
  }
}

# Select motif columns matching the output set
# from hierarchical selection
a3a2 <- a3a1 %>%
  select_(.dots=hierset)

# Run model with this set vs. all 5bp motifs

# a3a1a<-colSums(a3a1[,c(16,20:ncol(a3a1))])

# rf1<-rates_full[rates_full$Category2=="AT_CG",]
# rf1b<-rf1 %>%
#   group_by(Seq3) %>%
#   summarise(ll5=sum(logLik5), ll3=sum(logLik3)) %>%
#   mutate(L=2*(ll5-ll3), p=1-pchisq(L, 15)) %>%
#   arrange(-L)
