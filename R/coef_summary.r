coefs <- read.table("/net/bipolar/jedidiah/mutation/output/logmod_data/coefs/coefs_full.txt", header=F, stringsAsFactors=F)

names(coefs) <- c("Cov", "Est", "SE", "Z", "pval", "Sequence", "Category2")

rates <- read.table("/net/bipolar/jedidiah/mutation/output/7bp_1000k_rates.txt", header=T, stringsAsFactors=F)
rates$Sequence <- substr(rates$Sequence, 0, 7)
rates$Category2 <- gsub("cpg_", "", rates$Category2)

rates2 <- rates %>%
  group_by(Category2) %>%
  mutate(prop=num/sum(num))

# Weighted mean version
# mdat <- merge(coefs, rates2, by=c("Category2", "Sequence")) %>%
#   mutate(OR=exp(Est), w.OR=OR*prop) %>%
#   group_by(Category2, Cov) %>%
#   summarise(OR=sum(w.OR)) %>%
#   spread(key=Category2, value=OR)

# Summarise by median
coef_summary <- coefs %>%
  # filter(coefs$pval<0.01) %>%
  mutate(OR=exp(Est)) %>%
  group_by(Category2, Cov) %>%
  summarise(OR=median(OR)) %>%
  spread(key=Category2, value=OR)

# Round
csdf <- data.frame(coef_summary)
csdf[-1] <- round(csdf[,-1],4)
csdf
