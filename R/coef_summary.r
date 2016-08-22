coefs <- read.table("/net/bipolar/jedidiah/mutation/output/logmod_data/coefs/coefs_full.txt", header=F, stringsAsFactors=F)

names(coefs) <- c("Cov", "Est", "SE", "Z", "pval", "Sequence", "Category")

coefs$Category <- ifelse(substr(coefs$Sequence, 4, 5)=="CG", paste0("cpg_", coefs$Category), coefs$Category)

 coefs_dhs <- coefs %>%
  filter(Cov=="DHS", pval<0.05) %>%
  dplyr::select(Sequence, Category)

rates <- read.table("/net/bipolar/jedidiah/mutation/output/7bp_1000k_rates.txt", header=T, stringsAsFactors=F)
rates$Sequence <- substr(rates$Sequence, 0, 7)
rates$Category2 <- gsub("cpg_", "", rates$Category2)

rates2 <- rates %>%
  group_by(Category2) %>%
  mutate(prop=num/sum(num))

# Summarise by median
coef_summary <- coefs %>%
  filter(coefs$pval<0.05) %>%
  mutate(OR=exp(Est)) %>%
  group_by(Category2, Cov) %>%
  summarise(OR=median(OR)) %>%
  spread(key=Category2, value=OR)

# Round
csdf <- data.frame(coef_summary)
csdf[-1] <- round(csdf[,-1],4)
csdf

plotdat <- coefs %>%
  group_by(Cov, Category) %>%
  # group_by(Category, Sequence) %>%
  mutate(qval=p.adjust(pval, method="fdr")) %>%
  filter(qval<0.05) %>%
  filter(!grepl("Intercept|DP", Cov)) %>%
  mutate(n=n()) %>%
  filter(n>=10) %>%
  ungroup()

  plotdat <- coefs %>%
    group_by(Cov, Category) %>%
    # group_by(Category, Sequence) %>%
    mutate(qval=p.adjust(pval, method="fdr")) %>%
    filter(Cov=="DHS" & Category=="AT_TA" & Sequence=="TTAAAAA")

plotcts <- plotdat %>%
  group_by(Cov, Category) %>%
  # group_by(Category, Sequence) %>%
  summarise(n=n(), Estmax=max(Est)+0.4*abs(max(Est))) %>% #max(Est)+0.4*median(Est)) %>%
  mutate(nmotifs=ifelse(grepl("^cpg", Category),
    1024,
    ifelse(grepl("^GC", Category),
      3072, 4096)),
    label=n,
    label2=paste0(n, " (", round(n/nmotifs, 3)*100, "%)"))
plotcts$Category <- factor(plotcts$Category, levels=orderedcats)

orderedcats <- c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA", "cpg_GC_AT", "cpg_GC_CG", "cpg_GC_TA")
plotdat$Category <- factor(plotdat$Category, levels=orderedcats)

ggplot(plotdat, aes(fill=Category))+
  geom_hline(yintercept=0, linetype="dashed")+
  geom_label(data=plotcts,
    aes(x=Category, y=Estmax, label=label), alpha=0.1, size=5)+
  geom_violin(aes(x=Category, y=Est, fill=Category))+
  scale_fill_brewer(palette="Set1")+
      # box.padding=unit(1, 'lines'))+
  # geom_boxplot(width=0.2, outlier.shape = NA)+
  facet_wrap(~Cov, ncol=4, scales="free_y")+
  ylab("log odds of mutation")+
  theme_bw()+
  guides(fill = guide_legend(nrow = 3))+
  theme(legend.position=c(.75, .1),
  strip.text=element_text(size=16),
  axis.text.x=element_blank(),
  axis.title.x=element_blank(),
  axis.text.y=element_text(size=14),
  axis.title.y=element_text(size=16))
ggsave("/net/bipolar/jedidiah/mutation/images/coef_violin2.png", width=12, height=8)

coefs$Category <- factor(coefs$Category, levels=orderedcats)
coefs %>%
  filter(!grepl("Intercept|DP", Cov)) %>%
  filter(abs(Est)<5) %>%
ggplot(aes(fill=Category))+
  geom_hline(yintercept=0, linetype="dashed")+
  # geom_label(data=plotcts,
  #   aes(x=Category, y=Estmax, label=label), alpha=0.1, size=5)+
  geom_violin(aes(x=Category, y=Est, fill=Category))+
  scale_fill_brewer(palette="Set1")+
      # box.padding=unit(1, 'lines'))+
  # geom_boxplot(width=0.2, outlier.shape = NA)+
  facet_wrap(~Cov, ncol=4, scales="free_y")+
  ylab("log odds of mutation")+
  theme_bw()+
  guides(fill = guide_legend(nrow = 2))+
  theme(legend.position="bottom",
  strip.text=element_text(size=16),
  axis.text.x=element_blank(),
  axis.title.x=element_blank(),
  axis.text.y=element_text(size=14),
  axis.title.y=element_text(size=16))
ggsave("/net/bipolar/jedidiah/mutation/images/coef_violin_full.png", width=9, height=8)
