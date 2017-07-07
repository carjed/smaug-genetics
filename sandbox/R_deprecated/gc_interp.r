##############################################################################
# Interpolate GC effects for Jun
##############################################################################
plotdat_gc <- plotdat %>%
  filter(Cov=="GC" & Category=="A>G") %>%
  mutate(Odds=exp(Est), subtype=paste0(Category, "_", Sequence)) %>%
  dplyr::select(Category, Sequence, subtype, Odds) %>%
  filter(Odds<2)

plotdat_exp <- plotdat_gc %>%
  tidyr::expand(subtype, ind=seq(0:10)-1)

plotdat_gc <- merge(plotdat_exp, plotdat_gc, by="subtype") %>%
  mutate(Odds=Odds^ind, pctGC=ind*10, seq3=substr(Sequence, 3, 5))

# Get top 6 positive/negative effects
top <- plotdat_gc %>%
  ungroup() %>%
  filter(pctGC==100) %>%
  top_n(6, Odds) %>%
  dplyr::select(subtype)

bottom <- plotdat_gc %>%
  ungroup() %>%
  filter(pctGC==100) %>%
  top_n(-6, Odds) %>%
  dplyr::select(Sequence)

tb <- rbind(top, bottom)

plotdat_gc <- plotdat_gc %>%
  filter(subtype %in% tb$subtype)

ggplot(plotdat_gc, aes(x=pctGC, y=log(Odds), colour=subtype, group=subtype))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=0)+
  scale_x_continuous(expand=c(0,0.4),
    breaks=seq(0, 100, 10),
    labels=seq(0, 100, 10))+
  scale_colour_manual(values=iwhPalette[1:12])+
  ylab("log-odds of mutability")+
  xlab("%GC")+
  theme_bw()+
  theme(legend.position="bottom")

ggsave(paste0(parentdir, "/images/gc_effect.png"))

sitefile1 <- paste0(parentdir, "/output/logmod_data/motifs3/CGTAATT.txt")
sitefile2 <- paste0(parentdir, "/output/logmod_data/motifs3/TCGATTG.txt")

sites1 <- read.table(sitefile1, header=F, stringsAsFactors=F)
names(sites1) <- c("CHR", "POS", "Sequence", mut_cats, "DP")

sites1$GC <- gcCol(sites1,
  paste0(parentdir, "/reference_data/gc10kb.bed"))

sites2 <- read.table(sitefile2, header=F, stringsAsFactors=F)
names(sites2) <- c("CHR", "POS", "Sequence", mut_cats, "DP")

sites2$GC <- gcCol(sites2,
  paste0(parentdir, "/reference_data/gc10kb.bed"))

sites1$GP <- "CGT[A>G]ATT (+)"
sites2$GP <- "TCG[A>G]TTG (-)"

mround <- function(x,base){
  base*round(x/base)
}

# sites2 <- sites2 %>%
rbind(sites1, sites2) %>%
mutate(GC=mround(GC*100, 0.5)) %>%
group_by(GP, GC, AT_GC) %>%
tally() %>%
spread(AT_GC, n) %>%
setNames(c("GP", "GC", "nonmut", "nERVs")) %>%
na.omit() %>%
mutate(tot=nERVs+nonmut, pctmut=nERVs/tot) %>% #data.frame
filter(nERVs>2) %>%
ggplot(aes(x=GC, y=pctmut, colour=GP))+
  geom_point(aes(group=GP, size=nERVs))+
  geom_smooth(method="lm", aes(weight=tot))+
  facet_wrap(~GP)+
  # geom_line()+
  scale_colour_manual(values=iwhPalette[c(6,10)])+
  xlab("%GC")+
  ylab("ERV fraction")+
  theme_bw()
ggsave(paste0(parentdir, "/images/gc_effect_raw.png"))

# rbind(sites1, sites2) %>%
# # mutate(GC=mround(GC*100, 2)) %>%
# group_by(GP, GC, AT_GC) %>%
# filter(GC>0.3 & GC<0.6) %>%
# # tally() %>%
# # spread(AT_GC, n) %>%
# # setNames(c("GP", "GC", "nonmut", "mut")) %>%
# # na.omit() %>%
# # mutate(tot=mut+nonmut, pctmut=mut/tot) %>% #data.frame
# # filter(mut>3) %>%
# ggplot(aes(x=AT_GC, y=GC, fill=GP, group=AT_GC))+
#   # geom_point(alpha=0.3, shape=21)+
#   geom_violin()+
#   facet_wrap(~GP)+
#   # geom_line()+
#   scale_colour_manual(values=iwhPalette[c(6,10)])+
#   xlab("%GC")+
#   ylab("mutated")+
#   theme_bw()
# ggsave(paste0(parentdir, "/images/gc_effect_raw2.png"), width=8, height=6)
