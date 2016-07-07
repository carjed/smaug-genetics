
cat("Calculating likelihoods under different motif lengths...\n")
# Calculate likelihoods
rates_full$logLik1 <- dbinom(rates_full$num,
  rates_full$COUNT,
  rates_full$rel_prop1, log=T)

rates_full$logLik3 <- dbinom(rates_full$num,
  rates_full$COUNT,
  rates_full$rel_prop3, log=T)

rates_full$logLik5 <- dbinom(rates_full$num,
  rates_full$COUNT,
  rates_full$rel_prop5, log=T)

rates_full$logLik7 <- dbinom(rates_full$num,
  rates_full$COUNT,
  rates_full$rel_prop, log=T)

ra1<-rates_full %>%
  group_by(Category2) %>%
  summarise(L1=-2*sum(logLik1),
    L3=-2*sum(logLik3),
    L5=-2*sum(logLik5),
		L7=-2*sum(logLik7))

ra1a<-gather(ra1, model, log, L1:L7)

ra1a$k <- c(rep(1,9),
  c(rep(16,3), rep(12,3), rep(4,3)),
  c(rep(256,3), rep(192,3), rep(64,3)),
	c(rep(4096,3), rep(3072,3), rep(1024,3)))

ra1a$AIC <- 2*ra1a$k+ra1a$log

ra1c <- merge(ra1a, rates1, by="Category2")
# ra1c$BIC <- ra1c$k*log(ra1c$num1)+ra1c$log

ra1b <- gather(ra1c, stat, L, c(log, AIC))
levels(ra1b$model) <- c("1", "3", "5", "7")

levels(ra1b$stat) <- c("-2ln(L)", "AIC")
names(ra1b) <- c("Category", "Motif_Length", "k", "COUNT",
  "num1", "rel_prop", "Stat", "L")

# Plot AIC and -2log(L) for each category
cat("Plotting likelihood curves...\n")
ggplot(ra1b, aes(x=Motif_Length, y=L, group=Stat, colour=Stat))+
  scale_colour_brewer(palette="Set1")+
  geom_point()+
  geom_line()+
  facet_wrap(~Category, scales="free")+
  theme_bw()+
  xlab("Motif Length")+
  theme(axis.title.y=element_blank(),
    legend.title=element_blank())

ggsave("/net/bipolar/jedidiah/mutation/images/compare_AIC.png")

rf2<-rates_full %>%
  group_by(Category2, Seq3) %>%
  summarise(s3=-2*sum(logLik3), s5=-2*sum(logLik5)) %>%
  ungroup() %>%
  arrange(Category2, desc(s5)) %>%
  group_by(Category2) %>%
  mutate(D3_5=s5-s3,s3s=sum(s3)) %>%
  arrange(D3_5) %>%
  mutate(D3_5c=cumsum(D3_5),
    L=s3s+D3_5c,
    rk=rank(-L),
    rk2=max(rk)-rk+16*rk,
    rk3=rk2,
    gp="3>5")

rf3 <- rates_full %>%
  group_by(Category2, Seq5) %>%
  summarise(s5=-2*sum(logLik5), s7=-2*sum(logLik7)) %>%
  ungroup() %>%
  arrange(Category2, desc(s7)) %>%
  group_by(Category2) %>%
  mutate(D5_7=s7-s5,s5s=sum(s5)) %>%
  arrange(D5_7) %>%
  mutate(D5_7c=cumsum(D5_7),
    L=s5s+D5_7c,
    rk=rank(-L),
    rk2=max(rk)-rk+16*rk,
    rk3=min(rk2)+min(rk2)*rk2/max(rk2),
    gp="5>7")

names(rf3)[2] <- "Sequence"
names(rf2)[2] <- "Sequence"

rfc <- rbind(dplyr::select(rf2, Category2, Sequence, L, rk, rk2, rk3, gp),
  dplyr::select(rf3, Category2, Sequence, L, rk, rk2, rk3, gp))

write.table(ra1b, "/net/bipolar/jedidiah/mutation/ra1b.txt", col.names=T, row.names=F, quote=F, sep="\t")

l1ll <- ra1b %>%
  filter(Motif_Length %in% c("L1", "L3"), Stat=="log") %>%
  dplyr::mutate(Sequence=substr(Category, 1,1), rk=5, rk2=k,
    rk3=min(rk2)+min(rk2)*rk2/max(rk2), gp=Motif_Length) %>%
  dplyr::select(Category2=Category, Sequence, L, rk, rk2, rk3, gp)

rfc2 <- merge(rbind(l1ll, rfc), rates1, by="Category2")

rfc2$AIC <- 2*rfc2$rk2 + rfc2$L
rfc2$BIC <- rfc2$rk2*log(rfc2$num1) + rfc2$L

rfc3<-gather(rfc2, z=AIC:BIC)

ggplot(rfc3, aes(x=rk3, y=value, group=key, colour=key))+
  scale_x_continuous(expand = c(.15, .15))+
  scale_y_continuous(expand = c(.1, .1))+
  scale_colour_brewer(palette="Dark2")+
  geom_point(size=4, aes(alpha=0.4))+
  facet_wrap(~Category2, scales="free")+
  theme_bw()+
  ylab("-2log(L)")+
  theme(axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.y=element_text(size=16),
    axis.text.y=element_text(size=14),
    strip.text.x = element_text(size=16),
    legend.title=element_blank(),
    legend.text=element_text(size=16))

ggsave("/net/bipolar/jedidiah/mutation/images/compare_AIC-BIC_full.png",
  width=12, height=12)

ggplot(rfc2, aes(x=rk3, y=L))+
  scale_x_continuous(expand = c(.15, .15))+
  scale_y_continuous(expand = c(.1, .1))+
  scale_colour_brewer(palette="Set1")+
  geom_point(size=4, aes(colour=gp))+
  geom_line()+
  geom_text(size=6, angle=10,
    aes(label=ifelse(rk<5, as.character(Sequence), ''),
      hjust=rep(c(1,-0.5), length.out=length(Sequence)),
      vjust=0.5))+#rep(c(.5,-.5), length.out=length(Sequence))))+
  facet_wrap(~Category2, scales="free")+
  theme_bw()+
  ylab("-2log(L)")+
  theme(axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.y=element_text(size=16),
    axis.text.y=element_text(size=14),
    strip.text.x = element_text(size=16),
    legend.title=element_blank(),
    legend.text=element_text(size=16))

ggsave("/net/bipolar/jedidiah/mutation/images/compare_LL_full.png",
  width=12, height=12)
