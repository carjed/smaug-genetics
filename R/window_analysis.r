require(GRanges)
require(biovizBase)
require(ggbio)
data(hg19IdeogramCyto, package = "biovizBase")
# chr20<-sites %>% filter(CHR==20)

getOption("biovizBase")$cytobandColor

# chr20a <- sites %>%
#   mutate(Type=gsub("cpg_", "", Category2),
#     SEQA=substr(Sequence, cbp-i, cbp+i),
#     SEQB=substr(Sequence, cbp*3-i, cbp*3+i),
#     Motif=paste0(SEQA, "(", SEQB, ")"))

site_ranges <- GRanges(
  seqnames=paste0("chr",full_data$sites$CHR),
  ranges=IRanges(start=full_data$sites$POS, end=full_data$sites$POS)
)

s2 <- nearest(site_ranges, hg19IdeogramCyto)
sites <- cbind(full_data$sites, as.data.frame(hg19IdeogramCyto[s2])[-c(1:5)])

sites$name <- as.character(sites$name)
sites$gieStain <- as.character(sites$gieStain)
# sites$band <- paste0(sites$CHR,sites$name)

sites$bin <- ceiling(sites$POS/binw)

i<-1
cbp<-5
band_counts <- sites %>%
  mutate(Sequence=substr(Motif,1,9)) %>%
  mutate(Type=gsub("cpg_", "", Category2),
    SEQA=substr(Sequence, cbp-i, cbp+i),
    SEQB=substr(Sequence, cbp*3-i, cbp*3+i),
    Motif=paste0(SEQA, "(", SEQB, ")"),
    band=paste0(CHR, name)) %>%
  filter(band!="18p11.1") %>%
  group_by(CHR, band, gieStain, Type, Motif) %>%
  summarise(n=n())

window_counts <- sites %>%
  mutate(Type=gsub("cpg_", "", Category2),
    SEQA=substr(Sequence, cbp-i, cbp+i),
    SEQB=substr(Sequence, cbp*3-i, cbp*3+i),
    Motif=paste0(SEQA, "(", SEQB, ")"),
    band=paste0(CHR, name)) %>%
  filter(band!="18p11.1") %>%
  group_by(CHR, bin, Type, Motif) %>%
  summarise(n=n())

# chr20bins<-read.table("/net/bipolar/jedidiah/mutation/output/3bp_1000k_singletons_full/chr20.motif_counts_binned.txt", header=T, stringsAsFactors=F)
band_motif_counts <- read.table(paste0(parentdir,
    "/output/3bp_1000k_singletons_full/full_motif_counts_band.txt"),
  header=F, stringsAsFactors=F, skip=1)

names(band_motif_counts)<-c("CHR", "START", "END", "name", "gieStain", "Motif", "Count")

band_motif_counts <- band_motif_counts %>%
  mutate(CHR = gsub("chr", "", CHR),
    Length = END-START,
    band=paste0(CHR, name))

band_dat <- merge(band_counts, band_motif_counts, by=c("CHR", "band", "gieStain", "Motif"))
band_dat$ERV_rel_rate <- band_dat$n/band_dat$Count

pca_input <- band_dat %>%
# filter(Type=="AT_CG") %>%
  mutate(normrate=scale(ERV_rel_rate),
    subtype=paste0(Type, "_", Motif)) %>%
  dplyr::select(CHR, band, gieStain, START, Length, subtype, ERV_rel_rate) %>%
  spread(subtype, ERV_rel_rate)

total_band_counts <- band_counts %>%
  group_by(CHR, band) %>%
  summarise(nERVs=sum(n))

pca_input[is.na(pca_input)] <- 0
pcdat <- prcomp(scale(pca_input[,-c(1,2,3,4,5)]))
pcd2 <- cbind(pca_input[,1:5], as.data.frame(pcdat$x)) # %>% mutate(band=paste0(CHR,name))


pcd2 <- merge(pcd2, total_band_counts, by=c("CHR", "band"))

pcd2$gieStain <- as.factor(pcd2$gieStain)
levels(pcd2$gieStain) <- c("acen", "gneg", "gpos25", "gpos50", "gpos75", "gpos100", "gvar")


# tsnedat<-tsne(chr20_atcg_pca[,-c(1,2,3,4,5)], initial_dims=30, perplexity=30, max_iter=500)
# pcdat<-prcomp(scale(chr20_atcg_pca[,3:10]))
# tsnedat2<-tsne(pcdat$x, initial_dims=30, perplexity=30, max_iter=500)
tsnedat <- tsne(pcdat$x[,1:15], initial_dims=15, perplexity=10, max_iter=5000)
# cumsum((pcdat$sdev)^2) / sum(pcdat$sdev^2)
pcd2 <- cbind(pcd2, as.data.frame(tsnedat))

# pcd2.sav<-pcd2
pcd2.pc.clust <- pcd2 %>%
  # filter(gieStain=="acen" | grepl("9p2", band) | grepl("8p2", band) | grepl("^2q21", band) | grepl("^2q22", band))
  # filter(grepl("9p2", band) | grepl("8p2", band) | (PC2< -5 & PC1< -10))
  filter((PC2< -3.5 & PC1< 0))

clust1 <- pcd2 %>%
  filter((V1 > -25 & V1 < 25 & V2 > 65 & V2< 90)) %>%
  mutate(cluster="clust1")

clust2 <- pcd2 %>%
  filter((V1 > 15 & V1 < 30 & V2 > 50 & V2< 65)) %>%
  mutate(cluster="clust2")

clust3 <- pcd2 %>%
  filter((V1 > 40 & V1 < 75 & V2 > 10 & V2< 40)) %>%
  mutate(cluster="clust3")

clust4 <- pcd2 %>%
  filter((V1 > 62 & V1 < 100 & V2 > -50 & V2< -20)) %>%
  mutate(cluster="clust4")

clusterdat <- bind_rows(list(clust1, clust2, clust3, clust4))

pcd.path <- pcd2 %>%
  filter(CHR==8)# | (V2>30 & V1<0))

fullcols <- getOption("biovizBase")$cytobandColor

giemsacols <- fullcols[names(fullcols) %in% unique(pcd2$gieStain)]
giemsacols["gvar"] <- "#C6A6E6"
# giemsacols <- unlist(giemsacols)

ggplot()+
  geom_point(data=pcd2, aes(x=PC1, y=PC2, fill=gieStain, size=log10(nERVs)),
         colour="black",pch=21, alpha=0.8)+
  # geom_label_repel(data=pcd2.pc.clust, max.iter=5000, xlim=c(0, NA), ylim=c(NA, 5), fill="white", alpha=1, inherit.aes=F,
  #   aes(x=PC1, y=PC2, label = band, colour=factor(CHR)),
  # box.padding = unit(0.6, "lines")) +
  # scale_color_brewer(type='div', palette=2)+
  scale_fill_manual(values=giemsacols)+
  scale_colour_manual(values=myPaletteCatN(8))+
  # geom_path(data=pcd.path, aes(x=PC1, y=PC2), linetype=2)+
  theme_bw()+
  theme(#legend.title=element_blank(),
    # legend.key.size = unit(1.1, "cm"),
    axis.text.x=element_text(size=14),
    axis.text.y=element_text(size=14),
    axis.title.y = element_text(size=16),
    axis.title.x=element_text(size=16))+
  guides(colour=FALSE,
    fill = guide_legend(override.aes = list(size=6)))
ggsave(paste0(parentdir, "/images/chr20_bin_pca.png"), height=12, width=12)

ggplot(data, aes_string(x="V1", y="V2", color=var_cluster)) +
geom_point(size=0.25) +
guides(colour=guide_legend(override.aes=list(size=6))) +
xlab("") + ylab("") +
ggtitle("") +
theme_light(base_size=20) +
theme(axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      legend.direction = "horizontal",
      legend.position = "bottom",
      legend.box = "horizontal") +
  scale_colour_brewer(palette = palette)

  tsne_dist <- dist(pcd2[,c("V1", "V2")])
  cl <- hclust(tsne_dist)
cl2 <- kmeans(scale(pcd2[,c("V1", "V2")]), 4)

pcd2$tsne_clust <- factor(cl2$cluster)
  ## Creating hierarchical cluster model, and assigning the result to the data used to create the tsne
  fit_cluster_hierarchical=hclust(dist(scale(d_tsne_1)))

  ## setting 3 clusters as output
  pcd2$tsne_clust = factor(cutree(cl, k=4))

ggplot()+
  geom_point(data=pcd2, aes(x=V1, y=V2, fill=gieStain, size=log10(nERVs)),
         colour="black",
         pch=21, alpha=0.8)+
  # geom_label_repel(data=pcd2.tsne.clust1, max.iter=5000, xlim=c(NA, 0), ylim=c(NA, -55), fill="white", alpha=1, inherit.aes=F,
  #   aes(x=V1, y=V2, label = band, colour=factor(CHR)),
  # box.padding = unit(0.6, "lines")) +
  # geom_label_repel(data=pcd2.tsne.clust2, max.iter=5000, xlim=c(NA, -25), ylim=c(-60,-10), fill="white", alpha=1, inherit.aes=F,
  #   aes(x=V1, y=V2, label = band, colour=factor(CHR)),
  # box.padding = unit(0.1, "lines")) +
  # geom_label_repel(data=pcd2.cent.clust, max.iter=5000, xlim=c(NA, -50), ylim=c(30,NA), fill="white", alpha=1, inherit.aes=F,
  #   aes(x=V1, y=V2, label = band, colour=factor(CHR)),
  # box.padding = unit(0.6, "lines")) +
  # scale_color_brewer(type='div', palette=2)+
  scale_fill_manual(values=giemsacols)+
  # scale_colour_manual(values=c(iwhPalette, myPaletteCatN(8)))+
  # geom_path(data=pcd.path, aes(x=PC1, y=PC2), linetype=2)+
  theme_bw()+
  theme(#legend.title=element_blank(),
    # legend.key.size = unit(1.1, "cm"),
    axis.text.x=element_text(size=14),
    axis.text.y=element_text(size=14),
    axis.title.y = element_text(size=16),
    axis.title.x=element_text(size=16))+
  guides(colour=FALSE,
    fill = guide_legend(override.aes = list(size=6)))
ggsave(paste0(parentdir, "/images/chr20_bin_tsne.png"), height=12, width=12)

clusterdat$midpt <- clusterdat$START+clusterdat$Length/2

clust_markers <- GRanges(seqnames=paste0("chr", clusterdat$CHR),
                       ranges=IRanges(start=clusterdat$midpt,
                         end=clusterdat$midpt),
                       cluster=clusterdat$cluster,
                       yval=15)
data(hg19Ideogram, package = "biovizBase")
seqlengths(clust_markers) <- seqlengths(hg19Ideogram)[names(seqlengths(clust_markers))]
hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X", "Y")))

# clustchrs <- hg19
clustchrs<-hg19[seqnames(hg19) %in% unique(seqnames(clust_markers)),]
seqlengths(clustchrs) <- seqlengths(hg19Ideogram)[names(seqlengths(clustchrs))]

maskbed<-read.table(paste0(parentdir, "/reference_data/testcov.bed"), header=F, stringsAsFactors=F)
names(maskbed)<-c("CHR", "START", "END", "name", "gieStain", "V1", "V2", "V3", "pctmask")
maskbed<-maskbed %>% filter(!grepl("X|Y", CHR))
mask_markers <- GRanges(seqnames=maskbed$CHR,
                       ranges=IRanges(start=maskbed$START,
                         end=maskbed$END),
                        pctmask=maskbed$pctmask)
seqlengths(mask_markers) <- seqlengths(hg19Ideogram)[names(seqlengths(mask_markers))]

harbed<-read.table(paste0(parentdir, "/reference_data/2xHARs.hg19.sort.bed"), header=F, stringsAsFactors=F)
names(harbed)<-c("CHR", "START", "END", "ID", "VAL")
harbed<-harbed %>% filter(!grepl("X|Y", CHR))

# harbed <- harbed %>% arrange(desc(VAL)) %>% head(50)

har_markers <- GRanges(seqnames=harbed$CHR,
                       ranges=IRanges(start=harbed$START,
                         end=harbed$START),
                         ID=harbed$ID,
                         yval=15
                        )
seqlengths(har_markers) <- seqlengths(hg19Ideogram)[names(seqlengths(har_markers))]


p <- ggplot(clustchrs)+
  layout_karyogram(cytoband = TRUE)+
  geom_text(label=name)
  scale_fill_manual(values=c(giemsacols, stalk="#cd3333"))
p <- p+
  layout_karyogram(mask_markers,
    geom = "rect",
    ylim = c(10, 20),
    aes(alpha=pctmask, colour=NA), fill="black")+
  layout_karyogram(clust_markers,
    geom = "point",
    ylim = c(10, 20),
    aes(x=start, y=yval, color = cluster))+
  scale_colour_manual(values=myPaletteCatN(8)[4:8])+
  layout_karyogram(har_markers,
    geom = "point",
    ylim = c(0, 10),
    aes(x=start, y=yval, color = "black"))+
  guides(fill=FALSE,
    colour = guide_legend(override.aes=list(fill="white")))
p
ggsave(paste0(parentdir, "/images/chr20_cluster_karyo.png"), height=8, width=8)
# , size=6
# , xlim=c(-25, NA)


ggplot(chr20_atcg, aes(x=BIN, y=ERV_rel_rate, fill=Motif, colour=Motif, group=Motif))+
  geom_bar(stat="identity")+
  # geom_point()+
  # geom_line()+
  facet_wrap(~Motif)
ggsave(paste0(parentdir, "/images/chr20_rates_stacked.png"))


###############################################################################
# Analysis by 1Mb bin
###############################################################################
mct3a <- get_mct_b(bins1Mb)

# Count singletons per 3-mer subtype bin
bin_counts <- full_data$sites %>%
  mutate(Type=gsub("cpg_", "", Category2),
    SEQA=substr(Motif, cbp-i, cbp+i),
    SEQB=substr(Motif, cbp*3-i, cbp*3+i),
    Motif=paste0(SEQA, "(", SEQB, ")")) %>%
  group_by(CHR, BIN, Type, Motif) %>%
  summarise(n=n())

bin_counts <- merge(bin_counts, mct3a, by=c("CHR", "BIN", "Motif")) %>%
  mutate(ERV_rel_rate=n/nMotifs,
  subtype=paste0(Type, "_", Motif))


# Coerce from long to wide format
bins_wide <- bin_counts  %>%
  dplyr::select(CHR, BIN, subtype, ERV_rel_rate) %>%
  spread(subtype, ERV_rel_rate)

bins_wide[is.na(bins_wide)] <- 0

# Run PCA on bin matrix
bins_pca <- prcomp(scale(bins_wide[,-c(1:2)]))

# Cumulative variance explained by 96 PCs
# bin_cumvar <- cumsum((bins_pca$sdev)^2) / sum(bins_pca$sdev^2)

# Data frame containing bin IDs and components
bins_wide_pcs <- cbind(bins_wide[,c(1:2)], as.data.frame(bins_pca$x))

# Run NMF on bin matrix; get top signature for each bin
bins_nmf <- nmf(bins_wide[,-c(1:2)], 3, nrun=10)
bins_pred <- predict(bins_nmf, "rows", prob=T)

sigloads <- get_loads(bins_wide[,-c(1:2)], bins_nmf)
plot_loads(sigloads)
ggsave(paste0(parentdir, "/images/sigloads_bin.png"))

# Add predicted NMF class to data frame
# bins_wide_c <- cbind(bins_wide_c, sig=bins_pred[[1]])
bins_nmf_plotdat <- data.frame(bins_wide[,c(1:2)],
    basis(bins_nmf),
    sig=bins_pred[[1]]) %>%
  mutate(sum=X1+X2+X3,
    X1a=X1/sum,
    X2a=X2/sum,
    X3a=X3/sum)

# chrp_c <- merge(chrp_c, bins_nmf_plotdat, by=c("CHR", "BIN"))

bins_nmf_long <- bins_nmf_plotdat %>%
  gather(cluster, prob, X1a:X3a) %>%
  arrange(CHR, BIN) %>%
  group_by(cluster) %>%
  mutate(window=row_number(),
    CHR.BIN=paste0(CHR, "_", BIN),
    prob_m=rollmedian(prob, 5, fill=NA),
    # prob_m=prob,
    log_prob=log(prob_m)) %>%
  filter(log_prob > -5)

write.table(bins_nmf_long, "/net/csgsites/mutation/share/bin_sig_prob.txt", col.names=T, row.names=F, sep="\t", quote=F)
write.table(ind_nmf_long, "/net/csgsites/mutation/share/ing_sig_prob.txt", col.names=T, row.names=F, sep="\t", quote=F)

labs <- bins_nmf_long %>%
  group_by(CHR) %>%
  summarise(start=min(window),
    end=max(window),
    window=ceiling(start+(end-start)/2)) %>%
  mutate(fillcol=ifelse(CHR%%2==0, "even", "odd"))

iwhPalette5 <- c("#9555b4", "#98be57", "#514143", "#c7624e", "#92b5b5")
ggplot()+
  # facet_wrap(~CHR, ncol=1, scales="free_y", strip.position="right")+
  geom_rect(data=labs, aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=fillcol), alpha=0.5)+
  # geom_point(data=bins_nmf_long,
  #   aes(x=window, y=prob, colour=cluster, group=cluster), alpha=0.6)+
  geom_line(data=bins_nmf_long,
    aes(x=window, y=prob_m, colour=cluster, group=cluster), alpha=0.6)+
  # geom_bar(stat="identity")+
  # geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE)+
  # geom_smooth(method="loess", span=0.2)+
  scale_fill_manual(values=c("grey40", "grey50"), guide=FALSE)+
  scale_colour_manual(values=iwhPalette5)+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0), breaks=labs$window, labels=labs$CHR)+
  # ylab("log(cluster probability)")+
  ylab("signature contribution")+
  theme_bw()+
  # theme(axis.text.x=element_blank())+
  theme(legend.position="bottom",
    axis.title.x=element_blank())
ggsave(paste0(parentdir, "/images/by_bin2.png"), width=12, height=8)

bins_nmf_long %>%
  # filter(CHR==16) %>%
ggplot(aes(x=BIN, y=prob, colour=cluster, fill=cluster, group=cluster))+
  facet_wrap(~CHR, ncol=1, scales="free_y", strip.position="right")+
  # geom_point()+
  # geom_bar(stat="identity")+
  # geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE)+
  geom_smooth(method="loess", span=0.2)+
  theme_bw()+
  theme(axis.text.x=element_blank())+
  theme(legend.position="bottom")
ggsave(paste0(parentdir, "/images/by_bin.png"), width=12, height=18)

ggplot()+
  geom_point(data=bins_wide_c, aes(x=PC1, y=PC2),
    colour="black",pch=21, alpha=0.8)+
  theme_bw()+
  theme(axis.text.x=element_text(size=14),
    axis.text.y=element_text(size=14),
    axis.title.y = element_text(size=16),
    axis.title.x=element_text(size=16))+
  guides(colour=FALSE,
    fill = guide_legend(override.aes = list(size=6)))
ggsave(paste0(parentdir, "/images/bin_pca.png"), height=12, width=12)


clust_markers <- GRanges(seqnames=paste0("chr", pcd2$CHR),
                       ranges=IRanges(start=pcd2$BIN*binw-binw,
                         end=pcd2$BIN*binw),
                       sig=pcd2$sig,
                       yval=15)
data(hg19Ideogram, package = "biovizBase")
seqlengths(clust_markers) <- seqlengths(hg19Ideogram)[names(seqlengths(clust_markers))]
hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X", "Y")))

# clustchrs <- hg19
clustchrs<-hg19[seqnames(hg19) %in% unique(seqnames(clust_markers)),]
seqlengths(clustchrs) <- seqlengths(hg19Ideogram)[names(seqlengths(clustchrs))]

p <- ggplot(clustchrs)+
  layout_karyogram(cytoband = TRUE, label=name)+
  layout_karyogram(cytoband=FALSE, geom="text", label=name)
    # scale_fill_manual(values=c(giemsacols, stalk="#cd3333"))
p <- p+
  layout_karyogram(clust_markers,
    geom = "rect",
    size=6,
    ylim = c(10, 20),
    aes(x=start, color = factor(sig)))+
  scale_colour_manual(values=myPaletteCatN(8)[1:4])+
  guides(fill=FALSE,
    colour = guide_legend(override.aes=list(fill="white")))

ggsave(paste0(parentdir, "/images/nmf_cluster_karyo.png"), height=8, width=8)

# IDs PC1 vs PC2
ggplot()+
  geom_point(data=pcd2, aes(x=PC1, y=PC2),
         colour="black",pch=21, alpha=0.8)+
  theme_bw()+
  theme(axis.text.x=element_text(size=14),
    axis.text.y=element_text(size=14),
    axis.title.y = element_text(size=16),
    axis.title.x=element_text(size=16))+
  guides(colour=FALSE,
    fill = guide_legend(override.aes = list(size=6)))
ggsave(paste0(parentdir, "/images/id_pca.png"), height=12, width=12)
