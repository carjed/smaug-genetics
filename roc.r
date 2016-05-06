require(ggplot2)
require(dplyr)
require(tidyr)


##############################################################################
# Read and process data
##############################################################################
# chrp<-read.table("/net/bipolar/jedidiah/mutation/output/predicted/full/chr18_comb.txt", header=F)
# chrp<-read.table("/net/bipolar/jedidiah/mutation/output/predicted/chr18_1pct_mask.txt", header=F)

# Read data
cat("Reading data...\n")
chrpf<-read.table("/net/bipolar/jedidiah/mutation/output/predicted/full/rocdat_comb_3bp.txt", header=F)
names(chrpf)<-c("CHR", "POS", "MU", "OBS", "SEQ3", "MU3")

# Remove CpGs and sites with mu=0
chrpf<-chrpf[substr(chrpf$SEQ3, 2, 3)!="CG" & chrpf$MU>0,]

# Read DNMs
cat("Reading DNMs...\n")
dnms_full<-read.table("/net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs.txt", header=T, stringsAsFactors=F)
dnms_full<-dnms_full[,1:3]
names(dnms_full)<-c("ID", "CHR", "POS")

# Duplicate data, merge with DNMs to get ID
# chrpf<-chrp
cat("Annotating with ID...\n")
chrpf<-merge(chrpf, dnms_full, by=c("CHR", "POS"), all.x=T)
chrpf$ID[is.na(chrpf$ID)] <- "all"

# Subset to non-DNMs and DNMs
cat("Splitting by DNM status...\n")
chrpfa<-chrpf[chrpf$ID=="all",]
chrpfa<-chrpfa[sample(nrow(chrpfa), 1000000),]
chrpfdnm<-chrpf[chrpf$ID!="all",]

# Combine data
cat("Creating combined data...\n")
chrp<-rbind(chrpfdnm, chrpfa) %>% arrange(MU)
chrp$prop <- cumsum(chrp$OBS)/sum(chrp$OBS)


##############################################################################
# Permute data from the GoNL DNMs
# Each permutation simulates a single sample
# Will modify to extract median ROC curve for plotting
##############################################################################
nperm<-500

aucperm<-rep(0,nperm)
aucperm3<-aucperm
for(i in 1:nperm){
	cat("Permuting AUC", "(", i, "of", nperm, ")...\n")
	ndnms<-round(rnorm(1, 42.7, 10.3), 0)

	chrpsub<-rbind(chrpfdnm[sample(nrow(chrpfdnm), ndnms),], chrpfa) %>%
		arrange(MU)
	chrpsub$prop <- cumsum(chrpsub$OBS)/sum(chrpsub$OBS)

	nsamp<-50000

	chrpsub2<-chrpsub[sample(nrow(chrpsub), nsamp),] %>%
	  arrange(MU) %>%
	  mutate(ntile=ntile(MU, 1000))

	auctmp <- chrpsub2 %>% summarise(AUC=1-sum(prop)/nsamp)
	aucperm[i] <- auctmp$AUC

	# cat("Permuting 3-mer AUC", "(", i, "of", nperm, ")...\n")
	chrpsub3<-chrpsub %>% arrange(MU3)
	chrpsub3$prop <- cumsum(chrpsub3$OBS)/sum(chrpsub3$OBS)

	chrpsub3a<-chrpsub3[sample(nrow(chrpsub3), nsamp),] %>%
		arrange(MU3, prop) %>%
		mutate(ntile=ntile(MU3, 1000))

	auctmp3 <- chrpsub3a %>% summarise(AUC=1-sum(prop)/nsamp)
	aucperm3[i] <- auctmp3$AUC
}

##############################################################################
# Calculate AUC for each individual under the logit and 3-mer models
##############################################################################
ids<-unique(chrpfdnm$ID)
numind<-length(ids)
aucind<-rep(0, numind)
aucind3<-aucind

for(i in 1:numind){
	curid<-ids[i]
	cat("Calculating AUC for", curid, "(", i, "of", numind, ")...\n")
	datid<-chrpfdnm %>% filter(ID==curid)
	tmpdat<-rbind(datid, chrpfa) %>% arrange(MU)
	tmpdat$prop<-cumsum(tmpdat$OBS)/sum(tmpdat$OBS)

	nsamp<-50000

	tmpdat2<-tmpdat[sample(nrow(tmpdat), nsamp),] %>%
	  arrange(MU) %>%
	  mutate(ntile=ntile(MU, 1000))

	tmpauc <- tmpdat2 %>% summarise(AUC=1-sum(prop)/nsamp)
	aucind[i]<-tmpauc$AUC

	### Repeat with 3-mer model
	tmpdat3<-tmpdat %>% arrange(MU3)
	tmpdat3$prop<-cumsum(tmpdat3$OBS)/sum(tmpdat3$OBS)

	tmpdat3s<-tmpdat3[sample(nrow(tmpdat3), nsamp),] %>%
		arrange(MU3, prop) %>%
		mutate(ntile=ntile(MU3, 1000))

	tmpauc3 <- tmpdat3s %>% summarise(AUC=1-sum(prop)/nsamp)
	aucind3[i]<-tmpauc3$AUC
}

##############################################################################
# Function checks if elements in a exist in b
# Output is binary vector of length same as b
##############################################################################
toBin <- function(a,b){
	as.numeric(is.element(b,a))
}

##############################################################################
# Simulate mutations by randomly selecting sites and checking if
# Bernoulli(mu)==1; loop continues until we have simulated the same number
# of sites as in observed data
##############################################################################

# nsites <- sum(chrp$OBS)
nsim<-200

for(i in 1:nsim){
	cat("Running simulation", i, "of", nsim, "...\n")
	nsample <- 20000
	mutated <- c()
	nsites<-round(rnorm(1, 42.7, 10.3), 0)

	while(length(mutated) < nsites){
		rowind <- sample(nrow(chrp), nsample)
		row <- chrp[rowind,]
		mu <- row$MU

		batch <- sapply(mu, function(x) rbinom(1,1,x))
		mutated <- c(mutated, row$POS[which(as.logical(batch))])
	}

	# chrp$TMP <- toBin(mutated[1:nsites], chrp$POS)
	mus<-toBin(mutated, chrp$POS)
	nsimi<-sum(mus)
	chrp$last<-cumsum(mus)/nsimi

	colnames(chrp)[ncol(chrp)] <- paste0("simprop", i)
}

##############################################################################
# Subsample 100k sites, resort by mu, add percentile column, and coerce
# from wide to long format
#
# Can delete chrp after this step
##############################################################################
cat("Subsampling data...\n")

nsamp<-100000

chrp2<-chrp[sample(nrow(chrp), nsamp),] %>%
  arrange(MU) %>%
  mutate(ntile=ntile(MU, 1000)) %>%
  gather(group, val, .dots=grep("prop", names(chrp)))

auc <- chrp2 %>% group_by(group) %>% summarise(AUC=1-sum(val)/nsamp)

# chrp2<-chrp %>%
#   gather(group, val, .dots=grep("prop", names(chrp)))

cat("Cleaning up memory...\n")
# rm(chrp)
gc()

##############################################################################
# Further summarize data--get means by ntile
#
# Can delete chrp2 after this step
##############################################################################
chrp2_nt_mean<-chrp2 %>%
  group_by(group, ntile) %>%
  summarise(val=mean(val))

##############################################################################
# Get subsets of true and simulated data
##############################################################################
chrp2_nt_sim<-chrp2_nt_mean[chrp2_nt_mean$group!="prop",]
chrp2_nt_obs<-chrp2_nt_mean[chrp2_nt_mean$group=="prop",]

##############################################################################
# Calculate lower/upper 95% CIs from simulations
##############################################################################
chrp2_lb<-chrp2_nt_sim %>%
  group_by(ntile) %>%
  # summarise(lb=t.test(val)$conf.int[1]) %>%
	summarise(lb=max(val)) %>%
  mutate(group="lb") %>%
  select(group, ntile, lb)

chrp2_ub<-chrp2_nt_sim %>%
  group_by(ntile) %>%
  # summarise(ub=t.test(val)$conf.int[2]) %>%
	summarise(ub=min(val)) %>%
  mutate(group="ub") %>%
  select(group, ntile, ub)

##############################################################################
# Add bounds to data frame
##############################################################################
chrp2_nt_obs$lb <- chrp2_lb$lb
chrp2_nt_obs$ub <- chrp2_ub$ub
# names(chrp3m)[4:5]<-c("lb", "ub")

##############################################################################
# Repeat genome-wide with 3bp marginal rates
##############################################################################
# chrp3<-chrp %>% arrange(MU3)
# chrp3$prop <- cumsum(chrp3$OBS)/sum(chrp3$OBS)
#
# chrp3a<-chrp3[sample(nrow(chrp3), nsamp),] %>%
# 	arrange(MU3, prop) %>%
# 	mutate(ntile=ntile(MU3, 1000))
#
# chrp3b<-chrp3a %>%
#   group_by(ntile) %>%
#   summarise(val=mean(prop)) %>%
# 	mutate(group="3-mers") %>%
# 	dplyr::select(group, ntile, val)
#
# auc3 <- chrp3b %>% summarise(AUC=1-sum(val)/1000)

##############################################################################
# Combine AUC data to generate plots
##############################################################################

full_auc_dat<-rbind(data.frame(chrp2_nt_obs[,1:3]), data.frame(chrp3b))
names(full_auc_dat)[1]<-"Model"
full_auc_dat$Model<-ifelse(full_auc_dat$Model=="prop", "Logit", "3-mer")

full_auc_bounds<-chrp2_nt_obs[,c(2:5)]

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot()+
  geom_line(data=full_auc_dat,
		aes(x=(1000-ntile)/1000, y=1-val, group=Model, colour=Model), size=1.2)+
  geom_abline(intercept=0, slope=1)+
  geom_ribbon(data=full_auc_bounds,
		aes(x=(1000-ntile)/1000, y=1-val, ymin=1-ub, ymax=1-lb), alpha=0.2)+
	scale_colour_manual(values=cbbPalette[2:3])+
  coord_cartesian(xlim=c(0,1))+
	xlab("False Positive Rate")+
	ylab("True Positive Rate")+
  theme_bw()+
	theme(
		# legend.position="none",
		axis.title.x=element_text(size=20),
		axis.title.y=element_text(size=20))
  # scale_y_continuous(limits=c(0,1))+
  # scale_x_continuous(limits=c(0,1000), labels=c())
ggsave("/net/bipolar/jedidiah/mutation/images/pseudo_roc_curves.png", height=4, width=7.25)

auc$model<-"simulated"
auc$obs<-"simulated"

# aucind_df<-data.frame(AUC=aucind, gp="ind")
aucsim<-auc[-1,2:4]
auclogitobs<-data.frame(AUC=aucind, model="logit", obs="observed")
auc3merobs<-data.frame(AUC=aucind3, model="3-mer", obs="observed")
auclogitperm<-data.frame(AUC=aucperm, model="logit", obs="permuted")
auc3merperm<-data.frame(AUC=aucperm3, model="3-mer", obs="permuted")

plotaucfull<-rbind(aucsim, auclogitobs, auc3merobs, auclogitperm, auc3merperm)
plotaucnosim<-rbind(auclogitobs, auc3merobs, auclogitperm, auc3merperm)
plotaucobsonly<-rbind(auclogitobs, auc3merobs)
plotaucnoperm<-rbind(aucsim, auclogitobs, auc3merobs)

plotaucmean<-plotaucobsonly %>%
	group_by(model, obs) %>%
	summarise(AUC=mean(AUC))

plotaucobsonly$model<-relevel(plotaucobsonly$model, "3-mer")
plotaucfull$model<-relevel(plotaucfull$model, "3-mer")
plotaucnosim$model<-relevel(plotaucnosim$model, "3-mer")
plotaucnoperm$model<-relevel(plotaucnoperm$model, "3-mer")

ggplot()+
	geom_density(data=plotaucobsonly, aes(x=AUC, colour=model, linetype=obs))+
	geom_density(data=aucsim, aes(x=AUC), colour="black")+
	scale_colour_manual(values=cbbPalette[2:3])+
	geom_vline(data=plotaucmean,
		aes(xintercept=AUC, colour=model, linetype=obs))+
	# geom_vline(data=auc3merobs, aes(xintercept=mean(AUC)), linetype="longdash")+
	# stat_function(fun = dnorm,
	# 	colour = "red",
	# 	args = list(mean = mean(auc[-1,]$AUC), sd = sd(auc[-1,]$AUC)))+
	# xlim(mean(auc[-1,]$AUC)-4*sd(auc[-1,]$AUC), mean(auc[-1,]$AUC)+4*sd(auc[-1,]$AUC))+
	theme_classic()
ggsave("/net/bipolar/jedidiah/mutation/images/auc_hist.png", height=7, width=7)

aucind_df2<-data.frame(AUC=aucind, model="logit")
aucind3_df2<-data.frame(AUC=aucind3, model="3-mers")

aucind_comp<-rbind(aucind_df2, aucind3_df2)
aucind_comp$model<-relevel(aucind_comp$model, "3-mers")

plotaucfull$model<-as.factor(plotaucfull$model)
plotaucfull$model<-relevel(plotaucfull$model, "3-mer")

plotaucnoperm$model<-as.factor(plotaucnoperm$model)
plotaucnoperm$model<-relevel(plotaucnoperm$model, "3-mer")

ages<-read.table("/net/bipolar/jedidiah/mutation/reference_data/gonl_fam_age.txt", header=T)
names(ages)<-c("ID", "nDNM", "FatherAge", "MotherAge", "Coverage")
plotaucobsonly<-rbind(auclogitobs, auc3merobs)
plotaucobsonly$ID<-c(ids, ids)

plotaucobsonly<-merge(plotaucobsonly, ages, by="ID")

auc_age_cor<-plotaucobsonly %>%
	group_by(model) %>%
	summarise(cor=cor(AUC, FatherAge, method="spearman"),
		cor.p=cor.test(AUC, FatherAge, method="spearman")$p.value)

fao<-plotaucobsonly %>%
	filter(model=="logit") %>%
	filter(AUC>0.669 | AUC<0.597) %>%
	mutate(quart=ifelse(AUC>0.669, "yes", "no")) %>%
	dplyr::select(quart, FatherAge) %>%
	spread(quart, FatherAge)


ggplot(plotaucobsonly, aes(x=model, y=AUC, fill=model))+
	geom_violin()+
	geom_boxplot(width=0.3, fill="white", outlier.colour=NA)+
	geom_jitter(aes(colour=FatherAge), width=0.2, alpha=0.75)+
	geom_hline(yintercept=0.5, linetype="dashed")+
	# facet_wrap(~obs, ncol=1, drop=TRUE, scales="free_x")+
	# facet_grid(obs~., scales="free_x")+
	# scale_x_discrete(drop=TRUE)+
	# scale_y_discrete(drop=TRUE)+
	# scale_fill_discrete(drop=T)+
	# coord_flip()+
	# scale_colour_brewer(palette="Dark2")+
	scale_fill_manual(values=cbbPalette[c(2:3,1)], drop=TRUE)+
	# scale_colour_brewer()+
	scale_colour_gradientn(colours=brewer.pal(7,"PiYG"))+
	theme_classic()+
	theme(
		# legend.position="none",
	)
ggsave("/net/bipolar/jedidiah/mutation/images/ind_violin.png", width=6, height=4)

# ggplot(chrp2[chrp2$MU>0,], aes(x=log(MU)))+
# 	geom_histogram(bins=100)+
# 	xlab("log(relative rate)")+
# 	ylab("Frequency")+
# 	theme_bw()+
# 	theme(axis.text.y=element_blank(),
# 		axis.title.x=element_text(size=20),
# 		axis.title.y=element_text(size=20))
# ggsave("/net/bipolar/jedidiah/mutation/images/chr4_hist.png", width=10, height=3)
#
# ggplot(chrpm, aes(x=log(MU)))+
# 	geom_histogram(bins=32)+
# 	xlab("log(relative rate)")+
# 	ylab("Frequency")+
# 	theme_bw()+
# 	theme(axis.text.y=element_blank(),
# 		axis.title.x=element_text(size=20),
# 		axis.title.y=element_text(size=20))
# ggsave("/net/bipolar/jedidiah/mutation/images/chr4_hista.png")

# OLD VERSION--early attempt at pseudo-ROC curves
# ggplot(chrp3, aes(x=1000-ntile, y=1-val, group=group, colour=group))+
#   geom_line()
# ggsave("/net/bipolar/jedidiah/mutation/images/psuedo_roc_chr4a.png")


# OLD VERSION--uses true ROC calculations
# devtools::install_github("hadley/ggplot2")
# devtools::install_github("sachsmc/plotROC")
# library(plotROC)
#
# basicplot2<-ggplot(chrp2, aes(d = OBS, m = MU))+
# 	geom_roc()
#
# basicplot2+annotate("text", x = .75, y = .25,
#            label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2)))+
# 	style_roc()
# ggsave("/net/bipolar/jedidiah/mutation/images/chr18_roc2.png")
