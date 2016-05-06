require(ggplot2)
require(dplyr)
require(tidyr)

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
# Calculate AUC for each individual
##############################################################################
ids<-unique(chrpfdnm$ID)
numind<-length(ids)
aucind<-rep(0, length(ids))

for(i in 1:numind){
	curid<-ids[i]
	cat("Calculating AUC for", curid, "(", i, "of", numind, ")...\n")
	datid<-chrpfdnm %>% filter(ID==curid)
	tmpdat<-rbind(datid, chrpfa) %>% arrange(MU)
	tmpdat$prop<-cumsum(tmpdat$OBS)/sum(tmpdat$OBS)

	nsamp<-100000

	tmpdat2<-tmpdat[sample(nrow(tmpdat), nsamp),] %>%
	  arrange(MU) %>%
	  mutate(ntile=ntile(MU, 1000))

	tmpauc <- tmpdat2 %>% summarise(AUC=1-sum(prop)/nsamp)
	aucind[i]<-tmpauc$AUC
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
nsites <- sum(chrp$OBS)
nsim<-200

for(i in 1:nsim){
	cat("Running simulation ", i, "of ", nsim, "...\n")
	nsample <- 50000
	mutated <- c()

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
chrp3<-chrp2 %>%
  group_by(group, ntile) %>%
  summarise(val=mean(val))

##############################################################################
# Get subsets of true and simulated data
##############################################################################
chrp3a<-chrp3[chrp3$group!="prop",]
chrp3m<-chrp3[chrp3$group=="prop",]

##############################################################################
# Calculate lower/upper 95% CIs from simulations
##############################################################################
chrp3lb<-chrp3a %>%
  group_by(ntile) %>%
  # summarise(lb=t.test(val)$conf.int[1]) %>%
	summarise(lb=max(val)) %>%
  mutate(group="lb") %>%
  select(group, ntile, lb)

chrp3ub<-chrp3a %>%
  group_by(ntile) %>%
  # summarise(ub=t.test(val)$conf.int[2]) %>%
	summarise(ub=min(val)) %>%
  mutate(group="ub") %>%
  select(group, ntile, ub)

##############################################################################
# Add bounds to data frame
##############################################################################
chrp3m$lb <- chrp3lb$lb
chrp3m$ub <- chrp3ub$ub
# names(chrp3m)[4:5]<-c("lb", "ub")

##############################################################################
# Repeat with 3bp marginal rates
##############################################################################
# tris<-0
# if(tris){
# 	chrpm <- read.table("/net/bipolar/jedidiah/mutation/chr4.rates.txt")
# 	names(chrpm) <- c("POS", "MU")
#
# 	chrpm$OBS <- toBin(chrdnms$POS, chrpm$POS)
# 	chrpm <- chrpm %>% arrange(MU)
# 	chrpm$prop <- cumsum(chrpm$OBS)/sum(chrpm$OBS)
#
# 	chrpm2<-chrpm[sample(nrow(chrpm), 1000000),] %>%
# 	  arrange(MU) %>%
# 	  mutate(ntile=ntile(MU, 1000), group="marginal") %>%
# 	  group_by(group, ntile) %>%
# 	  summarise(val=mean(prop))
#
# 	chrpm2$lb <- chrp3m$lb
# 	chrpm2$ub <- chrp3m$ub
#
# 	chrp3m <- rbind(chrp3m, chrpm2)
# }

ggplot()+
  geom_line(data=chrp3m[1:1000,],
		aes(x=(1000-ntile)/1000, y=1-val, group=group, colour=group), size=1.2)+
  geom_abline(intercept=0, slope=1)+
  geom_ribbon(data=chrp3m[1:1000,],
		aes(x=(1000-ntile)/1000, y=1-val, ymin=1-ub, ymax=1-lb), alpha=0.2)+
  coord_cartesian(xlim=c(0,1))+
	xlab("False Positive Rate")+
	ylab("True Positive Rate")+
  theme_bw()+
	theme(legend.position="none",
		axis.title.x=element_text(size=20),
		axis.title.y=element_text(size=20))
  # scale_y_continuous(limits=c(0,1))+
  # scale_x_continuous(limits=c(0,1000), labels=c())
ggsave("/net/bipolar/jedidiah/mutation/images/pseudo_roc_curves.png", height=4, width=7.25)

ggplot(auc[-1,], aes(x=AUC))+
	geom_density()+
	geom_vline(data=auc[1,], aes(xintercept=AUC), linetype="longdash")+
	stat_function(fun = dnorm,
		colour = "red",
		args = list(mean = mean(auc[-1,]$AUC), sd = sd(auc[-1,]$AUC)))+
	xlim(mean(auc[-1,]$AUC)-4*sd(auc[-1,]$AUC), mean(auc[-1,]$AUC)+4*sd(auc[-1,]$AUC))+
	theme_bw()

ggsave("/net/bipolar/jedidiah/mutation/images/auc_hist.png", height=7, width=7)

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
