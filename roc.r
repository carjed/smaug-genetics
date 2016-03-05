require(ggplot2)
require(dplyr)
require(tidyr)

# Read DNM data
dnms <- read.table("/net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs.txt", header=T)
chrdnms <- dnms[dnms$CHROM==4,]

# Read predicted mutation rate data
chrp <- read.table("/net/bipolar/jedidiah/mutation/output/predicted/chr4_full_mask.txt", header=F)
names(chrp) <- c("POS", "MU")

# Function checks if elements in a exist in b
# Output is binary vector of length same as b
toBin <- function(a,b){
	as.numeric(is.element(b,a))
}

chrp$OBS <- toBin(chrdnms$POS, chrp$POS)

# Sort by rate and add column of cumulate proportion
chrp <- chrp %>% arrange(MU)
chrp$prop <- cumsum(chrp$OBS)/sum(chrp$OBS)

# Simulate mutations by randomly selecting a site and checking if
# Bernoulli(mu)==1; loop continues until we have simulated the same number
# of sites as in observed data
nsites <- sum(chrp$OBS)

for(i in 1:3){
	TMP <- rep(0,nrow(chrp))
	k <- 0
  counter <- 0

	while(k <= nsites){

    counter <- counter+1

		rowind <- sample(nrow(chrp), 1)
		row <- chrp[rowind,]
		mu <- row$MU
		mutated <- rbinom(1,1,mu)

		if(mutated==1){
			k <- k+1
			TMP[rowind] <- 1
			#chrp[rowind,]$TMP<-1
			cat("Simulation ", i, ": ", k, " of ", nsites, "\n")
		}
	}

  # Get cumulative proportion of simulated data
	nsim <- sum(TMP)
	tmpprop <- cumsum(TMP)/nsim

  # Add columns to main df
	chrp <- cbind(chrp, TMP)
	chrp <- cbind(chrp, tmpprop)

  # Rename
	colnames(chrp)[ncol(chrp)-1] <- paste0("SIMOBS", i)
	colnames(chrp)[ncol(chrp)] <- paste0("simprop", i)
}

chrp2<-chrp[sample(nrow(chrp), 1000000),] %>%
  arrange(MU) %>%
  mutate(ntile=ntile(MU, 1000)) %>%
  gather(group, val, c(prop, simprop, simprop2, simprop3))

# Plot pseudo-ROC curves
ggplot(chrp2, aes(x=1000-ntile, y=1-val, group=group, colour=group))+
  geom_line()
ggsave("/net/bipolar/jedidiah/mutation/images/psuedo_roc_chr4.png")

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
