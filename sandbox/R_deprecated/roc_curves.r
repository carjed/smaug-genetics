##############################################################################
# Get prediction curves under each model
##############################################################################
tmp <- 0
if(tmp){
	chrp_gfasc <- subByModel(chrp, "MU", "GFASC")
	head(chrp_gfasc) %>% data.frame()

	full_auc_dat <- subWrapper(chrp, sim=F)
	full_auc_dat %>%
		arrange(group, Category.x, ntile) %>%
		filter(Category.x=="GC_AT") %>%
		tail() %>%
		data.frame()

	outdat <- chrp %>%
		mutate_(mutmp="MU") %>%
		group_by(Category.x) %>%
		arrange(desc(mutmp)) %>%
		mutate(prop=cumsum(OBS)/sum(OBS)) %>%
		# dplyr::select(-mutmp) %>%
		# arrange_(mucol, "prop") %>%
		# mutate(mutmp=mucol, ntile=dense_rank(mutmp), group=groupname) %>%
		mutate(ntile=min_rank(desc(mutmp)), group="GFASC") %>%
		dplyr::select(-mutmp)
}

pred_curves <- 0
if(pred_curves){
	full_auc_dat <- subWrapper(chrp, sim=T)

	# Get table of AUC for each mutation type/model combo
	full_auc2 <- full_auc_dat %>%
	  group_by(group, Category.x) %>%
	  summarise(AUC=sum(prop)/n()) %>%
	  spread(Category.x, AUC)

	newrow <- c("Relative change in AUC", (full_auc2[5,-1]-0.5)/(full_auc2[3,-1]-0.5))
	full_auc2 <- data.frame(rbind(as.matrix(full_auc2), newrow))

	roc_plotdat <- full_auc_dat[sample(nrow(full_auc_dat), .1*nrow(full_auc_dat)),] %>%
		group_by(group, Category.x, ntile) %>%
		summarise(prop=min(prop), n=mean(n))

	# Plot everything together
	plotROC(roc_plotdat,
		paste0(parentdir, "/images/pseudo_roc_curves_all.png"))

	# Plot 3-mers only
	roc_plotdat_seq <- roc_plotdat %>%
		filter(grepl("^ERV 3", group))
	plotROC(roc_plotdat_seq,
		paste0(parentdir, "/images/pseudo_roc_curves_seq_3.png"))

	# Plot 3-mers vs 5-mers
	roc_plotdat_seq <- roc_plotdat %>%
		filter(grepl("^ERV 5|^ERV 3", group))
	plotROC(roc_plotdat_seq,
		paste0(parentdir, "/images/pseudo_roc_curves_seq_5_3.png"))

	# Plot 7-mers vs 5-mers vs 3-mers
	roc_plotdat_seq <- roc_plotdat %>%
		filter(grepl("^E", group))
	plotROC(roc_plotdat_seq,
		paste0(parentdir, "/images/pseudo_roc_curves_seq.png"))

	# Plot common vs A&V
	fig_a_dat <- roc_plotdat %>%
	  filter(grepl("^A|B", group))
	plotROC(fig_a_dat,
		paste0(parentdir, "/images/av_vs_common_roc.png"))

	# Plot common vs ERVs
	fig_b_dat <- roc_plotdat %>%
	  filter(grepl("^E|B", group))
	plotROC(fig_b_dat,
		paste0(parentdir, "/images/common_vs_erv_roc.png"))

	# Plot ERVs vs logit
	fig_k_dat <- roc_plotdat %>%
		filter(grepl("^G|^E", group))
	plotROC(fig_k_dat,
		paste0(parentdir, "/images/erv_vs_logit_roc.png"))

	# Plot ERVs vs logit vs simulated max
	fig_k2_dat <- roc_plotdat %>%
		filter(grepl("^G|^E|^m", group))
	plotROC(fig_k2_dat,
		paste0(parentdir, "/images/erv_vs_logit_max_roc.png"))
}

##############################################################################
# ROC on 7-mers with max contrast between ERVs and polymorphisms
##############################################################################
maxc_roc <- 0
if(maxc_roc){
	maxc <- read.table(paste0(parentdir, "/maxc_7bp.txt"), header=T, stringsAsFactors=F)
	maxc$Category <- gsub("cpg_", "", maxc$Category2)
	chrp_maxc <- chrp %>%
	  dplyr::filter(paste0(SEQ, ".", Category.x) %in%
									paste0(maxc$Sequence, ".", maxc$Category))

	maxc_auc_dat <- subWrapper(chrp_maxc, sim=F)

	maxc_auc <- maxc_auc_dat %>%
	  group_by(group, Category.x) %>%
		# filter(grepl("^B|^E", group)) %>%
	  summarise(AUC=sum(prop)/n()) %>%
	  spread(Category.x, AUC)

	newrow<-c("Relative change in AUC", (maxc_auc[2,-1]-0.5)/(maxc_auc[1,-1]-0.5))
	maxc_auc <- data.frame(rbind(as.matrix(maxc_auc), newrow))

	maxc_plotdat <- maxc_auc_dat[sample(nrow(maxc_auc_dat), .1*nrow(maxc_auc_dat)),]

	# fig_max_dat <- maxc_plotdat %>%
	# 	filter(grepl("^E|B", group))
	# Plot all models under max contrast
	fig_max_dat <- maxc_auc_dat %>%
		filter(grepl("^E|^B", group))
	plotROC(fig_max_dat, paste0(parentdir, "/images/pseudo_roc_curves_maxc.png"))
}

# Additional max contrast analyses, with features
maxc_extra1 <- 0
if(maxc_extra1){
	coefs <- read.table(
		paste0(parentdir, "/output/logmod_data/coefs/coefs_full.txt"),
		header=F, stringsAsFactors=F)

	names(coefs) <- c("Cov", "Est", "SE", "Z", "pval", "Sequence", "Category")

	# coefs$qval <- p.adjust(coefs$pval)
	splitCoefs <- function(x){
	  coefs %>%
	    # filter(Cov==x, pval<1e-08) %>%
			filter(Cov==x) %>%
			# mutate(qval=p.adjust(pval)) %>%
			# filter(qval<0.1) %>%
			mutate(ntile=ntile(pval, 100)) %>%
			filter(ntile<=10) %>%
			# filter(ntile<=25 | ntile>=75) %>%
			# filter(exp(Est)>1.5 | exp(Est)<0.67) %>%
	    dplyr::select(Sequence, Category)
	}

	coef_split <- lapply(unique(coefs$Cov), function(x) splitCoefs(x))
	names(coef_split) <- unique(coefs$Cov)

	matchMotifs <- function(x) {
	  chrp %>%
			# ungroup() %>%
			# mutate(Category.x=as.character(Category.x)) %>%
	    dplyr::filter((substr(SEQ,1,7) %in% x$Sequence) & (Category.x %in% x$Category))
	}

	maxc_feat <- data.frame(Sequence=plotdat$Sequence,
		Category=as.character(gsub("cpg_", "", plotdat$Category)))
	chrp_maxc_features<-matchMotifs(maxc_feat)
	maxc_auc_feat <- subWrapper(chrp_maxc_features, sim=F)
	cov_auc<-aucCalc(maxc_auc_feat)

	chrp_maxc_features <- lapply(maxc_feat, function(x) matchMotifs(x))

	maxc_auc_feat <- lapply(chrp_maxc_features, function(x) subWrapper(x, sim=F))
	# test <- bind_rows(maxc_auc_feat, .id="id")

	aucCalc <- function(x){
	  x %>%
	    dplyr::filter(grepl("^E|^G", group)) %>%
	    group_by(group, Category.x) %>%
	    summarise(AUC=1-sum(prop)/n())
	    # spread(Category.x, AUC)
	}

	cov_auc <- lapply(maxc_auc_feat, function(x) aucCalc(x))
	cov_auc_dat <- bind_rows(cov_auc, .id="id")
	cov_auc_dat <- cov_auc_dat[complete.cases(cov_auc_dat),]

	cov_auc_dat <- cov_auc_dat %>% filter(!grepl("Intercept|DP", id))

	ggplot(cov_auc_dat, aes(x=id, y=AUC, colour=group, fill=group, group=group))+
	  geom_bar(stat="identity", position="dodge")+
		# geom_point()+
	  facet_wrap(~Category.x)+
	  coord_cartesian(ylim=c(0.5, 0.85))+
		coord_flip()+
		# theme_bw()
	ggsave(paste0(parentdir, "/images/test_bar.png"))

	cov_erv <- cov_auc_dat %>% filter(grepl("^E", group))
	cov_logit <- cov_auc_dat %>% filter(!grepl("^E", group))
	cov_logit$diff <- (cov_logit$AUC-0.5)/(cov_erv$AUC-0.5)-1
	# t.test(cov_logit$diff)

	cov_logit$dir <- ifelse(cov_logit$diff<=0, "lower", "higher")
	ggplot(cov_logit, aes(x=id, y=diff, colour=dir, fill=dir))+
		geom_bar(stat="identity", position="dodge")+
		scale_fill_brewer(palette="Set1")+
		scale_colour_brewer(palette="Set1")+
		# geom_point()+
	  facet_wrap(~Category.x)+
	  coord_cartesian(ylim=c(-0.15, 0.15))+
		coord_flip()+
		geom_hline(yintercept=0, linetype="dashed")+
		# xlab("AUC_[]")+
		theme_bw()
	ggsave(paste0(parentdir, "/images/diff_bar.png"), height=6, width=12)
}

maxc_extra2 <- 0
if(maxc_extra2){
	maxc_sum <- maxc_auc_dat %>%
		# mutate(SEQ3=substr(SEQ,3,5))%>%
		# group_by(Category.x) %>%
		group_by(Category.x, SEQ) %>%
		summarise(n=sum(OBS)) %>%
		dplyr::filter(n>0)

	# chrp_maxc2 <- merge(maxc_sum[,1:2], chrp_maxc, by=c("Category.x"))
	chrp_maxc$SEQ3 <- substr(chrp_maxc$SEQ, 3, 5)

	chrp_maxc2 <- merge(maxc_sum[,1:2], chrp_maxc, by=c("Category.x", "SEQ"))
	# chrp_maxc2 <- merge(maxc_sum[,1:2],
		# chrp_maxc[!grepl("^gonl", chrp_maxc$ID),], by=c("Category.x", "SEQ"))
	maxc_auc_dat2 <- subWrapper(chrp_maxc2, sim=F)

	maxc_auc2 <- maxc_auc_dat2 %>%
	  # group_by(group, Category.x) %>%
		group_by(group, Category.x, SEQ) %>%
	  summarise(AUC=1-sum(prop)/n()) #%>%
	  #spread(Category.x, AUC)

	# maxc_auc2 %>% arrange(Category.x, group) %>% data.frame()
	maxc_auc2 %>% arrange(Category.x, SEQ, group) %>% head(20)
	maxc_auc2 %>% group_by(group) %>% summarise(AUC=mean(AUC))
	maxc_test<-maxc_auc2 %>%
		group_by(group, Category.x) %>%
		# summarise(AUC=mean(AUC)) %>%
		arrange(Category.x) %>%
		spread(group, AUC) %>%
		data.frame()

	names(maxc_test) <- c("Category.x", "SEQ", "AV", "BP", "ERV", "GFASC")
	maxc_test$diff <- maxc_test$ERV-maxc_test$BP
	maxc_test %>% group_by(Category.x)%>% summarise(diff=mean(diff), n=n())
	# not by seq
	# maxc_auc <- maxc_auc_dat %>%
	#   group_by(group, Category.x) %>%
	#   summarise(AUC=1-sum(prop)/n()) %>%
	#   spread(Category.x, AUC)
}

##############################################################################
# ROC by GoNL individual
##############################################################################
gonl_ind <- 0
if(gonl_ind){
	nperm<-258
	# aucperm<-rep(0,nperm)
	# aucperm3<-aucperm

	ids <- unique(chrpfdnm$ID)[grepl("gonl", unique(chrpfdnm$ID))]
	numind <- length(ids)
	# aucind<-rep(0, numind)
	aucind3 <- aucind

	auc_ind<-data.frame()
	for(i in 1:nperm){
		### Run permutations
		# cat("Permuting AUC", "(", i, "of", nperm, ")...\n")
		ndnms<-round(rnorm(1, 42.7, 10.3), 0)
		nsamp<-1000000

		curid<-ids[i]
		### Get empirical AUC
		cat("Calculating AUC for", curid, "(", i, "of", numind, ")...\n")
		datid<-chrpfdnm %>% filter(ID==curid)

		tmpdat<-rbind(datid, chrpfa) %>% arrange(MU)
		tmp_auc_dat <- subWrapper(tmpdat, sim=F)

		tmp_auc <- tmp_auc_dat %>%
		  group_by(group, Category.x) %>%
		  summarise(AUC=1-sum(prop)/n()) %>%
		  spread(Category.x, AUC)
		tmp_auc$ID <- curid

		auc_ind<-rbind(auc_ind, data.frame(tmp_auc))
	}


	# ggplot(auc_ind, aes(x=group, y=))
	# mapply(plotROC, maxc_auc_feat,
		# paste0(parentdir, "/images/", seq_along(max_auc_feat), "_roc.png"))
}
