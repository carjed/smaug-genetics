##############################################################################
# Function to get reverse complement
##############################################################################
revcomp <- function(DNAstr) {
	step1 <- chartr("ACGT","TGCA",DNAstr)
	step2 <- unlist(strsplit(step1, split=""))
	step3 <- rev(step2)
	step4 <- paste(step3, collapse="")
	return(step4)
}

##############################################################################
# Display plot with inset
##############################################################################
insetPlot <- function(main, inset, loc) {
	 print(main)
	 theme_set(theme_bw(base_size = 4))
	 print(inset, vp = loc)
	 theme_set(theme_bw())
 }

##############################################################################
# Get standard error estimate
##############################################################################
std <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

##############################################################################
# Compute standard error for correlations
##############################################################################
corSE<-function(corval, ct){
	sqrt((1-corval^2)/(ct-2))
}

##############################################################################
# Function to reverse sequence--used to correctly plot right flank in
# subsequence heatmaps
##############################################################################
reverse_chars <- function(string){
	string_split = strsplit(as.character(string), split = "")
	reversed_split = string_split[[1]][nchar(string):1]
	paste(reversed_split, collapse="")
}

##############################################################################
# Function to install and load packages
##############################################################################
usePackage <- function(p) {

	p <- as.character(substitute(p))

    if (!is.element(p, installed.packages()[,1])){
        install.packages(p, dep = TRUE)
		require(p, character.only = TRUE)
	} else {
		require(p, character.only = TRUE)
	}
}

##############################################################################
# Function to read file from disk
##############################################################################
make.data<-function(filename, chunksize, skiprows,...){
	conn<-NULL
	function(reset=FALSE){
		if(reset){
			if(!is.null(conn)) close(conn)
			conn<<-file(filename,open="r")
		} else{
			rval<-read.table(conn, nrows=chunksize, skip=skiprows,...)
			if ((nrow(rval)==0)) {
				close(conn)
				conn<<-NULL
				rval<-NULL
			}
			return(rval)
		}
	}
}

##############################################################################
# plot nbp heatmaps
##############################################################################
rrheat <- function(dat, f, levels, facetvar, nbp){
	p <- ggplot()+
	# log(v4*10000+1,2)
	# limits=c(min(dat$v4), max(dat$v4))
	geom_tile(data=dat, aes(x=v2a, y=v3, fill=v4))+
	# geom_text(data=dat, aes(x=v2a, y=v3, label=v4a, family="Courier", size=0.1))+
	geom_rect(data=f, size=0.8, colour="grey30",
		aes(xmin=xlo, xmax=xhi, ymin=ylo, ymax=yhi), fill=NA)+
	scale_fill_gradientn("Relative Rate\n",
		colours=myPalette((nbp-1)^4),
		trans="log",
		breaks=c(min(dat$v4), mean(dat$v4), max(dat$v4)),
		labels=c(round(min(dat$v4), 5),
	 		round(mean(dat$v4), 4),
			round(max(dat$v4), 3)),
		limits=c(min(dat$v4), max(dat$v4)))+
	xlab("5' flank")+
	ylab("3' flank")+
	theme(
		# legend.position="none",
		legend.title = element_text(size=18),
	  legend.text = element_text(size=16),
	  strip.text.x = element_text(size=20),
	  axis.title.x = element_text(size=20),
	  axis.title.y = element_text(size=20),
    axis.text.y = element_text(size=6, colour="black"),
	  axis.text.x = element_text(size=6, colour="black", angle=90, hjust=1))+
		# axis.text.y = element_blank(),
		# axis.text.x = element_blank())+
	scale_x_discrete(labels=levels)+
	facet_wrap(as.formula(paste("~", facetvar)), ncol=6, scales="free_x")

	return(p)
}

##############################################################################
# QQ plot in ggplot2 with qqline
##############################################################################
ggQQ <- function (vec) # argument: vector of numbers
{
  # following four lines from base R's qqline()
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]

  d <- data.frame(resids = vec)

  ggplot(d, aes(sample = resids)) +
	stat_qq() +
	geom_abline(slope = slope, intercept = int)

}

##############################################################################
# Negbin model processing functions
##############################################################################

# Run negbin model for each formula specified in formlist
runMod <- function(formlist, data){
	out <- lapply(formlist, function(x) glm.nb(x, data))
	names(out) <- names(formlist)
	return(out)
}

# Build list of fitted values (with CHR/BIN names) for each model
getFits <- function(modlist, data){
	out <- lapply(modlist,
		function(x){
			y <- fitted.values(x)
			names(y) <- paste0(data$CHR, ".", data$BIN)
			y
		})
	names(out) <- names(modlist)
	return(out)
}

# Build list of dataframes for each model
buildDF <- function(fitlist, data){
	out <- lapply(fitlist,
		function(x){
			data.frame(CHR, Category2=cat1, BIN,
				exp=x,
				obs=data$obs,
				# res=names(fitlist),
				stringsAsFactors=F)
		}
	)
	names(out) <- names(fitlist)
	return(out)
}

##############################################################################
# Append columns to windowed count data for all motif lengths
##############################################################################
getSubMotifs <- function(data, nts, b3){

	# nts <- ifelse(grepl("^AT", cat1), "A", "C")
	outdat <- data

	# Loop currently just runs for 7->5bp motif aggregation;
	# can run over 7->5->3 by setting last index to :1
	for(j in ((nbp-1)/2-1):2){

		# Specify iteration motif length
		mlength <- (j+1)*2+1

		# Define rule for substring evaluation
		# griddef <- paste(c(rep("bases", j), "nts", rep("bases", j)), collapse=",")

		griddef <- paste(c("bases", "bases", "nts", "b3", "bases"), collapse=",")

		# Evaluate substring rule and get vector of submotifs
		tris <- apply(eval(parse(text=paste("expand.grid(",griddef,")"))),
			1, paste, collapse="")

		# Loop through each substring and append column of
		# aggregated counts
		for(k in tris){
			# Generate regex string; j is fixed per iteration
			# (e.g., looking for internal 3-mers or 5-mers)
			# so we search for all 3-mers or 5-mers by allowing
			# any base preceding or following the internal motif
			# regtri <- paste0("^", "[A-Z]{", j, "}", i, "[A-Z]{", j, "}")
			regtri <- paste0("^[A-Z]", k, "[A-Z]")

			# Extract sequences matching this submotif
			z <- names(data)[grepl(regtri, names(data))]

			# Ensure motif match vector includes only sequences
			# corresponding to the appropriate motif length
			z <- z[nchar(head(gsub("_[A-Z]*", "", z)))==mlength]

			# Create column and append to df
			tripct <- data %>%
				dplyr::mutate_(.dots=setNames(paste(z, collapse="+"), k)) %>%
				dplyr::select_(.dots=k)
			outdat <- cbind(outdat, tripct)
		}
	}

	return(outdat)
}

##############################################################################
# Calculate standard error of auc curve
##############################################################################
seauc<-function(auc, n){
  q1<-auc/(2-auc)
  q2<-2*auc^2/(1+auc)
  num<-auc*(1-auc)+(n-1)*(q1-auc^2)+(n-1)*(q2-auc^2)
  denom<-n^2
  sqrt(num/denom)
}

##############################################################################
# Get recombination rate at each site
##############################################################################
rcrCol <- function(sites, file){
  feat_ranges <- bed_to_granges(file, header=T)
  site_ranges <- GRanges(seqnames=paste0("chr",sites$CHR),
                         ranges=IRanges(start=sites$POS, end=sites$POS))

  indices <- findOverlaps(site_ranges, feat_ranges, type="within", select="first")
  indices[is.na(indices)]<-0
  ind_df<-data.frame(POS=sites$POS, CHR=sites$CHR, indices)

  feat_df<-as.data.frame(feat_ranges)
  feat_df$indices<-seq_along(1:nrow(feat_df))
  rate_table <- merge(ind_df, feat_df, by="indices", all.x=T, incomparables=0) %>%
    arrange(CHR, POS)

  rates<-rate_table$id
  rates[is.na(rates)]<-0
  return(as.numeric(rates))
}

##############################################################################
# Get replication timing rate for each site
##############################################################################
repCol <- function(sites, repfile, binwidth){
  reptime <- read.table(repfile, header=F, stringsAsFactors=F, sep="\t")
  names(reptime) <- c("CHR", "END", "TIME")
  reptime2<-reptime %>%
    group_by(CHR) %>%
    arrange(CHR, END) %>%
    mutate(START=imputeStart(END))

  feat_ranges <- GRanges(seqnames=paste0("chr",reptime2$CHR),
                         ranges=IRanges(start=reptime2$START, end=reptime2$END), id=reptime2$TIME)
  site_ranges <- GRanges(seqnames=paste0("chr",sites$CHR),
                         ranges=IRanges(start=sites$POS, end=sites$POS))

  indices <- findOverlaps(site_ranges, feat_ranges, type="within", select="first")
  indices[is.na(indices)]<-0
  ind_df <- data.frame(POS=sites$POS, CHR=sites$CHR, indices)

  feat_df <- as.data.frame(feat_ranges)
  feat_df$indices <- seq_along(1:nrow(feat_df))
  rate_table <- merge(ind_df, feat_df, by="indices", all.x=T, incomparables=0) %>%
    arrange(CHR, POS)

  rates <- rate_table$id
  rates[is.na(rates)] <- 0
  return(as.numeric(rates))
  # return(as.integer(site_ranges %within% feat_ranges))
}

##############################################################################
# Check if site in bed file; returns 0 or 1
##############################################################################
binaryCol <- function(sites, bedfile){
  feat_ranges <- bed_to_granges(bedfile, header=F)
  site_ranges <- GRanges(seqnames=paste0("chr",sites$CHR),
                         ranges=IRanges(start=sites$POS, end=sites$POS))
  return(as.integer(site_ranges %within% feat_ranges))
}

##############################################################################
# Multiple plot function
#
# ggplot objects passed in ..., or to plotlist (as list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
##############################################################################
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##############################################################################
#' commandArgs parsing
#' return a named list of command line arguments
#'
#' Usage:
#' call the R script thus
#'   ./myfile.R --args myarg=something
#' or
#'   R CMD BATCH --args myarg=something myfile.R
#'
#' Then in R do
#'   myargs <- getArgs()
#' and myargs will be a named list
#' > str(myargs)
#' List of 2
#' $ file : chr "myfile.R"
#' $ myarg: chr "something"
#'
#' @title getArgs
#' @param verbose print verbage to screen
#' @param defaults a named list of defaults, optional
#' @return a named list
#' @author Chris Wallace
##############################################################################
getArgs <- function(verbose=FALSE, defaults=NULL) {
	myargs <- gsub("^--","",commandArgs(TRUE))
	setopts <- !grepl("=",myargs)
	if(any(setopts))
	myargs[setopts] <- paste(myargs[setopts],"=notset",sep="")
	myargs.list <- strsplit(myargs,"=")
	myargs <- lapply(myargs.list,"[[",2 )
	names(myargs) <- lapply(myargs.list, "[[", 1)

	## logicals
	if(any(setopts))
	myargs[setopts] <- TRUE

	## defaults
	if(!is.null(defaults)) {
		defs.needed <- setdiff(names(defaults), names(myargs))
		if(length(defs.needed)) {
		  myargs[ defs.needed ] <- defaults[ defs.needed ]
		}
	}

	## verbage
	if(verbose) {
		cat("read",length(myargs),"named args:\n")
		print(myargs)
	}
	myargs
}
