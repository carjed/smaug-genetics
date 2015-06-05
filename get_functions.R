# Function to get reverse complement
revcomp <- function(DNAstr) {
	step1 <- chartr("ACGT","TGCA",DNAstr)
	step2 <- unlist(strsplit(step1, split=""))
	step3 <- rev(step2)
	step4 <- paste(step3, collapse="")
	return(step4)
}

# Display plot with inset
insetPlot <- function(main, inset, loc) {
	 print(main)
	 theme_set(theme_bw(base_size = 4))
	 print(inset, vp = loc)
	 theme_set(theme_bw())
 }

# Get standard error estimate 
std <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

# Function to reverse sequence--used to correctly plot right flank in subsequence heatmaps
reverse_chars <- function(string){
	string_split = strsplit(as.character(string), split = "")
	reversed_split = string_split[[1]][nchar(string):1]
	paste(reversed_split, collapse="")
}

rrheat <- function(dat, levels, facetvar){
	p <- ggplot()+
	# log(v4*10000+1,2)
	# limits=c(min(dat$v4), max(dat$v4))
	geom_tile(data=dat, aes(x=v2a, y=v3, fill=v4))+
	# geom_text(data=dat, aes(x=v2a, y=v3, label=v4a, family="Courier", size=0.1))+
	geom_rect(data=f, size=1.4, colour="grey30", aes(xmin=xlo, xmax=xhi, ymin=ylo, ymax=yhi), fill=NA)+
	scale_fill_gradientn("Relative Rate\n", 
						 colours=myPalette((nbp-1)^4), 
						 trans="log", 
						 breaks=c(min(dat$v4), mean(dat$v4), max(dat$v4)),
						 labels=c(round(min(dat$v4), 5), round(mean(dat$v4), 4), round(max(dat$v4), 3)),
						 limits=c(min(dat$v4), max(dat$v4)))+
	xlab("5' flank")+
	ylab("3' flank")+
	theme(
		  # legend.position="none",
		  legend.title = element_text(size=18),
		  legend.text = element_text(size=16),
		  strip.text.x = element_text(size=40),
		  axis.title.x = element_text(size=20),
		  axis.title.y = element_text(size=20),
	      axis.text.y = element_text(size=16, colour="black"),
		  axis.text.x = element_text(size=16, colour="black"))+
	scale_x_discrete(labels=levels)+
	facet_wrap(as.formula(paste("~", facetvar)), ncol=1, scales="free_x")
	
	return(p)
}	

# QQ plot in ggplot2 with qqline
ggQQ <- function (vec) # argument: vector of numbers
{
  # following four lines from base R's qqline()
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]

  d <- data.frame(resids = vec)

  ggplot(d, aes(sample = resids)) + stat_qq() + geom_abline(slope = slope, intercept = int)

}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
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
