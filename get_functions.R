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

# rrheat <- function(dat, facetvar){
	# p <- ggplot()+
	# geom_tile(data=dat, aes(x=v2a, y=v3, fill=log(v4*10000+1,2)))+
	# geom_text(data=dat, aes(x=v2a, y=v3, label=round(v4,3), family="Courier", size=0.1))+
	# geom_rect(data=f, size=1.4, colour="grey30", aes(xmin=xlo, xmax=xhi, ymin=ylo, ymax=yhi), fill=NA)+
	# scale_fill_gradientn(colours=myPalette((ncol(pc1)-1)/6))+
	# xlab("Left flank")+
	# ylab("Right flank")+
	# theme(legend.position="none")+
	# scale_x_discrete(labels=levs_a)+
	# facet_wrap(~facetvar, ncol=1)
	
	# return(p)
# }	
