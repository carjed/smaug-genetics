options(warn=-1)
sink("R.log")

{
	# Load packages
	suppressMessages(require(ggplot2))
	suppressMessages(require(plyr))
	suppressMessages(require(reshape2))
	suppressMessages(require(RColorBrewer))
	suppressMessages(require(grid))

	args<-commandArgs(TRUE)

	# Read args
	chr<-as.character(args[1])
	macl<-as.character(args[2])
	binw<-as.numeric(args[3])
	cpg_flag<-as.character(args[4])
	summ<-as.character(args[5])
	adj<-as.numeric(args[6])
	hot_flag<-as.character(args[7])
	imgdir<-args[8]
	bin1<-as.character(args[9])
	# bin2<-as.character(args[10])

	if (macl=="singletons") mac<-"Singleton"
	if (macl=="doubletons") mac<-"Doubleton"

	# Define palettes
	myPaletteCat <- colorRampPalette(rev(brewer.pal(8, "Dark2")), space="Lab")
	myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
	myPaletteB <- colorRampPalette(rev(brewer.pal(9, "Blues")), space="Lab")
	myPaletteR <- colorRampPalette(rev(brewer.pal(9, "Reds")), space="Lab")
	myPaletteG <- colorRampPalette(rev(brewer.pal(9, "Greens")), space="Lab")
	rb<-c(myPaletteB(6)[1:3],myPaletteR(6)[1:3])
	g<-myPaletteG(6)[1:3]
}