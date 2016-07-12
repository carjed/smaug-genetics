#! /usr/bin/env Rscript

##########################################################################################
# Read singleton counts file and processed .ped
##########################################################################################
scount<-read.table("full.singletons.count")
names(scount)<-c("Singletons", "IID")
data2<-read.table("ped_out.ped", header=TRUE)
ped<-merge(data2[order(data2$IID) , ], scount, by="IID")

##########################################################################################
# Boxplot for singleton counts by coverage
##########################################################################################
require(ggplot2)
ggplot(ped_out, aes(x=factor(Coverage.Int), y=Singletons))+geom_boxplot()+xlab("Coverage")
ggsave("singletons_vs_coverage.png")