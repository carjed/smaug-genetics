#! /usr/bin/env Rscript

##########################################################################################
# script processes .ped file for downstream analysis
#
# need to reorganize or split into two scripts
# so the case/control files can be used right away 
# (without having to use singleton counts file)
##########################################################################################

##########################################################################################
#read ped file (may be unordered; also need to confirm bad samples are dropped)
# order ped file and merge with per-subject singleton counts; 
# subjects in full.singletons.count that are not in the ped are excluded
##########################################################################################
data2<-read.table("freeze4.20131216.v11.ped", header=TRUE)
ped<-data2[order(data2$IID) , ]
ped$Coverage.Int<-floor(ped$coverage)
ind<-sort(unique(ped$Coverage.Int))
ped_cases<-subset(ped, ped$Case.Control=="Case")
ped_controls<-subset(ped, ped$Case.Control=="Control")

##########################################################################################
#loop to match cases and controls by each coverage value
##########################################################################################
ped_out<-data.frame()
for (i in ind) {
	ped_cases_sub<-ped_cases[ped_cases$Coverage.Int==i,]
	ped_controls_sub<-ped_controls[ped_controls$Coverage.Int==i,]
	n<-nrow(ped_cases_sub)
	m<-nrow(ped_controls_sub)
	if (n==m) {
		ped_out<-rbind(ped_out, ped_cases_sub, ped_controls_sub)
	} else if (n>m) { 
		ped_cases_sub2<-ped_cases_sub[order(sample(ped_cases_sub$IID,m)),]
		ped_out<-rbind(ped_out, ped_cases_sub2, ped_controls_sub)
	} else if (m>n) {
		ped_controls_sub2<-ped_controls_sub[order(sample(ped_controls_sub$IID,n)),]
		ped_out<-rbind(ped_out, ped_cases_sub, ped_controls_sub2)
	}
}

##########################################################################################
#output indices for cases and controls to be used in BCFTools
##########################################################################################
write.table(ped_out, file="ped_out.ped", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(c(as.vector(ped_out[ped_out$Case.Control=="Case",]$IID), as.vector(ped_out[ped_out$Case.Control=="Control",]$IID)), file="subjects.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(as.vector(ped_out[ped_out$Case.Control=="Case",]$IID), file="cases.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(as.vector(ped_out[ped_out$Case.Control=="Control",]$IID), file="controls.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)




 
