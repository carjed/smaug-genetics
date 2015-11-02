require(speedglm)

make.data<-function(filename, chunksize,...){
	conn<-NULL
	function(reset=FALSE){
		if(reset){
			if(!is.null(conn)) close(conn)
			conn<<-file(filename,open="r")
		} else{
			rval<-read.table(conn, nrows=chunksize,...)
			if ((nrow(rval)==0)) {
				close(conn)
				conn<<-NULL
				rval<-NULL
			}
			return(rval)
		}
	}
}

fullfile<-"/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/full_sort.txt"
danames<-c("CHR", "BIN", "POS", "Sequence", "mut", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12")
da<-make.data(fullfile,chunksize=1000000,col.names=danames)

log_mod<-bigglm(mut~factor(Sequence)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12, data=da, family=binomial(), maxit=25)
log_mod<-shglm(mut~factor(Sequence)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12, datafun=da, family=binomial(), fitted=T)


da2<-read.table("/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/full_sort.txt", header=F, stringsAsFactors=F, nrows=1000000)
names(da2)<-danames

da3<-da2[sample(nrow(da2), 20000),]

da3<-read.table("/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/full_sort.txt", header=F, stringsAsFactors=F, nrows=500000)
names(da3)<-danames

log_mod_null<-speedglm(mut~Sequence+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12, data=da3, family=binomial(), fitted=T)
log_mod_null2<-glm(mut~Sequence+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12, data=da3, family=binomial)

# log_mod_null$coefficients<-log_mod$coefficients
# log_mod_null2$coefficients<-log_mod$coefficients

# original fitted vals from speedglm
lmnp<-head(plogis(log_mod_null$linear.predictors))

# original fitted vals from glm
lmnp2<-head(log_mod_null2$fitted.values)


lmnc<-head(log_mod_null$coefficients)
lmnc2<-head(log_mod_null$coefficients)

# replace coefficients with full model (returns NAs for speedglm object)
log_mod_null$coefficients<-log_mod$coefficients
log_mod_null3<-log_mod_null2
log_mod_null3$coefficients<-log_mod$coefficients

# preds<-plogis(predict(log_mod_null3, da3))

chunksize<-1000000
numsites<-1146613132
numchunks<-ceil(numsites/chunksize)

fullfile2<-"/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/AT_CG_full.txt"

for(i in 1:numchunks){

	tmp<-read.table(fullfile2, header=F, stringsAsFactors=F, nrows=chunksize, skip=((i-1)*chunksize))
	names(tmp)<-danames
	
	preds<-plogis(predict(log_mod_null3, tmp))
	
	out<-cbind(tmp[,1:3], preds)
	
	outpath<-paste0("/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/", i, "_test_out.txt")
	
	write.table(out, outpath, col.names=F, row.names=F, quote=F, sep="\t")
}

# predicted vals from manual matrix method
# build matrix
m <- model.matrix(mut~factor(Sequence)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12, da3)[,]

# get predicted vals
p2 <- coef(log_mod) %*% t(m)

# convert to probabilities
p2 <- plogis(as.vector(p2))
