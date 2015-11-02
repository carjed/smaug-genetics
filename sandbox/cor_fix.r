source("mod_shell.r")
warnings()
mod.corr <- compare.all %>% 
filter(exp>100) %>%
group_by(Category2, res) %>%
summarise(num=length(exp), 
cor=cor(exp, obs, method="pearson"),
cor.p=cor.test(exp, obs, method="pearson")$p.value)
mod.corr$SE <- corSE(mod.corr$cor, mod.corr$num)
limits <- aes(ymax = mod.corr$cor + mod.corr$SE, ymin=mod.corr$cor - mod.corr$SE)
dodge <- position_dodge(width=0.9)
ggplot(mod.corr, aes(x=Category2, y=cor, fill=res))+
geom_bar(stat="identity", position=dodge)+
scale_colour_manual("Predictor", values=myPaletteCat(8)[5:8])+
scale_fill_manual("Predictor", values=myPaletteCat(8)[5:8])+
xlab("Category")+
ylab("Correlation with observed count")+
geom_errorbar(limits, position=dodge, width=0.25)+
theme_bw()+
theme(legend.title = element_text(size=18),
  legend.text = element_text(size=16),
  axis.title.x = element_text(size=20),
  axis.title.y = element_text(size=20),
  axis.text.y = element_text(size=16), 
  axis.text.x = element_text(size=16, angle = 45,  vjust=1, hjust=1.01))
modelbar <- paste0(parentdir, "/images/gw_5bp_vs_mod.png")  
ggsave(modelbar, width=7, height=7)
head(compare.all)
head(aggcat)
head(agg_5bp_100k)
head(summagg2)
sa2<-summagg2 %>% group_by(CHR, BIN, Category2) %>% summarise(AT=mean(AT), CG=mean(CG), obs=sum(obs, na.rm=T))
head(sa2)
187/446801
sa3<-sa2 %>% group_by(Category2) %>% summarise(AT=sum(AT), CG=sum(CG), obs=sum(obs))
head(sa3)
dim(sa3)
sa3
substr(sa3$Category2, nchar(sa3$Category2), -3)
substr(sa3$Category2, nchar(sa3$Category2)-3, nchar(sa3$Category2))
substr(sa3$Category2, nchar(sa3$Category2)-5, nchar(sa3$Category2)-3)
substr(sa3$Category2, nchar(sa3$Category2)-5, nchar(sa3$Category2)-2)
substr(sa3$Category2, nchar(sa3$Category2)-4, nchar(sa3$Category2)-3)
sa3$rel_prop<-ifelse(substr(sa3$Category2, nchar(sa3$Category2)-4, nchar(sa3$Category2)-3)=="AT", sa3$obs/sa3$AT, sa3$obs/sa3$CG)
sa3
5202194/1100439365
head(summagg2)
sa2<-summagg2 %>% 
group_by(CHR, BIN, Category2) %>% 
summarise(AT=mean(AT), CG=mean(CG), ct=sum(Count, na.rm=T), obs=sum(obs, na.rm=T))
sa3<-sa2 %>% 
group_by(Category2) %>% 
summarise(AT=sum(AT), CG=sum(CG), obs=sum(obs))
head(sa2)
sa2[1:9,]
head(sa3)
sa3<-sa2 %>% 
group_by(Category2) %>% 
summarise(ct=sum(ct), obs=sum(obs), rel_prop=obs/ct)
head(sa3)
sa3
head(agg_cov)
head(agg_5bp_100k)
head(s3)
head(sa3)
0.003290183*446791
head(agg_5bp_100k)
a3<-merge(agg_5bp_100k, sa3[,c(1,2,4)], by="Category2")
head(a3)
a3$exp1<-a3$nmotifs*a3$rel_prop
head(a3)
mc2<-a3 %>% group_by(Category2) %>% summarise(cor1=cor(exp, obs, method="pearson"), cor2=cor(exp1, obs, method="pearson"))
head(mc2)
mc2
head(compare.all)
head(mod.corr)
mod.corr
data.frame(mod.corr)
filter(mc2, CHR==1, BIN==1)
?filter
filter(mc2, CHR==1 & BIN==1)
filter(mc2, CHR==1 && BIN==1)
head(mc2)
filter(a3, CHR==1 && BIN==1)
dim(a3)
head(a3)
af<-filter(a3, CHR==1 & BIN==1)
dim(af)
head(af)
head(agg_5bp_100k)
a4<-agg_5bp_100k %>% group_by(Category2) %>% summarise(cor=cor(exp, obs, method="pearson"))
dim(a4)
a4
head(compare.all)
head(a3)
if(negbin_model){
cat("Initializing negbin regression model...\n")
ptm <- proc.time()
source("negbin_mod.r")
tottime<-(proc.time()-ptm)[3]
cat("Done. Finished in", tottime, "seconds \n")
}
head(agg_cov)
acv2<-agg_cov %>% group_by(Category2) %>% summarise(cor=cor(exp, obs, method="pearson")
)
head(acv2)
acv2<-agg_5bp_100k %>% group_by(Category2) %>% summarise(cor=cor(exp, obs, method="pearson"))
head(acv2)
head(mut_cov)
head(agg_cov)
head(agg_5bp_100k)
filter(agg_5bp_100k, CHR==10 && BIN==10)
filter(agg_cov, CHR==10 && BIN==10)
dim(agg_cov)
dim(agg_5bp_100k)
head(aggcat)
mut_cats
i<-6
cat1 <- mut_cats[i]
aggcat <- agg_cov[agg_cov$Category2==mut_cats[i],]
head(aggcat)
head(agg_5bp_100k)
filter(agg_5bp_100k, CHR==10 && BIN==10)
c1b1<-filter(agg_5bp_100k, CHR==10 && BIN==10)
c1b1<-filter(agg_5bp_100k, CHR==10 && BIN==10 && Category2=="GC_TA")
cor(c1b1$exp, c1b1$obs)
head(c1b1)
c1b1<-filter(agg_5bp_100k, CHR==10 & BIN==10 & Category2=="GC_TA")
head(c1b1)
c1b1<-filter(Category2=="GC_TA")
c1b1<-filter(agg_5bp_100k, Category2=="GC_TA")
head(c1b1)
cor(c1b1$exp, c1b1$obs)
head(agg_cov)
c1b2<-filter(agg_cov, Category2=="GC_TA")
head(c1b2)
cor(c1b2$exp, c1b2$obs)
dim(c1b1)
dim(c1b2)
head(c1b2)
head(c1b1)
c1b2<-arrange(c1b2, CHR, BIN)
head(c1b1)
head(c1b2)
c1b1[1:17,]
c1b2[1:17,1:8]
tail(c1b1)
tail(c1b2)
cor(c1b2$exp, c1b2$obs)

z1<-c1b2$exp
z2<-c1b2$obs
cor(z1,z2)
c1b1<-c1b1[,1:5]
c1b2<-c1b2[,1:5]
head(c1b1)
head(c1b2)
m2<-merge(c1b1, c1b2)
head(m2)
dim(m2)
cor(m2$exp, m2$obs)
fun.12 <- function(x.1,x.2,...){
     x.1p <- do.call("paste", x.1)
     x.2p <- do.call("paste", x.2)
     x.1[! x.1p %in% x.2p, ]
}
c1a<-fun.12(c1b1, c1b2)
dim(c1a)
head(c1a)
cor(c1a$exp, c1a$obs)
fun.12 <- function(x.1,x.2,...){
     x.1p <- do.call("paste", x.1)
     x.2p <- do.call("paste", x.2)
     x.1[x.1p %in% x.2p, ]
}
c1b<-fun.12(c1b1, c1b2)
dim(c1b)
cor(c1b$exp, c1b$obs)
420/2800
.15*.94+.85*.13
head(z1)
head(z2)
savehistory("cor_fix.r")
