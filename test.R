#Test for Version control in R
#Wei
#20170810

library(gplots)
library(ggplot2)

#read table
All<-read.table("../Transcriptome.csv",sep=",",header=T,na.strings = "NA")

#function
AVG_Gene_expression<-function(DF){
  NEWDF<-DF[,1:6]
  T0_WT<-apply(DF[,7:9],1,mean)
  T4_WT<-apply(DF[,10:13],1,mean)
  T8_WT<-apply(DF[,14:17],1,mean)
  T12_WT<-apply(DF[,18:21],1,mean)
  T16_WT<-apply(DF[,22:25],1,mean)
  T20_WT<-apply(DF[,26:29],1,mean)
  T0_KO<-apply(DF[,30:32],1,mean)
  T4_KO<-apply(DF[,33:36],1,mean)
  T8_KO<-apply(DF[,37:40],1,mean)
  T12_KO<-apply(DF[,41:43],1,mean)
  T16_KO<-apply(DF[,44:46],1,mean)
  T20_KO<-apply(DF[,47:50],1,mean)

  NEWDF<-data.frame(NEWDF,T0_WT,T4_WT,T8_WT,T12_WT,T16_WT,T20_WT,T0_KO,
                    T4_KO,T8_KO,T12_KO,T16_KO,T20_KO)
  return(NEWDF)
}

time_plot<-function(x){
  GE_value<-x
  cond<-rep(c("WT","KO"),each=6)
  Time<-rep(seq(from= 0, to= 20,by = 4),2)
  DF<-data.frame(GE_value,cond,Time)
  GGP<-ggplot(aes(y=GE_value,x=Time),data=DF)+geom_point()+geom_line()+facet_wrap(~cond)
  print(GGP)
}

#AVG for gene expression
AVG_ALL<-AVG_Gene_expression(All)

time_plot(as.numeric(AVG_ALL[167,7:18]))

AVG_ALL[167,]

TimeGE<-t(apply(AVG_ALL[, 7:18],1,function(x) x/sum(x)))
head(TimeGE)
plot(TimeGE[1,1:12])

dim(as.data.frame((as.matrix(AVG_ALL[, 7:18]))))

#rank the pvalue
WT_NOTKO<-which(AVG_ALL[,3]<0.05 & AVG_ALL[,5]>0.05)
KO_NOTWT<-which(AVG_ALL[,3]>0.05 & AVG_ALL[,5]<0.05)

WT_NOTKO_200<-WT_NOTKO[order(AVG_ALL[WT_NOTKO,5],decreasing = T)[1:200]]
KO_NOTWT_200<-KO_NOTWT[order(AVG_ALL[KO_NOTWT,3],decreasing = T)[1:200]]

WT_200<-order(AVG_ALL[,3])[1:200]
KO_200<-order(AVG_ALL[,5])[1:200]

length(WT_KO)

WT_KO<-Reduce(union,list(WT_200,KO_200,WT_NOTKO_200,KO_NOTWT_200))
length(WT_KO)

#K-means
GE_Cluster <- kmeans(TimeGE[c(WT_NOTKO_200),], 5, nstart = 5)
heatmap.2(TimeGE[c(WT_NOTKO_200),][which(GE_Cluster$cluster==5),],Colv = F,trace='none',
          col = blues9,scale = "row")
heatmap.2(TimeGE[c(WT_NOTKO_200),1:6],Colv = F,trace='none',
          col = blues9,scale = "row")
