#install.packages("meta")
#Get SNPs into one df
library(data.table)
library(dplyr)
sig<-fread("ColonSigmoid.csv")
sig$dataset<-"sigmoid"
trans<-fread("ColonTransverse.csv")
trans$dataset<-"transverse"

data<-rbind(sig,trans)
data$b<-data$slope
data$se<-data$slope_se
data$SNP<-data$variant_id
data<-select(data,phenotype_id,SNP,b,se,dataset)
write.csv(data,"ColonCombined.csv",row.names = FALSE,quote=FALSE,sep=",")

list<-data.frame(unique(data$phenotype_id))
list$file<-c(rep(1:100,times=260,each=1),1:7)
data<-merge(data,list,by.x="phenotype_id",by.y="unique.data.phenotype_id.")
list<-unique(data$file)

for (i in 1:length(list)){
genes<-data[which(data$file %in% list[i]),]
write.csv(genes,paste("genes",i,".csv",sep=""),row.names=FALSE,quote=FALSE,sep=",")
}
