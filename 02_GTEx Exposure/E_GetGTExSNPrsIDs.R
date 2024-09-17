sink("GTExSNPs.txt")

library(data.table)
library(dplyr)
library(tidyr)


ref<-fread("GTEx-Colon/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")
ref<-select(ref,variant_id,rs_id_dbSNP151_GRCh38p7)

#Transverse
paste("Transverse")
x<-dir(pattern="_colontransverse_eqtl.csv")
x<-substr(x,1,15)
x<-unique(x)

for (i in 1:length(x)){
  y<-dir(pattern=x[i])
  print(i)
  print(x[i])
  dat<-fread(grep("transverse_eqtl.csv",y,value=TRUE))
  dat<-merge(dat,ref)
  write.csv(dat, paste(x[i], "colontransverse_eqtl_snps.csv", sep="_"),row.names=FALSE)
}

#Sigmoid
paste("Sigmoid")
x<-dir(pattern="_colonsigmoid_eqtl.csv")
x<-substr(x,1,15)
x<-unique(x)

for (i in 1:length(x)){
  y<-dir(pattern=x[i])
  print(x[i])
  dat<-fread(grep("sigmoid_eqtl.csv",y,value=TRUE))
  dat<-merge(dat,ref)
  write.csv(dat, paste(x[i], "colonsigmoid_eqtl_snps.csv", sep="_"),row.names=FALSE)
}

#Total
paste("Total")
x<-dir(pattern="_colontotal_eqtl.csv")
x<-substr(x,1,15)
x<-unique(x)

for (i in 1:length(x)){
  y<-dir(pattern=x[i])
  print(x[i])
  dat<-fread(grep("total_eqtl.csv",y,value=TRUE))
  dat<-merge(dat,ref)
  write.csv(dat, paste(x[i], "colontotal_eqtl_snps.csv", sep="_"),row.names=FALSE)
}

sink()