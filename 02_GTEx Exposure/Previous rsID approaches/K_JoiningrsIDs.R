sink('Colonjoining.txt')


library(data.table)
library(dplyr)

#Sigmoid
SNPs<-fread("biomartResultssigmoid.txt")

x<-dir(pattern="_colonsigmoid_eqtl.csv")
x<-substr(x,1,15)
x<-unique(x)

for (i in 1:length(x)){
  y<-dir(pattern=x[i])
  print(i)
  print(x[i])
  dat<-fread(grep("sigmoid_eqtl.csv",y,value=TRUE))
  dat<-merge(dat,SNPs,by.x=c("chr_name","chrom_start"),by.y=c("chrom_strand","chrom_start"))
  fwrite(dat, paste(x[i], "colonsigmoid_eqtl.csv",sep="_"))
}

#transverse
SNPs<-fread("biomartResultstransverse.txt")

x<-dir(pattern="_colontransverse_eqtl.csv")
x<-substr(x,1,15)
x<-unique(x)

for (i in 1:length(x)){
  y<-dir(pattern=x[i])
  print(i)
  print(x[i])
  dat<-fread(grep("transverse_eqtl.csv",y,value=TRUE))
  dat<-merge(dat,SNPs,by.x=c("chr_name","chrom_start"),by.y=c("chrom_strand","chrom_start"))
  fwrite(dat, paste(x[i], "colontransverse_eqtl.csv",sep="_"))
}

sink()