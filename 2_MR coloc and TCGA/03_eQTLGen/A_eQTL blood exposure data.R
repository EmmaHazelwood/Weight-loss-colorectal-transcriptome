library(MRInstruments)
library(TwoSampleMR)
library(ggplot2)
library(dplyr)
library(data.table)
library(ieugwasr)

file<-fread("data/ABHD11_paper/Differentially_expressed_genes.csv")

file$id<-paste("eqtl-a-",file$ensembl_gene_id,sep="")

file1<-file[1:80,]
file2<-file[81:160,]
file3<-file[161:200,]
file4<-file[201:nrow(file),]

exposure_dat<-extract_instruments(outcomes=file1$id, p1=5e-8)
exposure_dat2<-extract_instruments(outcomes=file2$id, p1=5e-8)
exposure_dat3<-extract_instruments(outcomes=file3$id, p1=5e-8)
exposure_dat4<-extract_instruments(outcomes=file4$id, p1=5e-8)

exposure_dat<-rbind(exposure_dat,exposure_dat2,exposure_dat3,exposure_dat4)

#Find gr37 start and end points
library(biomaRt)
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
gene = data.frame(file$ensembl_gene_id)
genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","start_position","chromosome_name"), values=gene, mart=ensembl,useCache = FALSE)
genes$cisstart <- genes$start_position-500000
genes$cisend <- genes$start_position+500000

exposure_dat$ensembl_gene_id<-substr(exposure_dat$id.exposure,8,23)
  
test<-merge(exposure_dat, genes,by="ensembl_gene_id")
exposure_dat_keep<-test[which((test$chr.exposure==test$chromosome_name) & (test$pos.exposure<=test$cisend) & (test$pos.exposure>=test$cisstart)),]

list<-exposure_dat_keep$id.exposure
length(unique(list))

write.csv(exposure_dat_keep,file='eQTLGenGenes.csv',row.names = FALSE,quote=FALSE)

