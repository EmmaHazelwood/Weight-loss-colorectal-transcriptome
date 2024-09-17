library(MRInstruments)
library(TwoSampleMR)
library(ggplot2)
library(dplyr)
library(data.table)
library(ieugwasr)


file<-read.csv("results/MiniProject2/01_Gene Expression/Genesforloop.csv")

file$id<-paste("eqtl-a-",file$ensembl_gene_id,sep="")

file1<-file[1:80,]
file2<-file[81:160,]
file3<-file[161:200,]
file3.5<-file[201:240,]
file4<-file[241:320,]
file5<-file[321:400,]
file6<-file[401:465,]

exposure_dat<-extract_instruments(outcomes=file1$id, p1=1e-06)
exposure_dat2<-extract_instruments(outcomes=file2$id, p1=1e-06)
exposure_dat3<-extract_instruments(outcomes=file3$id, p1=1e-06)
exposure_dat3.5<-extract_instruments(outcomes=file3.5$id, p1=1e-06)
exposure_dat4<-extract_instruments(outcomes=file4$id, p1=1e-06)
exposure_dat5<-extract_instruments(outcomes=file5$id, p1=1e-06)
exposure_dat6<-extract_instruments(outcomes=file6$id, p1=1e-06)
exposure_dat<-rbind(exposure_dat,exposure_dat2,exposure_dat3,exposure_dat3.5,exposure_dat4,exposure_dat5,exposure_dat6)

#Read in one gene that didn't work
#x<-fread("ENSG00000131475SNPs.csv")
#names(x)[names(x) == "CHROM"] <- "chr.exposure"
#names(x)[names(x) == "n"] <- "samplesize.exposure"
#names(x)[names(x) == "Effect_size"] <- "beta.exposure"
#names(x)[names(x) == "p"] <- "pval.exposure"
#names(x)[names(x) == "POS"] <- "pos.exposure"
#names(x)[names(x) == "SE"] <- "se.exposure"
#x$id.exposure <- "eqtl-a-ENSG00000131475"
#names(x)[names(x) == "ID"] <- "SNP"
#names(x)[names(x) == "REF"] <- "effect_allele.exposure"
#names(x)[names(x) == "ALT"] <- "other_allele.exposure"
#x$eaf.exposure <- 1-x$AAF
#x$exposure <- "|| id:eqtl-a-ENSG00000131475"
#x$mr_keep.exposure <- "TRUE"
#x$pval_origin.exposure <- "reported"
#x$data_source.exposure<-"igd"
#x<-select(x,chr.exposure,samplesize.exposure,beta.exposure,pval.exposure,pos.exposure,se.exposure,id.exposure,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,exposure,mr_keep.exposure,pval_origin.exposure,data_source.exposure)
#clump
#x_clumped<-ld_clump(data.frame(rsid=x$id.exposure, pval=x$pval.exposure))

#write.csv(x,file="ENSG00000131475.csv",row.names = FALSE,quote = FALSE)
#exposure_dat<-rbind(exposure_dat,x)

#Find gr37 start and end points
library(biomaRt)
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
gene = data.frame(file$ensembl_gene_id)
genes <- getBM(filters="ensembl_gene_id", attributes=c("ensembl_gene_id","start_position","chromosome_name"), values=gene, mart=ensembl)
genes$cisstart <- genes$start_position-1000000
genes$cisend <- genes$start_position+1000000

exposure_dat$ensembl_gene_id<-substr(exposure_dat$id.exposure,8,23)
  
test<-merge(exposure_dat, genes,by="ensembl_gene_id")
exposure_dat_keep<-test[which((test$chr.exposure==test$chromosome_name) & (test$pos.exposure<=test$cisend) & (test$pos.exposure>=test$cisstart)),]

list<-exposure_dat_keep$id.exposure
length(unique(list))

write.csv(exposure_dat_keep,file='eQTLGenGenes.csv',row.names = FALSE,quote=FALSE)

#Working out why 3.5 doesn't work
#for 3.5, 6 (eqtl-a-ENSG00000131475) doesn't seem to work so remove it

for (i in (1:nrow(file3.5))){
  fileid<-file3.5[i,]
  print(i)
  print(fileid$id)
  exposure_dat<-extract_instruments(outcomes=fileid$id, p1=1e-06)
}
