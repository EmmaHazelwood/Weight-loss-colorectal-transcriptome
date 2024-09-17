sink('biomartScriptResults.txt')
library(data.table)
library(dplyr)
library(tidyr)

chooseCRANmirror(ind=92)
92


x<-dir(pattern="_colonsigmoid_eqtl.csv")
x<-substr(x,1,15)
x<-unique(x)

#First try with myvariant

SNPs=data.frame(matrix(ncol=1,nrow=0))

for (i in 1:length(x)){
  y<-dir(pattern=x[i])
  print(x[i])
  dat<-fread(grep("sigmoid_eqtl.csv",y,value=TRUE))
  dat<-select(dat,chr,chrom_start,effect_allele.exposure,other_allele.exposure)
  dat <- dat %>% dplyr::mutate(hgvsid = paste0("chr", chr, ":g.", chrom_start, effect_allele.exposure, ">", other_allele.exposure))
  dat <- dat %>% dplyr::mutate(hgvsid2 = paste0("chr", chr, ":g.", chrom_start, other_allele.exposure, ">", effect_allele.exposure))
  SNP1<-dat$hgvsid
  SNP2<-dat$hgvsid2
  dat<-append(SNP1,SNP2)
  dat<-data.frame(dat)
  SNPs<-rbind(SNPs,dat)
}
SNPs<-distinct(SNPs)
fwrite(SNPs,"SNPsformyvariant.csv")
SNPs<-fread("SNPsformyvariant.csv")

variants <- myvariant::getVariants(SNPs$dat, fields = "dbsnp.rsid",assembly="hg38")

one<-variants$query
two<-variants$dbsnp.rsid

df<-data.frame(one,two)
df<-na.omit(df)
df<-df %>% separate(one,c("A","B","C","D"))
df$C<-substr(df$C,1,nchar(df$C)-1)
df$A<- gsub("chr", "", df$A)
df$chrom_strand<-df$A
df$chrom_start<-df$C
df$refsnp_id<-df$two
df<-dplyr::select(df,chrom_strand,chrom_start,refsnp_id)
df$chrom_strand<-as.numeric(df$chrom_strand)
df$chrom_start<-as.numeric(df$chrom_start)
df$allele<-"NA"

#For any which aren't found use biomart
#Get list of SNPs needed
SNPs=data.frame(matrix(ncol=4,nrow=0))

for (i in 1:length(x)){
  y<-dir(pattern=x[i])
  print(x[i])
  dat<-fread(grep("sigmoid_eqtl.csv",y,value=TRUE))
  dat<-dplyr::select(dat,chr,chrom_start,effect_allele.exposure,other_allele.exposure)
  SNPs<-rbind(SNPs,dat,use.names=FALSE)
}
SNPs<-distinct(SNPs)

SNPs<-merge(SNPs,df,by.x=c("X1","X2"),by.y=c("chrom_strand","chrom_start"),all=TRUE)
colnames(SNPs)<-c("chr_name","start","effect_allele.exposure","other_allele.exposure","SNP","allele")

SNPs<-SNPs[-which(SNPs$SNP!="NA"),]
SNPs$end<-SNPs$start
SNPs<-dplyr::select(SNPs,chr_name,start,end)

library(biomaRt)
mart<-useMart(biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp")

results=data.frame(matrix(ncol=5,nrow=0))

for (i in seq (1,nrow(SNPs))){
  tryCatch({
  print(i)
  print(SNPs[i,2])
  ens=getBM(attributes = c('refsnp_id','allele','chrom_start','chrom_strand'), filters = c('chr_name','start','end'), values = as.list(SNPs[i,]), mart = mart,useCache = FALSE)
  results=rbind(results,ens)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

results<-results[,c(4,3,1,2)]
results<-rbind(results,df)

write.table(results, file = "biomartResults.txt", sep = "\t",
            row.names = F)


sink()