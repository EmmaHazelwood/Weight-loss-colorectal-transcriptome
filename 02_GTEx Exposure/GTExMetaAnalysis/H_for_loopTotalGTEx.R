sink('ColonTotalLoopResults.txt')

chooseCRANmirror(ind=92)
92

library(tidyr)
library(data.table)
library(dplyr)
#load in colon cancer file
#colon<-fread("GTExMetaAnalysisResults.txt")
#colnames(colon)<-colon[1,]
#colon<-colon[which(colon$test=="fixed"),]
#fwrite(colon,"GTExMetaAnalysisResultsFixed.csv",row.names = FALSE,quote = FALSE)
colon<-fread("GTExMetaAnalysisResultsFixed.csv")
file<-fread("Genesforloop.csv")
colon<-colon[which(colon$pval<=0.000005),]
colon$phenotype_id<-substr(colon$phenotype_id,1,15)
file$id<-substr(file$id,1,15)
colon<-colon[colon$phenotype_id %in% file$id,]
#separate variant_id column to get chr and SNP position
chromosome<-data.frame(colon$variant_id)
chromosome$colon.variant_id<- as.character(chromosome$colon.variant_id)
chromosome_1<-chromosome %>% separate(colon.variant_id, c("CHR", "BP", "other_allele.exposure", "effect_allele.exposure"))
colon_snp<-cbind(colon, chromosome_1)
colon_snp$CHR<- as.numeric(gsub("chr", "", colon_snp$CHR)) #remove the "chr" letters
colon_snp$pos<-paste(colon_snp$CHR, colon_snp$BP,sep=":" ) #  
#rename column headings in prep for TwoSampleMR package

colon_snp<-rename(colon_snp, c("beta.exposure"="beta", "se.exposure"="se",
                             "chr_name"= "CHR", "chrom_start"="BP"))

file<-distinct(file)
file<-file[file$id %in% colon_snp$phenotype_id,]

for (i in 1:nrow(file)){
print(i)
exposure_dat<-colon_snp[which(colon_snp$phenotype_id %in% file[i,1]),]
exposure_dat$chr.exposure<-as.numeric(exposure_dat$chr.exposure)
exposure_dat$chrom_start<-as.numeric(exposure_dat$chrom_start)
exposure_dat$chr<-as.numeric(file[i,2]) #gtex chr
exposure_dat$cis_start<-as.numeric(file[i,4]) # gtex cis SNP start position
exposure_dat$cis_end<-as.numeric(file[i,5]) # gtex cis SNP end position
exposure_dat$cis_end<-as.numeric(exposure_dat$cis_end)
exposure_dat$cis_start<-as.numeric(exposure_dat$cis_start)
exposure_dat$chr<-as.numeric(exposure_dat$chr)
#filter out for cis SNPs
exposure_dat_keep<-exposure_dat[which((exposure_dat$chr_name=exposure_dat$chr) & (exposure_dat$chrom_start <=exposure_dat$cis_end)& (exposure_dat$chrom_start >=exposure_dat$cis_start)),]
exposure_id<-file[i,1]
fwrite(exposure_dat_keep, paste(exposure_id, "colontotal_eqtl.csv", sep="_")) #adjust file name for types of colon cancer e.g. colon sigmoid or colon transverse otherwise you'll save over the other with just colon_eqtl
}