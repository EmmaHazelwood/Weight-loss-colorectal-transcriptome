sink('ColonLoopResults.txt')

chooseCRANmirror(ind=92)
92

library(tidyr)
library(data.table)
library(dplyr)
#BiocManager::install("biomaRt")
#n
library(biomaRt)

#Transverse
#load in colon cancer file
colon<-fread("ColonTransverse.csv")

#Get genes needed
file<-fread("Genesforloop.csv")

colon$gene<-substr(colon$phenotype_id,1,15)

colon<-colon[colon$gene %in% file$ensembl_gene_id]

#separate variant_id column to get chr and SNP position
chromosome<-data.frame(colon$variant_id)
chromosome$colon.variant_id<- as.character(chromosome$colon.variant_id)
chromosome_1<-chromosome %>% separate(colon.variant_id, c("CHR", "BP", "other_allele.exposure", "effect_allele.exposure"))
colon_snp<-cbind(colon, chromosome_1)
colon_snp$CHR<- as.numeric(gsub("chr", "", colon_snp$CHR)) #remove the "chr" letters
colon_snp$pos<-paste(colon_snp$CHR, colon_snp$BP,sep=":" ) #  
#rename column headings in prep for TwoSampleMR package

colon_snp<-rename(colon_snp, c("beta.exposure"="slope", "se.exposure"="slope_se",
                               "chr_name"= "CHR", "chrom_start"="BP", "eaf.exposure"="maf"))


colon_snp$gene<-substr(colon_snp$phenotype_id,1,15)
file<-file[file$ensembl_gene_id %in% colon_snp$gene]

for (i in 1:nrow(file)){
print(i)
exposure_dat<-colon_snp[which(colon_snp$gene %in% file[i,1]),]
exposure_dat$chr<-as.numeric(file[i,2]) # gtex cis SNP chr
exposure_dat$cis_start<-as.numeric(file[i,4]) # gtex cis SNP start position
exposure_dat$cis_end<-as.numeric(file[i,5]) # gtex cis SNP end position
#filter out for cis SNPs
exposure_dat_keep<-exposure_dat[which((exposure_dat$chr_name=exposure_dat$chr) & (exposure_dat$chrom_start <=exposure_dat$cis_end)& (exposure_dat$chrom_start >=exposure_dat$cis_start)),]
#Keep SNPs associated at P <0.000005
exposure_dat_keep<-exposure_dat_keep[which(exposure_dat_keep$pval_nominal<=0.000005),]
exposure_id<-file[i,1]
if (nrow(exposure_dat_keep)>0){
fwrite(exposure_dat_keep, paste("p0.000005/",exposure_id, "_colontransverse_eqtl.csv", sep="")) #adjust file name for types of colon cancer e.g. colon transverse or colon transverse otherwise you'll save over the other with just colon_eqtl
}}

#Sigmoid
#load in colon cancer file
colon<-fread("ColonSigmoid.csv")

#Get genes needed
file<-fread("Genesforloop.csv")

colon$gene<-substr(colon$phenotype_id,1,15)

colon<-colon[colon$gene %in% file$ensembl_gene_id]

#separate variant_id column to get chr and SNP position
chromosome<-data.frame(colon$variant_id)
chromosome$colon.variant_id<- as.character(chromosome$colon.variant_id)
chromosome_1<-chromosome %>% separate(colon.variant_id, c("CHR", "BP", "other_allele.exposure", "effect_allele.exposure"))
colon_snp<-cbind(colon, chromosome_1)
colon_snp$CHR<- as.numeric(gsub("chr", "", colon_snp$CHR)) #remove the "chr" letters
colon_snp$pos<-paste(colon_snp$CHR, colon_snp$BP,sep=":" ) #  
#rename column headings in prep for TwoSampleMR package

colon_snp<-rename(colon_snp, c("beta.exposure"="slope", "se.exposure"="slope_se",
                               "chr_name"= "CHR", "chrom_start"="BP", "eaf.exposure"="maf"))


colon_snp$gene<-substr(colon_snp$phenotype_id,1,15)
file<-file[file$ensembl_gene_id %in% colon_snp$gene]

for (i in 1:nrow(file)){
  print(i)
  exposure_dat<-colon_snp[which(colon_snp$gene %in% file[i,1]),]
  exposure_dat$chr<-as.numeric(file[i,2]) # gtex cis SNP chr
  exposure_dat$cis_start<-as.numeric(file[i,4]) # gtex cis SNP start position
  exposure_dat$cis_end<-as.numeric(file[i,5]) # gtex cis SNP end position
  #filter out for cis SNPs
  exposure_dat_keep<-exposure_dat[which((exposure_dat$chr_name=exposure_dat$chr) & (exposure_dat$chrom_start <=exposure_dat$cis_end)& (exposure_dat$chrom_start >=exposure_dat$cis_start)),]
  #Keep SNPs associated at P <0.000005
  exposure_dat_keep<-exposure_dat_keep[which(exposure_dat_keep$pval_nominal<=0.000005),]
  exposure_id<-file[i,1]
  if (nrow(exposure_dat_keep)>0){
    fwrite(exposure_dat_keep, paste("p0.000005/",exposure_id, "_colonsigmoid_eqtl.csv", sep="")) #adjust file name for types of colon cancer e.g. colon transverse or colon transverse otherwise you'll save over the other with just colon_eqtl
  }}

#Meta-analysis results
colon<-fread("GTEx-Colon/GTExMetaAnalysisResultsFixed.csv")
colon$gene<-substr(colon$phenotype_id,1,15)

#Get genes needed
file<-fread("Genesforloop.csv")

colon<-colon[colon$gene %in% file$ensembl_gene_id,]

#separate variant_id column to get chr and SNP position
chromosome<-data.frame(colon$variant_id)
chromosome$colon.variant_id<- as.character(chromosome$colon.variant_id)
chromosome_1<-chromosome %>% separate(colon.variant_id, c("CHR", "BP", "other_allele.exposure", "effect_allele.exposure"))
colon_snp<-cbind(colon, chromosome_1)
colon_snp$CHR<- as.numeric(gsub("chr", "", colon_snp$CHR)) #remove the "chr" letters
colon_snp$pos<-paste(colon_snp$CHR, colon_snp$BP,sep=":" ) #  
#rename column headings in prep for TwoSampleMR package
maf<-fread("GTEx-Colon/ColonSigmoid.csv")
maf<-dplyr::select(maf,variant_id,maf)
colon_snp<-merge(colon_snp,maf,by="variant_id")

colon_snp<-dplyr::rename(colon_snp, c("beta.exposure"="beta", "se.exposure"="se",
                               "chr_name"= "CHR", "chrom_start"="BP", "eaf.exposure"="maf"))

file<-file[file$ensembl_gene_id %in% colon_snp$gene]


for (i in 1:nrow(file)){
  print(i)
  exposure_dat<-colon_snp[which(colon_snp$gene == paste(file[i,1])),]
  exposure_dat$chr<-as.numeric(file[i,2]) # gtex cis SNP chr
  exposure_dat$cis_start<-as.numeric(file[i,4]) # gtex cis SNP start position
  exposure_dat$cis_end<-as.numeric(file[i,5]) # gtex cis SNP end position
  exposure_dat$chr_name<-as.numeric(exposure_dat$chr_name)
  #filter out for cis SNPs
  exposure_dat_keep<-exposure_dat[which((exposure_dat$chr_name==exposure_dat$chr) & (exposure_dat$chrom_start <=exposure_dat$cis_end)& (exposure_dat$chrom_start >=exposure_dat$cis_start)),]
  #Keep SNPs associated at P <0.000005
  exposure_dat_keep<-exposure_dat_keep[which(exposure_dat_keep$pval<=0.000005),]
  exposure_id<-file$ensembl_gene_id[i]
  if (nrow(exposure_dat_keep)>0){
    write.csv(exposure_dat_keep, paste(exposure_id, "colontotal_eqtl.csv", sep="_")) #adjust file name for types of colon cancer e.g. colon transverse or colon transverse otherwise you'll save over the other with just colon_eqtl
  }}



sink()