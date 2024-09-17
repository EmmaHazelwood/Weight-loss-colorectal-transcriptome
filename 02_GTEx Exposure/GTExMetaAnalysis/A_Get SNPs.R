library(tidyr)
library(data.table)
#load in colon cancer file
colon<-fread("GTExOverall.csv")
#separate SNP column to get chr and SNP position
chromosome<-data.frame(colon$SNP)
chromosome$colon.SNP<- as.character(chromosome$colon.SNP)
chromosome_1<-chromosome %>% separate(colon.SNP, c("CHR", "BP", "other_allele.exposure", "effect_allele.exposure"))
colon_snp<-cbind(colon, chromosome_1)
colon_snp$CHR<- as.numeric(gsub("chr", "", colon_snp$CHR)) #remove the "chr" letters
colon_snp$pos<-paste(colon_snp$CHR, colon_snp$BP,sep=":" ) #
#rename column headings in prep for TwoSampleMR package
library(dplyr)
colon_snp<-rename(colon_snp, c("beta.exposure"="slope", "se.exposure"="slope_se",
                               "chr_name"= "CHR", "chrom_start"="BP", "eaf.exposure"="maf"))

file<-fread("Genesforloop.csv")
#Only get genes which are in both data frames
genes <- file$id
GTEx <- colon_snp$phenotype_id
keep <- intersect(genes,GTEx)
file <- unique(file[file$id %in% keep ,])
colon_snp <- unique(colon_snp[colon_snp$phenotype_id %in% keep ,])

for (i in 1:nrow(file)){
  print(i)
  exposure_dat<-colon_snp[which(colon_snp$phenotype_id %in% file[i,1]),]
  exposure_dat$chr<-file[i,2] #gtex chr
  exposure_dat$chr<-as.numeric(exposure_dat$chr)
  exposure_dat$cis_start<-file[i,4] # gtex cis SNP start position
  exposure_dat$cis_start<-as.numeric(exposure_dat$cis_start)
  exposure_dat$cis_end<-file[i,5] # gtex cis SNP end position
  exposure_dat$cis_end<-as.numeric(exposure_dat$cis_end)
  #filter out for cis SNPs
  exposure_dat_keep<-exposure_dat[which((exposure_dat$chr_name=exposure_dat$chr) & (exposure_dat$chrom_start <=exposure_dat$cis_end)& (exposure_dat$chrom_start >=exposure_dat$cis_start)),]
  #Keep SNPs associated at P <0.05
  exposure_dat_keep<-exposure_dat_keep[which(exposure_dat_keep$pval_nominal<=0.000001),]
  exposure_id<-file[i,1]
  write.csv(exposure_dat_keep, paste(exposure_id, "colonoverall_eqtl.csv", sep="_")) #adjust file name for types of colon cancer e.g. colon sigmoid or colon transverse otherwise you'll save over the other with just colon_eqtl
}

