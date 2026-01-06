
library(tidyr)
library(data.table)
library(dplyr)


colon<-fread("data/GTEx-Colon/GTExMetaAnalysisResultsFixed.csv")
file<-fread("data/ABHD11_paper/Differentially_expressed_genes.csv")
snps<-fread("data/GTEx-Colon/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")


colon<-colon[which(colon$pval<=5e-8),]
colon$phenotype_id<-substr(colon$phenotype_id,1,15)
colon<-colon[colon$phenotype_id %in% file$ensembl_gene_id,]
#separate variant_id column to get chr and SNP position
chromosome<-data.frame(colon$variant_id)
chromosome$colon.variant_id<- as.character(chromosome$colon.variant_id)
chromosome_1<-chromosome %>% separate(colon.variant_id, c("CHR", "BP", "other_allele.exposure", "effect_allele.exposure"))
colon_snp<-cbind(colon, chromosome_1)
colon_snp$CHR<- as.numeric(gsub("chr", "", colon_snp$CHR)) #remove the "chr" letters

#Get rsID
snps<-dplyr::select(snps,variant_id,rs_id_dbSNP151_GRCh38p7)
colon_snp<-merge(colon_snp,snps)

colon_snp<-dplyr::select(colon_snp,phenotype_id,rs_id_dbSNP151_GRCh38p7,CHR,BP,effect_allele.exposure,other_allele.exposure,beta,se,pval)

colnames(colon_snp)<-c("exposure","SNP","chromosome.exposure","position.exposure","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","pval.exposure")

file<-distinct(file)
file<-file[file$ensembl_gene_id %in% colon_snp$exposure,]

for (i in 1:nrow(file)){
  print(i)
  exposure_id<-file$ensembl_gene_id[i]
  
  print(exposure_id)
  
  exposure_dat<-colon_snp[which(colon_snp$exposure %in% file$ensembl_gene_id[i]),]
  
  exposure_dat<-exposure_dat[exposure_dat$chromosome.exposure==as.numeric(file$chromosome_name[i])]
  exposure_dat<-exposure_dat[exposure_dat$position.exposure>=as.numeric(file$cis_start_position[i])]
  exposure_dat<-exposure_dat[exposure_dat$position.exposure<=as.numeric(file$cis_end_position[i])]
  
  
  
  fwrite(exposure_dat, paste(exposure_id, "colontotal_eqtl.csv", sep="_"))
}



