sink('ColonSigmoidLoopResults.txt')

chooseCRANmirror(ind=92)
92

library(tidyr)
library(data.table)
#load in colon cancer file
colon<-fread("barcuvaseq.eqtls.allpairs.sumstats.txt")

library(dplyr)
colon<-rename(colon, c("beta.exposure"="slope", "se.exposure"="slope_se",
                               "chr_name"= "chr", "chrom_start"="pos", "eaf.exposure"="maf"))

colon$gene_id<-substr(colon$gene_id,1,15)

file<-fread("01_Gene Expression/Genesforloop.csv")
file$gene_id<-file$ensembl_gene_id
file<-distinct(file)
colon_snp <- merge(file, colon, by="gene_id")
#Only get genes which are in both data frames
genes <- file$gene_id
BARC <- colon_snp$gene_id
keep <- intersect(genes,BARC)
file <- unique(file[file$gene_id %in% keep ,])
colon_snp <- unique(colon_snp[colon_snp$gene_id %in% keep ,])

for (i in 1:nrow(file)){
  print(i)
  exposure_dat<-colon_snp[which(colon_snp$gene_id %in% file[i,1]),]
  exposure_dat$chr<-file[i,2] #BARC chr
  exposure_dat$chr<-as.numeric(exposure_dat$chr)
  exposure_dat$cis_start<-file[i,4] # BARC cis SNP start position
  exposure_dat$cis_start<-as.numeric(exposure_dat$cis_start)
  exposure_dat$cis_end<-file[i,5] # BARC cis SNP end position
  exposure_dat$cis_end<-as.numeric(exposure_dat$cis_end)
  #filter out for cis SNPs
  exposure_dat_keep<-exposure_dat[which((exposure_dat$chr_name=exposure_dat$chr) & (exposure_dat$chrom_start <=exposure_dat$cis_end)& (exposure_dat$chrom_start >=exposure_dat$cis_start)),]
  #Keep SNPs associated at P <0.0000005
  exposure_dat_keep<-exposure_dat_keep[which(exposure_dat_keep$pval_nominal<=0.000005),]
  exposure_id<-file[i,1]
  write.csv(exposure_dat_keep, paste(exposure_id, "colonBARC_eqtl.csv", sep="_")) #adjust file name for types of colon cancer e.g. colon sigmoid or colon transverse otherwise you'll save over the other with just colon_eqtl
}
