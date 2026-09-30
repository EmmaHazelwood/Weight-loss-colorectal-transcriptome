library(tidyr)
library(data.table)
library(dplyr)

# BarcUVa-Seq eQTL summary statistics
colon <- fread("data/BARCUVa-Seq/barcuvaseq.eqtls.allpairs.sumstats.txt")

# maf is the minor allele frequency, carried as eaf.exposure but not used downstream
colon <- rename(colon, c("beta.exposure"="slope", "se.exposure"="slope_se",
                         "chr_name"="chr", "chrom_start"="pos", "eaf.exposure"="maf"))

# Remove Ensembl version suffixes
colon$gene_id <- substr(colon$gene_id, 1, 15)

# Differentially expressed genes
file <- fread("data/ABHD11_paper/Differentially_expressed_genes.csv")

file <- distinct(file)

colon_snp <- colon[colon$gene_id %in% file$ensembl_gene_id, ]

# Genes present in both
file <- file[file$ensembl_gene_id %in% colon_snp$gene_id, ]

for (i in 1:nrow(file)){
  print(i)
  exposure_id <- file$ensembl_gene_id[i]
  
  print(exposure_id)
  
  exposure_dat <- colon_snp[which(colon_snp$gene_id %in% file$ensembl_gene_id[i]), ]
  
  # cis region
  exposure_dat <- exposure_dat[exposure_dat$chr_name == as.numeric(file$chromosome_name[i])]
  exposure_dat <- exposure_dat[exposure_dat$chrom_start >= as.numeric(file$cis_start_position[i])]
  exposure_dat <- exposure_dat[exposure_dat$chrom_start <= as.numeric(file$cis_end_position[i])]
  
  # P < 5e-8
  exposure_dat <- exposure_dat[which(exposure_dat$pval_nominal <= 5e-8), ]
  
  if (nrow(exposure_dat) > 0){
  write.csv(exposure_dat, paste(exposure_id, "colonBARC_eqtl.csv", sep="_"))
  }
}