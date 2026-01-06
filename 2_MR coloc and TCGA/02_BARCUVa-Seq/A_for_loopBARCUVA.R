sink('ColonSigmoidLoopResults.txt')

chooseCRANmirror(ind=92)
92


library(tidyr)
library(data.table)
library(dplyr)

# Load in BARCUVa-Seq eQTL data
colon <- fread("data/BARCUVa-Seq/barcuvaseq.eqtls.allpairs.sumstats.txt")

# Rename columns to match standard format
colon <- rename(colon, c("beta.exposure"="slope", "se.exposure"="slope_se",
                         "chr_name"="chr", "chrom_start"="pos", "eaf.exposure"="maf"))

# Standardize gene ID format (remove version numbers)
colon$gene_id <- substr(colon$gene_id, 1, 15)

# Load the same differential expression file as GTEx script
file <- fread("data/ABHD11_paper/Differentially_expressed_genes.csv")

# Remove duplicates
file <- distinct(file)

# Filter for genes present in BARCUVa-Seq data
colon_snp <- colon[colon$gene_id %in% file$ensembl_gene_id, ]

# Only keep genes which are in both data frames
file <- file[file$ensembl_gene_id %in% colon_snp$gene_id, ]

for (i in 1:nrow(file)){
  print(i)
  exposure_id <- file$ensembl_gene_id[i]
  
  print(exposure_id)
  
  # Extract data for current gene
  exposure_dat <- colon_snp[which(colon_snp$gene_id %in% file$ensembl_gene_id[i]), ]
  
  # Filter for cis region using chromosome and positions from file
  exposure_dat <- exposure_dat[exposure_dat$chr_name == as.numeric(file$chromosome_name[i])]
  exposure_dat <- exposure_dat[exposure_dat$chrom_start >= as.numeric(file$cis_start_position[i])]
  exposure_dat <- exposure_dat[exposure_dat$chrom_start <= as.numeric(file$cis_end_position[i])]
  
  # Keep SNPs associated at P < 5e-6
  exposure_dat <- exposure_dat[which(exposure_dat$pval_nominal <= 5e-8), ]
  
  if (nrow(exposure_dat>0)){
  
  # Write output file
  write.csv(exposure_dat, paste(exposure_id, "colonBARC_eqtl.csv", sep="_"))
  }
}