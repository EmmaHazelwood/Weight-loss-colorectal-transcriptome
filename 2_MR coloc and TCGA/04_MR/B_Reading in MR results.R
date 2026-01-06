library(data.table)
library(tidyverse)

# Read all MR result files
print("Reading MR result files...")
tbl_fread <- list.files(pattern = "_MR_Results.csv$") %>% 
  map_df(~fread(.))
tbl_fread <- as.data.frame(tbl_fread)

print(paste("Total results loaded:", nrow(tbl_fread)))
print(paste("Unique genes:", length(unique(tbl_fread$id.exposure))))
print(paste("Outcome types:", paste(unique(tbl_fread$id.outcome), collapse = ", ")))
print(paste("Exposure types:", paste(unique(tbl_fread$type), collapse = ", ")))

# Load gene annotation file (if available)
# This should map Ensembl IDs to gene names
gene_annotation_file <- "data/ABHD11_paper/Differentially_expressed_genes.csv"

if (file.exists(gene_annotation_file)) {
  print("Loading gene annotations...")
  genes <- fread(gene_annotation_file)
  
  # Merge with gene annotations
  # Assumes gene file has columns: ensembl_gene_id, gene_name (or similar)
  tbl_fread <- merge(genes, tbl_fread, 
                     by.x = "ensembl_gene_id", 
                     by.y = "id.exposure", 
                     all.y = TRUE)
  print(paste("Merged with gene annotations:", nrow(tbl_fread), "results"))
} else {
  print("Gene annotation file not found, skipping gene name merge")
  print(paste("Expected location:", gene_annotation_file))
}

# Calculate Benjamini-Hochberg adjusted p-values
print("Calculating FDR-adjusted p-values...")

# Adjust p-values globally across all tests 
tbl_fread$bh_p <- p.adjust(tbl_fread$pval, method = "BH")

# Sort by adjusted p-value
tbl_fread <- tbl_fread[order(tbl_fread$bh_p), ]

fwrite(tbl_fread,"MR_combined_results.csv")


# Average F stats and r2 --------------------------------------------------

f_stats<-fread("instrument_strength_summary.csv")

min(f_stats$F_statistic)
max(f_stats$F_statistic)
mean(f_stats$F_statistic)
median(f_stats$F_statistic)

length(unique(f_stats$gene))
