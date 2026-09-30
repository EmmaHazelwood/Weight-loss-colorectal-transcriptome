library(data.table)
library(tidyverse)

mr_dir <- "data/ABHD11_paper/MR_Results"
coloc_dir <- "data/ABHD11_paper/Coloc"

# Read MR results
cat("Reading MR results...\n")
mr_results <- fread(file.path(mr_dir, "MR_combined_results.csv"))

# Filter for significant MR results (FDR < 0.05)
mr_sig <- mr_results %>%
  filter(bh_p < 0.05)

cat(paste("Found", nrow(mr_sig), "significant MR results (FDR < 0.05)\n"))

# Colocalisation results
cat("Reading colocalisation results...\n")
coloc_results <- fread(file.path(coloc_dir, "pwcoco_combined_results.csv"))

coloc_results <- coloc_results %>%
  mutate(
    tissue_outcome = str_remove(analysis, "^ENSG[0-9]+_")
  ) %>%
  separate(tissue_outcome, into = c("outcome", "type"), sep = "_", extra = "merge")

coloc_results$type <- sub("\\_coloc.coloc$", "", coloc_results$type)

coloc_results<-coloc_results[rev(order(coloc_results$H4)),]

coloc_results<-coloc_results %>% group_by(Dataset1, Dataset2) %>% slice(1)

cat("Creating summary table...\n")

summary_table <- mr_sig %>%
  mutate(
    gene = ensembl_gene_id,
    dataset = case_when(
      grepl("GTEx", type, ignore.case = TRUE) ~ "GTEx",
      grepl("BARC", type, ignore.case = TRUE) ~ "BARC",
      grepl("eQTLGen", type, ignore.case = TRUE) ~ "eQTLGen",
      TRUE ~ type
    )
  )

# Add colocalisation results
summary_table<-merge(summary_table,coloc_results,by.x = c("ensembl_gene_id","id.outcome","type"), by.y = c("gene","outcome","type"))
  
summary_table <- summary_table %>%
  select(
    ensembl_gene_id,
    ILMN_Gene,
    dataset,
    id.outcome,
    method,
    nsnp,
    b,
    se,
    pval,
    bh_p,
    H0,
    H1,
    H2,
    H3,
    H4
  ) 


summary_table<-summary_table[rev(order(summary_table$H4)),]

output_file <- file.path(mr_dir, "MR_Summary_with_Coloc_Stats.csv")
fwrite(summary_table, output_file)
cat(paste("\nSaved full summary table to:", output_file, "\n"))
