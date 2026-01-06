library(data.table)
library(tidyverse)

# Set working directories
mr_dir <- "data/ABHD11_paper/MR_Results"
coloc_dir <- "data/ABHD11_paper/Coloc"

# Read MR results
cat("Reading MR results...\n")
mr_results <- fread(file.path(mr_dir, "MR_combined_results.csv"))

# Filter for significant MR results (FDR < 0.05)
mr_sig <- mr_results %>%
  filter(bh_p < 0.05)

cat(paste("Found", nrow(mr_sig), "significant MR results (FDR < 0.05)\n"))

# Read colocalization results
cat("Reading colocalization results...\n")
coloc_results <- fread(file.path(coloc_dir, "pwcoco_combined_results.csv"))

coloc_results <- coloc_results %>%
  mutate(
    tissue_outcome = str_remove(analysis, "^ENSG[0-9]+_")
  ) %>%
  separate(tissue_outcome, into = c("outcome", "type"), sep = "_", extra = "merge")

coloc_results$type <- sub("\\_coloc.coloc$", "", coloc_results$type)

coloc_results<-coloc_results[rev(order(coloc_results$H4)),]

coloc_results<-coloc_results %>% group_by(Dataset1, Dataset2) %>% slice(1)

# Read r2 and F statistics
#cat("Reading r2 and F statistics...\n")
#stats_results <- fread(file.path(mr_dir, "r2_and_F_stats_results.csv"))

# Remove any NA rows from stats
#stats_results <- stats_results %>%
#  filter(!is.na(Gene))

# Ensure numeric columns are properly formatted
#stats_results <- stats_results %>%
#  mutate(
#    r2 = as.numeric(r2),
#    `F stat` = as.numeric(`F stat`)
#  )

# Create the comprehensive summary table
cat("Creating summary table...\n")

# Join MR results with stats
summary_table <- mr_sig %>%
  # Extract dataset type from id.exposure or type column
  mutate(
    gene = ensembl_gene_id,
    dataset = case_when(
      grepl("GTEx", type, ignore.case = TRUE) ~ "GTEx",
      grepl("BARC", type, ignore.case = TRUE) ~ "BARC",
      grepl("eQTLGen", type, ignore.case = TRUE) ~ "eQTLGen",
      TRUE ~ type
    )
  )

# Join with stats
#summary_table <- summary_table %>%
#  left_join(
#    stats_results,
#    by = c("gene" = "Gene", "dataset" = "Dataset")
#  )

# Join with colocalization results
# Match on gene and outcome

summary_table<-merge(summary_table,coloc_results,by.x = c("ensembl_gene_id","id.outcome","type"), by.y = c("gene","outcome","type"))
  
# Select and order columns for the final summary
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
#    r2,
#    F_stat = `F stat`,
    H0,
    H1,
    H2,
    H3,
    H4
  ) 


summary_table<-summary_table[rev(order(summary_table$H4)),]

# Save full summary table
output_file <- file.path(mr_dir, "MR_Summary_with_Coloc_Stats.csv")
fwrite(summary_table, output_file)
cat(paste("\nSaved full summary table to:", output_file, "\n"))
