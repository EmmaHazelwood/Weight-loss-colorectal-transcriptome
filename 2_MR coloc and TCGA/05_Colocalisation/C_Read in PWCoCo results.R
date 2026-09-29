
library(data.table)

coloc_files <- list.files(pattern = "_coloc.coloc", full.names = TRUE)

# Keep only gene-outcome-dataset trios with MR evidence (FDR-P < 0.05) in the current MR results; older runs stay in the folder
mr_results <- fread("data/ABHD11_paper/MR_Results/MR_combined_results.csv")
sig_trios <- unique(mr_results[bh_p < 0.05, paste(ensembl_gene_id, id.outcome, type, sep = "_")])
trio <- sub("_coloc\\.coloc$", "", basename(coloc_files))
if (any(!trio %in% sig_trios)) print(paste("Skipping", sum(!trio %in% sig_trios), "PWCoCo results without current MR evidence:", paste(trio[!trio %in% sig_trios], collapse = ", ")))
if (any(!sig_trios %in% trio)) print(paste("No PWCoCo results yet for:", paste(setdiff(sig_trios, trio), collapse = ", ")))
coloc_files <- coloc_files[trio %in% sig_trios]

results_list <- lapply(coloc_files, function(file) {
  base_name <- gsub("_coloc$", "", basename(file))
  dt <- fread(file)
  dt[, analysis := base_name]
  return(dt)
})

combined_results <- rbindlist(results_list, fill = TRUE)

combined_results<-combined_results[rev(order(combined_results$H4)),]
combined_results$gene<-substr(combined_results$Dataset1,1,15)

fwrite(combined_results, "pwcoco_combined_results.csv")