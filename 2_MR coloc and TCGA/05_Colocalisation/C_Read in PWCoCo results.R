
library(data.table)

coloc_files <- list.files(pattern = "_coloc.coloc", full.names = TRUE)

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