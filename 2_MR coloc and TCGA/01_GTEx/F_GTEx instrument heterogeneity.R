library(data.table)

# Heterogeneity (I2) between GTEx transverse and sigmoid colon for each GTEx instrument,
# taken from the fixed-effect meta-analysis in C_GTEx meta-analysis.R

unclumped <- rbindlist(lapply(list.files(pattern = "_colontotal_eqtl\\.csv$"), fread))
clumped <- rbindlist(lapply(list.files(pattern = "_clumped\\.csv$"), fread), fill = TRUE)

# GTEx variant IDs (chr_pos_ref_alt_b38)
unclumped[, variant_id := paste0("chr", chromosome.exposure, "_", position.exposure, "_",
                                 other_allele.exposure, "_", effect_allele.exposure, "_b38")]
instruments <- merge(clumped[, .(exposure, SNP, beta = beta.exposure, se = se.exposure, pval = pval.exposure)],
                     unique(unclumped[, .(exposure, SNP, variant_id)]),
                     by = c("exposure", "SNP"))

# Pre-filter the (~26 GB) meta-analysis file with grep rather than reading it into R
meta_file <- "data/GTEx-Colon/GTExMetaAnalysisResultsFixed.csv"
ids_file <- tempfile(fileext = ".txt")
writeLines(unique(instruments$variant_id), ids_file)
meta_cols <- names(fread(meta_file, nrows = 0))
meta <- fread(cmd = paste("LC_ALL=C grep -F -w -f", shQuote(ids_file), shQuote(meta_file)),
              header = FALSE, col.names = meta_cols)
meta <- meta[variant_id %in% instruments$variant_id, .(phenotype_id, variant_id, I2)]
meta[, phenotype_id := substr(phenotype_id, 1, 15)]

het <- merge(instruments, meta, by.x = c("exposure", "variant_id"), by.y = c("phenotype_id", "variant_id"))
if (nrow(het) < nrow(instruments)) warning(nrow(instruments) - nrow(het), " instruments not found in the meta-analysis results")

# meta stores I2 as a proportion
het[, I2 := 100 * I2]

fwrite(het[, .(exposure, variant_id, SNP, beta, se, pval, I2)], "GTEx_instrument_heterogeneity.csv")

print(summary(het$I2))
