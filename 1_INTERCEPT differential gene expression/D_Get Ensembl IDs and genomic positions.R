library(data.table)
library(biomaRt)
library(tidyverse)

merged <- fread("Probes_with_Entrez_ID.csv")
# Ensembl release 115
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://e115.ensembl.org")
genes <- getBM(filters="entrezgene_id", 
               attributes=c("ensembl_gene_id","entrezgene_id","external_gene_name",
                            "chromosome_name","start_position","end_position"), 
               values=unique(na.omit(merged$Entrez_Gene_ID)), 
               mart=ensembl,
               useCache = FALSE)
merged_new <- merge(merged, genes, by.x="Entrez_Gene_ID", by.y="entrezgene_id")

#Save full results for supplementary table
genes_collapsed <- genes %>%
  group_by(entrezgene_id) %>%
  summarise(
    ensembl_gene_id    = paste(unique(ensembl_gene_id), collapse = ","),
    external_gene_name = paste(unique(external_gene_name), collapse = ","),
    chromosome_name    = paste(unique(chromosome_name), collapse = ","),
    start_position     = paste(unique(start_position), collapse = ","),
    end_position       = paste(unique(end_position), collapse = ","),
    .groups = "drop"
  )

merged_full <- merge(merged, genes_collapsed, by.x="Entrez_Gene_ID", by.y="entrezgene_id",all.x = T,allow.cartesian = F)
fwrite(merged_full,"differential_expression_analysis_full_results.csv")

# Exclude sex chromosomes and the MHC region (chr6:28,477,797-33,448,354, GRCh38)
filtered_genes <- merged_new[
  !chromosome_name %in% c("X", "Y") &
    !(chromosome_name == "6" & 
        start_position >= 28477797 & 
        end_position <= 33448354)
]

filtered_genes$cis_start_position<-filtered_genes$start_position-500000
filtered_genes$cis_end_position<-filtered_genes$end_position+500000


# Write filtered results
fwrite(filtered_genes, "Differentially_expressed_genes.csv")
