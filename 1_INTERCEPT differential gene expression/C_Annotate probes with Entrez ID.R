
library(dplyr)
library(data.table)
library(biomaRt)

res <- fread("Differentially_expressed_probes.csv")

ref <- read.delim("data/Ensembl/GPL10558-50081.txt", quote="")

res$ID <- res$Probe_Id
res$FC <- exp(res$logFC)

mergeall <- merge(ref, res, by="ID")

merged <- dplyr::select(mergeall, ID, ILMN_Gene, Entrez_Gene_ID, Chromosome, Synonyms, FC, P.Value)

write.csv(merged, "Probes_with_Entrez_ID.csv", row.names = FALSE, quote = FALSE)
