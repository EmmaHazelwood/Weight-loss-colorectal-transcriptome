library(dplyr)
library(readxl)
library(data.table)

gene_annotation_file <- "ABHD11_paper/Differentially_expressed_genes.csv"
genes <- fread(gene_annotation_file)
genes<-dplyr::distinct(dplyr::select(genes,ensembl_gene_id,external_gene_name))


# Supplementary table 1 - differential expression analysis ----------------

st_1<-fread("ABHD11 paper/Analysis/Results.csv")

#Annotate with Entrez gene ID
st_1_ref <- read.delim("Ensembl/GPL10558-50081.txt", quote="")

st_1<-merge(st_1,st_1_ref,by.x="Probe_Id",by.y="ID")

#Get Ensembl ID
# Entrez to Ensembl mapping from the Bioconductor annotation package, as the Ensembl BioMart service is unreliable;
# cached so reruns give the same mapping. Called via namespaces so AnnotationDbi::select does not mask dplyr::select.
bm_cache <- "Ensembl/entrez_to_ensembl_orgHsegdb.csv"
if (file.exists(bm_cache)) {
  bm <- fread(bm_cache)
} else {
  bm <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = as.character(unique(na.omit(st_1$Entrez_Gene_ID))),
                              columns = "ENSEMBL", keytype = "ENTREZID")
  bm <- as.data.table(bm)[!is.na(ENSEMBL)]
  setnames(bm, c("ENTREZID", "ENSEMBL"), c("entrezgene_id", "ensembl_gene_id"))
  bm[, entrezgene_id := as.integer(entrezgene_id)]
  fwrite(bm, bm_cache)
}

st_1<-merge(st_1,bm,by.x="Entrez_Gene_ID",by.y="entrezgene_id",all.x=T)

st_1$FC <- 2^st_1$logFC

st_1<-dplyr::select(st_1,Symbol.x,Synonyms,ensembl_gene_id,Probe_Id,FC,logFC,`P.Value`)

colnames(st_1)<-c("Gene","Synonyms","Ensembl gene ID","Illumina microarray probe ID","Fold change","log2 fold change","P-value")

st_1$`Fold change`<-round(st_1$`Fold change`,digits=2)
st_1$`log2 fold change`<-round(st_1$`log2 fold change`,digits=3)
st_1$`P-value`<-signif(st_1$`P-value`,digits=3)

st_1 <- st_1 %>%
  mutate(across(where(is.character), ~ na_if(.x, "")))

st_1<-st_1[order(st_1$Gene),]

# Supplementary table 2 - F stats and r2 ----------------------------------

st_2<-fread("MR_Results/instrument_strength_summary.csv")

st_2<-merge(st_2,genes,by.x = "gene",by.y = "ensembl_gene_id",all.x=T)

st_2<-dplyr::select(st_2,external_gene_name,gene,dataset,n_snps,r2,F_statistic)

colnames(st_2)<-c("Gene","Ensembl gene ID","eQTL Dataset","Number of SNPs","r2","F-statistic")

st_2$r2<-round(st_2$r2,digits = 2)
st_2$`F-statistic`<-round(st_2$`F-statistic`,digits = 0)

st_2$`eQTL Dataset`<-gsub("GTExColon","GTEx",gsub("BARC","BarcUVa-Seq",st_2$`eQTL Dataset`))

# Supplementary table 3 - heterogeneity between GTEx colon tissues --------

st_3<-fread("GTEx/GTEx_instrument_heterogeneity.csv")

# Keep only the GTEx instruments used in the MR; the GTEx folder also holds clumped files for genes not carried into the MR
mr_instruments <- unique(fread("MR_Results/supplementary_harmonised_data.csv")[Dataset == "GTExColon", .(Gene, SNP)])
st_3<-merge(st_3,mr_instruments,by.x = c("exposure","SNP"),by.y = c("Gene","SNP"))

st_3<-merge(st_3,genes,by.x = "exposure",by.y = "ensembl_gene_id",all.x=T)

st_3<-dplyr::select(st_3,external_gene_name,exposure,SNP,variant_id,beta,se,pval,I2)

colnames(st_3)<-c("Gene","Ensembl gene ID","SNP","GTEx variant ID","Beta (meta-analysis)","SE (meta-analysis)","P-value (meta-analysis)","I2 (%)")

st_3$`Beta (meta-analysis)`<-round(st_3$`Beta (meta-analysis)`,digits=2)
st_3$`SE (meta-analysis)`<-round(st_3$`SE (meta-analysis)`,digits=2)
st_3$`P-value (meta-analysis)`<-signif(st_3$`P-value (meta-analysis)`,digits=3)
st_3$`I2 (%)`<-round(st_3$`I2 (%)`,digits=1)

st_3<-st_3[order(st_3$Gene,st_3$SNP),]


# Supplementary table 4 - harmonised summary genetic data -----------------

st_4<-fread("MR_Results/supplementary_harmonised_data.csv")

st_4<-merge(st_4,genes,by.x = "Gene",by.y = "ensembl_gene_id",all.x=T)

st_4<-dplyr::select(st_4,external_gene_name,Gene,Dataset,CRC_Subtype,SNP,Effect_Allele,Other_Allele,Beta_Exposure,SE_Exposure,Beta_Outcome,SE_Outcome)

colnames(st_4)<-c("Gene","Ensembl gene ID","eQTL Dataset","CRC Subtype","SNP", "Effect_Allele","Other_Allele","Beta_Exposure","SE_Exposure","Beta_Outcome","SE_Outcome")

st_4$`eQTL Dataset`<-gsub("GTExColon","GTEx",gsub("BARC","BarcUVa-Seq",st_4$`eQTL Dataset`))

st_4$`CRC Subtype`<-stringr::str_to_title(st_4$`CRC Subtype`)

st_4$Beta_Exposure<-round(st_4$Beta_Exposure,digits=2)
st_4$SE_Exposure<-round(st_4$SE_Exposure,digits=2)
st_4$Beta_Outcome<-round(st_4$Beta_Outcome,digits=2)
st_4$SE_Outcome<-round(st_4$SE_Outcome,digits=2)


# Supplementary table 5 - MR results --------------------------------------

st_5<-fread("MR_Results/MR_combined_results.csv")
st_5$OR<-exp(st_5$b)

st_5$LCI<-exp(st_5$b-1.96*st_5$se)
st_5$UCI<-exp(st_5$b+1.96*st_5$se)

st_5<-dplyr::select(st_5,external_gene_name,ensembl_gene_id,type,id.outcome,method,nsnp,OR,LCI,UCI,pval,bh_p)
st_5<-dplyr::distinct(st_5)

colnames(st_5)<-c("Gene","Ensembl gene ID","eQTL Dataset","CRC Subtype","Method", "Number of SNPs","Odds ratio","Lower 95% confidence interval","Upper 95% confidence interval","P-value","FDR-corrected P-value")

st_5$`eQTL Dataset`<-gsub("GTExColon","GTEx",gsub("BARC","BarcUVa-Seq",st_5$`eQTL Dataset`))

st_5$`CRC Subtype`<-stringr::str_to_title(st_5$`CRC Subtype`)

st_5$`Odds ratio`<-round(st_5$`Odds ratio`,digits=2)
st_5$`Lower 95% confidence interval`<-round(st_5$`Lower 95% confidence interval`,digits=2)
st_5$`Upper 95% confidence interval`<-round(st_5$`Upper 95% confidence interval`,digits=2)
st_5$`P-value`<-signif(st_5$`P-value`,digits=3)
st_5$`FDR-corrected P-value`<-signif(st_5$`FDR-corrected P-value`,digits=3)


# Supplementary table 6 - genetic colocalisation results ------------------

st_6<-fread("Coloc/pwcoco_combined_results.csv")

tissue_outcome <- sub("^ENSG[0-9]+_", "", st_6$analysis)
st_6$outcome <- sub("_.*$", "", tissue_outcome)
st_6$type <- sub("_coloc\\.coloc$", "", sub("^[^_]+_", "", tissue_outcome))

st_6<-merge(st_6,genes,by.x="gene",by.y="ensembl_gene_id")

st_6<-dplyr::select(st_6,external_gene_name,gene,type,outcome,SNP1,SNP2,nsnps,H0,H1,H2,H3,H4)

colnames(st_6)<-c("Gene","Ensembl gene ID","eQTL Dataset","CRC Subtype","SNP1","SNP2","Number of SNPs","H0","H1","H2","H3","H4")

st_6$`eQTL Dataset`<-gsub("GTExColon","GTEx",gsub("BARC","BarcUVa-Seq",st_6$`eQTL Dataset`))

st_6$`CRC Subtype`<-stringr::str_to_title(st_6$`CRC Subtype`)

st_6$H0<-round(st_6$H0,digits=2)
st_6$H1<-round(st_6$H1,digits=2)
st_6$H2<-round(st_6$H2,digits=2)
st_6$H3<-round(st_6$H3,digits=2)
st_6$H4<-round(st_6$H4,digits=2)

st_6<-st_6[order(st_6$Gene,st_6$`eQTL Dataset`,st_6$`CRC Subtype`,-st_6$H4),]


# Supplementary table 7 - TCGA analysis -----------------------------------

st_7 <- data.frame(
  Site = c(
    "Colon", "Colon",
    "Rectum", "Rectum"
  ),
  Gene = c(
    "CHMP2A", "ABHD11",
    "CHMP2A", "ABHD11"
  ),
  Fold_change_normal_vs_tumour = c(
    1.28, 1.90,
    1.35, 7.87
  ),
  P_value_normal_vs_tumour = c(
    1.41e-6, 4.39e-25,
    2.58e-12, 1.47e-60
  )
)

colnames(st_7)<-gsub("_"," ",colnames(st_7))



# Put into Excel format ---------------------------------------------------

library(openxlsx)

wb <- createWorkbook()

tables <- list(
  "Supplementary table 1" = st_1,
  "Supplementary table 2" = st_2,
  "Supplementary table 3" = st_3,
  "Supplementary table 4" = st_4,
  "Supplementary table 5" = st_5,
  "Supplementary table 6" = st_6,
  "Supplementary table 7" = st_7
)

for (sheet_name in names(tables)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, tables[[sheet_name]])
}

saveWorkbook(wb, file = "ABHD11 paper/Supplementary_tables.xlsx", overwrite = TRUE)
