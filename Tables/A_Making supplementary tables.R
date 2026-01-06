library(dplyr)
library(readxl)
library(data.table)

gene_annotation_file <- "data/ABHD11_paper/Differentially_expressed_genes.csv"
genes <- fread(gene_annotation_file)
genes<-dplyr::select(genes,ensembl_gene_id,external_gene_name)


# Supplementary table 1 - differential expression analysis ----------------

st_1<-fread("Analysis/differential_expression_analysis_full_results.csv")

st_1<-dplyr::select(st_1,external_gene_name,Synonyms,ensembl_gene_id,ID,FC,`P.Value`)

colnames(st_1)<-c("Gene","Synonyms","Ensembl gene ID","Illumina microarray probe ID","Fold change","P-value")

st_1$`Fold change`<-round(st_1$`Fold change`,digits=2)
st_1$`P-value`<-signif(st_1$`P-value`,digits=3)

# Supplementary table 2 - F stats and r2 ----------------------------------

st_2<-fread("MR_Results/instrument_strength_summary.csv")

st_2<-merge(st_2,genes,by.x = "gene",by.y = "ensembl_gene_id",all.x=T)

st_2<-dplyr::select(st_2,gene,external_gene_name,dataset,n_snps,r2,F_statistic)

colnames(st_2)<-c("Gene","Ensembl gene ID","eQTL Dataset","Number of SNPs","r2","F-statistic")

st_2$r2<-round(st_2$r2,digits = 2)
st_2$`F-statistic`<-round(st_2$`F-statistic`,digits = 0)

st_2$`eQTL Dataset`<-gsub("GTExColon","GTEx",st_2$`eQTL Dataset`)

# Supplementary table 3 - harmonised summary genetic data -----------------

st_3<-fread("MR_Results/supplementary_harmonised_data.csv")

st_3<-merge(st_3,genes,by.x = "Gene",by.y = "ensembl_gene_id",all.x=T)

st_3<-dplyr::select(st_3,external_gene_name,Gene,Dataset,CRC_Subtype,SNP,Effect_Allele,Other_Allele,Beta_Exposure,SE_Exposure,Beta_Outcome,SE_Outcome)

colnames(st_3)<-c("Gene","Ensembl gene ID","eQTL Dataset","CRC Subtype","SNP", "Effect_Allele","Other_Allele","Beta_Exposure","SE_Exposure","Beta_Outcome","SE_Outcome")

st_3$`eQTL Dataset`<-gsub("GTExColon","GTEx",st_3$`eQTL Dataset`)

st_3$`CRC Subtype`<-stringr::str_to_title(st_3$`CRC Subtype`)


# Supplementary table 4 - MR results --------------------------------------

st_4<-fread("MR_Results/MR_combined_results.csv")
st_4$OR<-exp(st_4$b)

st_4$LCI<-exp(st_4$b-1.96*st_4$se)
st_4$UCI<-exp(st_4$b+1.96*st_4$se)

st_4<-dplyr::select(st_4,ILMN_Gene,ensembl_gene_id,type,id.outcome,method,nsnp,OR,LCI,UCI,pval,bh_p)

colnames(st_4)<-c("Gene","Ensembl gene ID","eQTL Dataset","CRC Subtype","Method", "Number of SNPs","Odds ratio","Lower 95% confidence interval","Upper 95% confidence interval","P-value","FDR-corrected P-value")

st_4$`eQTL Dataset`<-gsub("GTExColon","GTEx",st_4$`eQTL Dataset`)

st_4$`CRC Subtype`<-stringr::str_to_title(st_4$`CRC Subtype`)

st_4$`Odds ratio`<-round(st_4$`Odds ratio`,digits=2)
st_4$`Lower 95% confidence interval`<-round(st_4$`Lower 95% confidence interval`,digits=2)
st_4$`Upper 95% confidence interval`<-round(st_4$`Upper 95% confidence interval`,digits=2)
st_4$`P-value`<-signif(st_4$`P-value`,digits=3)
st_4$`FDR-corrected P-value`<-signif(st_4$`FDR-corrected P-value`,digits=3)


# Supplementary table 5 - genetic colocalisation results ------------------

st_5<-fread("Coloc/pwcoco_combined_results.csv")

st_5 <- st_5 %>%
  mutate(
    tissue_outcome = stringr::str_remove(analysis, "^ENSG[0-9]+_")
  ) %>%
  tidyr::separate(tissue_outcome, into = c("outcome", "type"), sep = "_", extra = "merge")

st_5$type <- sub("\\_coloc.coloc$", "", st_5$type)

st_5<-merge(st_5,genes,by.x="gene",by.y="ensembl_gene_id")

st_5<-dplyr::select(st_5,external_gene_name,gene,type,outcome,SNP1,SNP2,nsnps,H0,H1,H2,H3,H4)

colnames(st_5)<-c("Gene","Ensembl gene ID","eQTL Dataset","CRC Subtype","SNP1","SNP2","Number of SNPs","H0","H1","H2","H3","H4")

st_5$`eQTL Dataset`<-gsub("GTExColon","GTEx",st_5$`eQTL Dataset`)

st_5$`CRC Subtype`<-stringr::str_to_title(st_5$`CRC Subtype`)

st_5$H0<-round(st_5$H0,digits=2)
st_5$H1<-round(st_5$H1,digits=2)
st_5$H2<-round(st_5$H2,digits=2)
st_5$H3<-round(st_5$H3,digits=2)
st_5$H4<-round(st_5$H4,digits=2)


# Supplementary table 6 - TCGA analysis -----------------------------------

st_6 <- data.frame(
  Site = c(
    "Colon", "Colon", "Colon",
    "Rectum", "Rectum", "Rectum"
  ),
  Gene = c(
    "CHMP2A", "ABHD11", "ATP5MC2",
    "CHMP2A", "ABHD11", "ATP5MC2"
  ),
  Fold_change_normal_vs_tumour = c(
    1.28, 1.90, 1.46, 
    1.35, 7.87, 1.28
  ),
  P_value_normal_vs_tumour = c(
    1.41e-6, 4.39e-25, 2.38e-14, 
    2.58e-12, 1.47e-60, 1.31e-7
  )
)

colnames(st_6)<-gsub("_"," ",colnames(st_6))



# Put into Excel format ---------------------------------------------------

library(openxlsx)

wb <- createWorkbook()

tables <- list(
  "Supplementary table 1" = st_1,
  "Supplementary table 2" = st_2,
  "Supplementary table 3" = st_3,
  "Supplementary table 4" = st_4,
  "Supplementary table 5" = st_5,
  "Supplementary table 6" = st_6
)

for (sheet_name in names(tables)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, tables[[sheet_name]])
}

saveWorkbook(wb, file = "Supplementary_tables.xlsx", overwrite = TRUE)
