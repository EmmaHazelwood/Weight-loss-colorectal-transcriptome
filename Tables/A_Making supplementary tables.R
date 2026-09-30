library(dplyr)
library(readxl)
library(data.table)

gene_annotation_file <- "ABHD11_paper/Differentially_expressed_genes.csv"
genes <- fread(gene_annotation_file)
genes<-dplyr::distinct(dplyr::select(genes,ensembl_gene_id,external_gene_name))


# Supplementary table 1 - data sources ----------------------------------

st_1 <- data.frame(
  "Data source" = c(
    "INTERCEPT trial",
    "GTEx (version 8)",
    "BarcUVa-Seq",
    "eQTLGen (phase I)",
    "Colorectal cancer GWAS (overall)",
    "Colorectal cancer GWAS (by anatomical subsite)",
    "The Cancer Genome Atlas (TCGA) and GTEx, via TNMplot"
  ),
  "Data type" = c(
    "Gene expression microarray (Illumina HumanHT-12 v4 Expression BeadChip) before and after an eight-week low-energy liquid diet",
    "cis-eQTL summary statistics (European-ancestry subset)",
    "cis-eQTL summary statistics",
    "cis-eQTL summary statistics",
    "GWAS summary statistics for colorectal cancer",
    "GWAS summary statistics for colon, proximal colon, distal colon and rectal cancer",
    "RNA-Seq gene expression"
  ),
  "Tissue" = c(
    "Rectal mucosa (biopsies)",
    "Transverse and sigmoid colon, combined by fixed-effect meta-analysis",
    "Normal colon mucosa (ascending 31%, transverse 32%, descending 37%)",
    "Whole blood (80% of samples) and peripheral blood mononuclear cells (20%)",
    "Not applicable",
    "Not applicable",
    "Colon and rectal adenocarcinoma (TCGA); unpaired normal tissue from GTEx donors"
  ),
  "Population and ancestry" = c(
    "Adults with obesity (BMI \u2265 30 kg/m\u00B2), London, UK, recruited 2013\u20132014; full trial 90% White ethnicity",
    "Post-mortem donors, United States; European ancestry",
    "Individuals without colonic lesions at screening colonoscopy, Catalonia, Spain; almost all European ancestry",
    "37 datasets; predominantly European ancestry",
    "European ancestry; includes UK Biobank",
    "European ancestry; excludes UK Biobank",
    "Not restricted by ancestry"
  ),
  "Sample size" = c(
    "8 participants (of 20 in the full trial)",
    "Transverse colon n = 306; sigmoid colon n = 274",
    "n = 445",
    "Up to 31,684",
    "78,473 cases; 107,143 controls",
    "Colon 28,736; proximal colon 14,416; distal colon 12,879; rectum 14,150 cases; 43,099 controls",
    "Colon: 469 tumour, 315 normal; rectum: 166 tumour, 243 normal"
  ),
  "Adjustment" = c(
    "Age, sex, within-participant correlation and array quality (paired analysis)",
    "Top five genotype principal components, PEER factors, whole-genome sequencing platform, library construction protocol and sex",
    "Age, sex, sequencing batch, anatomical location, two genetic principal components and PEER factors",
    "Cohort-specific (see original publication)",
    "Varied by analytical unit (see original publication)",
    "Age, sex, study- or genotyping-project-specific covariates and genetic principal components",
    "None (unpaired comparison of TCGA tumour and GTEx normal tissue; Mann\u2013Whitney U test)"
  ),
  "Use in this study" = c(
    "Identification of genes differentially expressed following weight loss",
    "Genetic instruments for colon gene expression (MR); colocalisation",
    "Genetic instruments for colon gene expression (MR); colocalisation",
    "Genetic instruments for blood gene expression (MR); colocalisation",
    "MR outcome; colocalisation",
    "MR outcome; colocalisation",
    "Comparison of tumour and normal tissue expression for prioritised genes"
  ),
  "Access" = c(
    "Differential expression results in Supplementary table 2",
    "GTEx Portal (https://gtexportal.org)",
    "Colon Transcriptome Explorer version 2.0 (https://barcuvaseq.org/cotrex)",
    "eQTLGen Consortium (https://www.eqtlgen.org); IEU OpenGWAS (dataset IDs eqtl-a-[Ensembl gene ID])",
    "GWAS Catalog (GCST90255675)",
    "Genetics and Epidemiology of Colorectal Cancer Consortium (GECCO), on application",
    "TNMplot (https://tnmplot.com)"
  ),
  "Reference" = c(
    "Beeken et al., Obesity 2017 (ISRCTN35702367)",
    "GTEx Consortium, Science 2020; Barbeira et al., Genome Biol 2021",
    "D\u00EDez-Obrero et al., Cell Mol Gastroenterol Hepatol 2021",
    "V\u00F5sa et al., Nat Genet 2021",
    "Fernandez-Rozadilla et al., Nat Genet 2023",
    "Huyghe et al., Nat Genet 2019; Huyghe et al., Gut 2021",
    "Cancer Genome Atlas Network, Nature 2012; GTEx Consortium, Nat Genet 2013; Bartha and Gy\u0151rffy, Int J Mol Sci 2021"
  ),
  check.names = FALSE
)

# Supplementary table 2 - differential expression analysis ----------------

st_2<-fread("ABHD11 paper/Analysis/Results.csv")

#Annotate with Entrez gene ID
st_2_ref <- read.delim("Ensembl/GPL10558-50081.txt", quote="")

st_2<-merge(st_2,st_2_ref,by.x="Probe_Id",by.y="ID")

# Entrez to Ensembl mapping (org.Hs.eg.db), cached for reproducibility; called via namespaces to avoid masking dplyr::select
bm_cache <- "Ensembl/entrez_to_ensembl_orgHsegdb.csv"
if (file.exists(bm_cache)) {
  bm <- fread(bm_cache)
} else {
  bm <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = as.character(unique(na.omit(st_2$Entrez_Gene_ID))),
                              columns = "ENSEMBL", keytype = "ENTREZID")
  bm <- as.data.table(bm)[!is.na(ENSEMBL)]
  setnames(bm, c("ENTREZID", "ENSEMBL"), c("entrezgene_id", "ensembl_gene_id"))
  bm[, entrezgene_id := as.integer(entrezgene_id)]
  fwrite(bm, bm_cache)
}

st_2<-merge(st_2,bm,by.x="Entrez_Gene_ID",by.y="entrezgene_id",all.x=T)

st_2$FC <- 2^st_2$logFC

st_2<-dplyr::select(st_2,Symbol.x,Synonyms,ensembl_gene_id,Probe_Id,FC,logFC,`P.Value`)

colnames(st_2)<-c("Gene","Synonyms","Ensembl gene ID","Illumina microarray probe ID","Fold change","log2 fold change","P-value")

st_2$`Fold change`<-round(st_2$`Fold change`,digits=2)
st_2$`log2 fold change`<-round(st_2$`log2 fold change`,digits=3)
st_2$`P-value`<-signif(st_2$`P-value`,digits=3)

st_2 <- st_2 %>%
  mutate(across(where(is.character), ~ na_if(.x, "")))

st_2<-st_2[order(st_2$Gene),]

# Supplementary table 3 - F stats and r2 ----------------------------------

st_3<-fread("MR_Results/instrument_strength_summary.csv")

st_3<-merge(st_3,genes,by.x = "gene",by.y = "ensembl_gene_id",all.x=T)

st_3<-dplyr::select(st_3,external_gene_name,gene,dataset,n_snps,r2,F_statistic)

colnames(st_3)<-c("Gene","Ensembl gene ID","eQTL Dataset","Number of SNPs","r2","F-statistic")

st_3$r2<-round(st_3$r2,digits = 2)
st_3$`F-statistic`<-round(st_3$`F-statistic`,digits = 0)

st_3$`eQTL Dataset`<-gsub("GTExColon","GTEx",gsub("BARC","BarcUVa-Seq",st_3$`eQTL Dataset`))

# Supplementary table 4 - heterogeneity between GTEx colon tissues --------

st_4<-fread("GTEx/GTEx_instrument_heterogeneity.csv")

# GTEx instruments used in the MR
mr_instruments <- unique(fread("MR_Results/supplementary_harmonised_data.csv")[Dataset == "GTExColon", .(Gene, SNP)])
st_4<-merge(st_4,mr_instruments,by.x = c("exposure","SNP"),by.y = c("Gene","SNP"))

st_4<-merge(st_4,genes,by.x = "exposure",by.y = "ensembl_gene_id",all.x=T)

st_4<-dplyr::select(st_4,external_gene_name,exposure,SNP,variant_id,beta,se,pval,I2)

colnames(st_4)<-c("Gene","Ensembl gene ID","SNP","GTEx variant ID","Beta (meta-analysis)","SE (meta-analysis)","P-value (meta-analysis)","I2 (%)")

st_4$`Beta (meta-analysis)`<-round(st_4$`Beta (meta-analysis)`,digits=2)
st_4$`SE (meta-analysis)`<-round(st_4$`SE (meta-analysis)`,digits=2)
st_4$`P-value (meta-analysis)`<-signif(st_4$`P-value (meta-analysis)`,digits=3)
st_4$`I2 (%)`<-round(st_4$`I2 (%)`,digits=1)

st_4<-st_4[order(st_4$Gene,st_4$SNP),]


# Supplementary table 5 - harmonised summary genetic data -----------------

st_5<-fread("MR_Results/supplementary_harmonised_data.csv")

st_5<-merge(st_5,genes,by.x = "Gene",by.y = "ensembl_gene_id",all.x=T)

st_5<-dplyr::select(st_5,external_gene_name,Gene,Dataset,CRC_Subtype,SNP,Effect_Allele,Other_Allele,Beta_Exposure,SE_Exposure,Beta_Outcome,SE_Outcome)

colnames(st_5)<-c("Gene","Ensembl gene ID","eQTL Dataset","CRC Subtype","SNP", "Effect_Allele","Other_Allele","Beta_Exposure","SE_Exposure","Beta_Outcome","SE_Outcome")

st_5$`eQTL Dataset`<-gsub("GTExColon","GTEx",gsub("BARC","BarcUVa-Seq",st_5$`eQTL Dataset`))

st_5$`CRC Subtype`<-stringr::str_to_title(st_5$`CRC Subtype`)

st_5$Beta_Exposure<-round(st_5$Beta_Exposure,digits=2)
st_5$SE_Exposure<-round(st_5$SE_Exposure,digits=2)
st_5$Beta_Outcome<-round(st_5$Beta_Outcome,digits=2)
st_5$SE_Outcome<-round(st_5$SE_Outcome,digits=2)


# Supplementary table 6 - MR results --------------------------------------

st_6<-fread("MR_Results/MR_combined_results.csv")
st_6$OR<-exp(st_6$b)

st_6$LCI<-exp(st_6$b-1.96*st_6$se)
st_6$UCI<-exp(st_6$b+1.96*st_6$se)

st_6<-dplyr::select(st_6,external_gene_name,ensembl_gene_id,type,id.outcome,method,nsnp,OR,LCI,UCI,pval,bh_p)
st_6<-dplyr::distinct(st_6)

colnames(st_6)<-c("Gene","Ensembl gene ID","eQTL Dataset","CRC Subtype","Method", "Number of SNPs","Odds ratio","Lower 95% confidence interval","Upper 95% confidence interval","P-value","FDR-corrected P-value")

st_6$`eQTL Dataset`<-gsub("GTExColon","GTEx",gsub("BARC","BarcUVa-Seq",st_6$`eQTL Dataset`))

st_6$`CRC Subtype`<-stringr::str_to_title(st_6$`CRC Subtype`)

st_6$`Odds ratio`<-round(st_6$`Odds ratio`,digits=2)
st_6$`Lower 95% confidence interval`<-round(st_6$`Lower 95% confidence interval`,digits=2)
st_6$`Upper 95% confidence interval`<-round(st_6$`Upper 95% confidence interval`,digits=2)
st_6$`P-value`<-signif(st_6$`P-value`,digits=3)
st_6$`FDR-corrected P-value`<-signif(st_6$`FDR-corrected P-value`,digits=3)


# Supplementary table 7 - genetic colocalisation results ------------------

st_7<-fread("Coloc/pwcoco_combined_results.csv")

tissue_outcome <- sub("^ENSG[0-9]+_", "", st_7$analysis)
st_7$outcome <- sub("_.*$", "", tissue_outcome)
st_7$type <- sub("_coloc\\.coloc$", "", sub("^[^_]+_", "", tissue_outcome))

st_7<-merge(st_7,genes,by.x="gene",by.y="ensembl_gene_id")

st_7<-dplyr::select(st_7,external_gene_name,gene,type,outcome,SNP1,SNP2,nsnps,H0,H1,H2,H3,H4)

colnames(st_7)<-c("Gene","Ensembl gene ID","eQTL Dataset","CRC Subtype","SNP1","SNP2","Number of SNPs","H0","H1","H2","H3","H4")

st_7$`eQTL Dataset`<-gsub("GTExColon","GTEx",gsub("BARC","BarcUVa-Seq",st_7$`eQTL Dataset`))

st_7$`CRC Subtype`<-stringr::str_to_title(st_7$`CRC Subtype`)

st_7$H0<-round(st_7$H0,digits=2)
st_7$H1<-round(st_7$H1,digits=2)
st_7$H2<-round(st_7$H2,digits=2)
st_7$H3<-round(st_7$H3,digits=2)
st_7$H4<-round(st_7$H4,digits=2)

st_7<-st_7[order(st_7$Gene,st_7$`eQTL Dataset`,st_7$`CRC Subtype`,-st_7$H4),]


# Supplementary table 8 - TCGA analysis -----------------------------------

st_8 <- data.frame(
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

colnames(st_8)<-gsub("_"," ",colnames(st_8))



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
  "Supplementary table 7" = st_7,
  "Supplementary table 8" = st_8
)

for (sheet_name in names(tables)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, tables[[sheet_name]])
}

saveWorkbook(wb, file = "ABHD11 paper/Supplementary_tables.xlsx", overwrite = TRUE)
