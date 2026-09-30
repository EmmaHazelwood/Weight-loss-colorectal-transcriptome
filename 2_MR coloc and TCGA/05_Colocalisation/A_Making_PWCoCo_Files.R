sink("PWCoCo_File_Creation.txt")

library(data.table)
library(ieugwasr)
library(dplyr)
library(TwoSampleMR)
library(tidyverse)

sample_sizes <- list(
  overall = list(cases = 78473, controls = 107143, total = 185616),
  colon = list(cases = 28736, controls = 43099, total = 71835),
  proximal = list(cases = 14416, controls = 43099, total = 57515),
  distal = list(cases = 12879, controls = 43099, total = 55978),
  rectal = list(cases = 14150, controls = 43099, total = 57249)
)

outcome_files <- list(
  overall = "data/Fernandez-Rozadilla/GCST90255675_buildGRCh37.tsv",
  colon = "data/gecco/annotated/colon_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
  proximal = "data/gecco/annotated/proximal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
  distal = "data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
  rectal = "data/gecco/annotated/rectal_CRC_GWAS_noUKBio_summary_stats_annotated.txt"
)

exposure_dirs <- list(
  GTExColon = "data/ABHD11_paper/GTEx",
  BARC = "data/ABHD11_paper/BARCUVa-Seq"
)

mr_results <- fread("data/ABHD11_paper/MR_Results/MR_combined_results.csv")
sig_results <- mr_results %>% filter(bh_p < 0.05)

# eQTLGen
eqtl_results <- sig_results %>% filter(type == "eQTLGen")
# Only create files for analyses without existing PWCoCo results
eqtl_results <- eqtl_results[!file.exists(paste0(eqtl_results$ensembl_gene_id, "_", eqtl_results$id.outcome, "_eQTLGen_coloc.coloc")), ]

if (nrow(eqtl_results) > 0) {
  for (i in 1:nrow(eqtl_results)) {
    print(i)
    row <- eqtl_results[i, ]
    gene_id <- row$ensembl_gene_id
    outcome_type <- row$id.outcome
    print(gene_id)
    
    tryCatch({
      clumped_file <- "data/ABHD11_paper/eQTLGen/eQTLGenGenes.csv"
      if (!file.exists(clumped_file)) next
      
      clumped_data <- fread(clumped_file)
      gene_clumped <- clumped_data %>% filter(ensembl_gene_id == gene_id)
      
      if (nrow(gene_clumped) == 0) next
      
      lead_snp <- gene_clumped %>% arrange(pval.exposure) %>% head(1)
      lead_chr <- lead_snp$chr.exposure
      lead_pos <- lead_snp$pos.exposure
      
      window_start <- lead_pos - 500000
      window_end <- lead_pos + 500000
      
      eqtl_id <- paste0("eqtl-a-", gene_id)
      variants_query <- paste0(lead_chr, ":", window_start, "-", window_end)
      snps <- associations(id = eqtl_id, variants = variants_query)
      
      if (!is.data.frame(snps) || nrow(snps) == 0) {
        message("No eQTLGen data returned for ", eqtl_id, " - check OpenGWAS authentication (OPENGWAS_JWT)")
        next
      }
      
      exposure_dat <- format_data(
        snps,
        type = "exposure",
        snp_col = "rsid",
        beta_col = "beta",
        se_col = "se",
        effect_allele_col = "ea",
        other_allele_col = "nea",
        eaf_col = "eaf",
        pval_col = "p",
        chr_col = "chr",
        pos_col = "position"
      )
      
      if (outcome_type == "overall") {
        outcome_dat <- read_outcome_data(
          snps = snps$rsid,
          filename = outcome_files[[outcome_type]],
          sep = "\t",
          snp_col = "variant_id",
          beta_col = "beta",
          se_col = "standard_error",
          effect_allele_col = "effect_allele",
          other_allele_col = "other_allele",
          pval_col = "p_value"
        )
      } else {
        outcome_dat <- read_outcome_data(
          snps = snps$rsid,
          filename = outcome_files[[outcome_type]],
          sep = " ",
          snp_col = "SNP",
          beta_col = "Effect",
          se_col = "StdErr",
          effect_allele_col = "Allele1",
          other_allele_col = "Allele2",
          eaf_col = "Freq1"
        )
      }
      
      if (nrow(outcome_dat) == 0) next
      
      dat <- harmonise_data(exposure_dat, outcome_dat, action = 3)
      dat <- dat[dat$mr_keep == TRUE, ]
      
      if (nrow(dat) == 0) next
      
      exposure_out <- dat %>%
        dplyr::select(
          SNP = SNP,
          A1 = effect_allele.exposure,
          A2 = other_allele.exposure,
          A1_freq = eaf.exposure,
          beta = beta.exposure,
          se = se.exposure,
          p = pval.exposure
        )
      
      exposure_out$n <- median(snps$n)
      
      if (outcome_type == "overall") {
        outcome_out <- dat %>%
          dplyr::select(
            SNP = SNP,
            A1 = effect_allele.outcome,
            A2 = other_allele.outcome,
            A1_freq = eaf.exposure,
            beta = beta.outcome,
            se = se.outcome,
            p = pval.outcome
          )
      } else {
        outcome_out <- dat %>%
          dplyr::select(
            SNP = SNP,
            A1 = effect_allele.outcome,
            A2 = other_allele.outcome,
            A1_freq = eaf.outcome,
            beta = beta.outcome,
            se = se.outcome,
            p = pval.outcome
          )
      }
      outcome_out$n <- sample_sizes[[outcome_type]]$total
      outcome_out$case <- sample_sizes[[outcome_type]]$cases
      
      fwrite(exposure_out, paste0(gene_id, "_", outcome_type, "_eQTLGen_exposure.csv"))
      fwrite(outcome_out, paste0(gene_id, "_", outcome_type, "_eQTLGen_outcome.csv"))
      
    }, error = function(e) {
      print(paste("ERROR:", e$message))
    })
  }
}

# BARC
barc_results <- sig_results %>% filter(type == "BARC")
# Only create files for analyses without existing PWCoCo results
barc_results <- barc_results[!file.exists(paste0(barc_results$ensembl_gene_id, "_", barc_results$id.outcome, "_BARC_coloc.coloc")), ]

if (nrow(barc_results) > 0) {
  barc_full <- fread("data/BARCUVa-Seq/barcuvaseq.eqtls.allpairs.sumstats.txt")
  barc_full$phenotype_id <- substr(barc_full$gene_id, 1, 15)
  
  for (i in 1:nrow(barc_results)) {
    print(i)
    row <- barc_results[i, ]
    gene_id <- row$ensembl_gene_id
    outcome_type <- row$id.outcome
    print(gene_id)
    
    tryCatch({
      clumped_file <- file.path(exposure_dirs$BARC, paste0(gene_id, "_BARC_eqtl_clumped.csv"))
      if (!file.exists(clumped_file)) next
      
      clumped_data <- fread(clumped_file)
      
      if (nrow(clumped_data) == 0) next
      
      lead_snp <- clumped_data %>% arrange(pval.exposure) %>% head(1)
      lead_variant_id <- lead_snp$SNP
      
      idx <- which(barc_full$phenotype_id == gene_id)
      barc_window <- barc_full[idx, ]
      
      chr <- barc_window$chr[barc_window$variant_id == lead_variant_id]
      pos <- barc_window$pos[barc_window$variant_id == lead_variant_id]
      
      window_start <- pos - 500000
      window_end <- pos + 500000
      
      barc_window <- barc_window %>%
        filter(chr %in% chr, 
               pos >= window_start, 
               pos <= window_end)
      
      barc_window <- data.frame(barc_window)
      
      # BarcUVa-Seq reports only the minor allele frequency, so A1 frequencies come from the colorectal cancer GWAS
      exposure_dat <- format_data(
        barc_window,
        type = "exposure",
        header = TRUE,
        phenotype_col = "phenotype_id",
        snp_col = "variant_id",
        beta_col = "slope",
        se_col = "slope_se",
        effect_allele_col = "ALT",
        other_allele_col = "REF",
        pval_col = "pval_nominal",
        chr_col = "chr",
        pos_col = "pos"
      )
      
      if (outcome_type == "overall") {
        # The overall GWAS has no allele frequencies, so take them from the colon GWAS
        colon_dat <- read_outcome_data(
          snps = barc_window$variant_id,
          filename = outcome_files[["colon"]],
          sep = " ",
          snp_col = "SNP",
          beta_col = "Effect",
          se_col = "StdErr",
          effect_allele_col = "Allele1",
          other_allele_col = "Allele2",
          eaf_col = "Freq1"
        )
        
        if (nrow(colon_dat) == 0) next
        
        dat_colon <- harmonise_data(exposure_dat, colon_dat, action = 3)
        dat_colon <- dat_colon[dat_colon$mr_keep == TRUE, ]
        
        if (nrow(dat_colon) == 0) next
        
        outcome_dat <- read_outcome_data(
          snps = barc_window$variant_id,
          filename = outcome_files[[outcome_type]],
          sep = "\t",
          snp_col = "variant_id",
          beta_col = "beta",
          se_col = "standard_error",
          effect_allele_col = "effect_allele",
          other_allele_col = "other_allele",
          pval_col = "p_value"
        )
        
        if (nrow(outcome_dat) == 0) next
        
        dat <- harmonise_data(exposure_dat, outcome_dat, action = 3)
        dat <- dat[dat$mr_keep == TRUE, ]
        
        if (nrow(dat) == 0) next
        
        dat <- dat %>%
          left_join(
            dat_colon %>% dplyr::select(SNP, eaf.outcome, effect_allele.exposure, other_allele.exposure),
            by = "SNP",
            suffix = c("", ".colon")
          ) %>%
          mutate(
            A1_freq = case_when(
              effect_allele.exposure == effect_allele.exposure.colon ~ eaf.outcome.colon,
              effect_allele.exposure == other_allele.exposure.colon ~ 1 - eaf.outcome.colon,
              TRUE ~ NA_real_
            )
          )
      } else {
        outcome_dat <- read_outcome_data(
          snps = barc_window$variant_id,
          filename = outcome_files[[outcome_type]],
          sep = " ",
          snp_col = "SNP",
          beta_col = "Effect",
          se_col = "StdErr",
          effect_allele_col = "Allele1",
          other_allele_col = "Allele2",
          eaf_col = "Freq1"
        )
        
        if (nrow(outcome_dat) == 0) next
        
        dat <- harmonise_data(exposure_dat, outcome_dat, action = 3)
        dat <- dat[dat$mr_keep == TRUE, ]
        
        if (nrow(dat) == 0) next
        
        dat$A1_freq <- dat$eaf.outcome
      }
      
      exposure_out <- dat %>%
        dplyr::select(
          SNP = SNP,
          A1 = effect_allele.exposure,
          A2 = other_allele.exposure,
          A1_freq = A1_freq,
          beta = beta.exposure,
          se = se.exposure,
          p = pval.exposure
        )
      
      exposure_out$n <- 445  # BarcUVa-Seq individuals in the eQTL analysis
      
      outcome_out <- dat %>%
        dplyr::select(
          SNP = SNP,
          A1 = effect_allele.outcome,
          A2 = other_allele.outcome,
          A1_freq = A1_freq,
          beta = beta.outcome,
          se = se.outcome,
          p = pval.outcome
        )
      
      outcome_out$n <- sample_sizes[[outcome_type]]$total
      outcome_out$case <- sample_sizes[[outcome_type]]$cases
      
      fwrite(exposure_out, paste0(gene_id, "_", outcome_type, "_BARC_exposure.csv"))
      fwrite(outcome_out, paste0(gene_id, "_", outcome_type, "_BARC_outcome.csv"))
      
    }, error = function(e) {
      print(paste("ERROR:", e$message))
    })
  }
}

# GTEx Colon
gtex_results <- sig_results %>% filter(type == "GTExColon")
# Only create files for analyses without existing PWCoCo results
gtex_results <- gtex_results[!file.exists(paste0(gtex_results$ensembl_gene_id, "_", gtex_results$id.outcome, "_GTExColon_coloc.coloc")), ]

if (nrow(gtex_results) > 0) {
  gtex_full <- fread("data/GTEx-Colon/GTExMetaAnalysisResultsFixed.csv")
  snp_info <- fread("data/GTEx-Colon/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt", select = c("variant_id", "rs_id_dbSNP151_GRCh38p7", "chr", "variant_pos", "ref", "alt"))
  
  gtex_full$phenotype_id<-substr(gtex_full$phenotype_id,1,15)
  gtex_full<-gtex_full[gtex_full$phenotype_id %in% gtex_results$ensembl_gene_id,]
  
  snp_info<-dplyr::select(snp_info,variant_id,rs_id_dbSNP151_GRCh38p7,chr,variant_pos,ref,alt)
  gtex_full<-merge(gtex_full,snp_info)
  
  gtex_full$chr<-gsub("chr","",gtex_full$chr)
  
  gtex_full<-dplyr::select(gtex_full,phenotype_id,rs_id_dbSNP151_GRCh38p7,chr,variant_pos,alt,ref,beta,se,pval)
  colnames(gtex_full)<-c("exposure","SNP","chromosome.exposure","position.exposure","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","pval.exposure")
  
  for (i in 1:nrow(gtex_results)) {
    print(i)
    row <- gtex_results[i, ]
    gene_id <- row$ensembl_gene_id
    outcome_type <- row$id.outcome
    print(gene_id)
    
    tryCatch({
      clumped_file <- file.path(exposure_dirs$GTExColon, paste0(gene_id, "_colontotal_eqtl_clumped.csv"))
      if (!file.exists(clumped_file)) next
      
      clumped_data <- fread(clumped_file)
      
      if (nrow(clumped_data) == 0) next
      
      lead_snp <- clumped_data %>% arrange(pval.exposure) %>% head(1)
      lead_variant_id <- lead_snp$SNP
      
      gtex_window<-gtex_full[gtex_full$exposure==gene_id,]
      
      if (nrow(gtex_window) == 0) next
      
      chr <- gtex_full$chromosome.exposure[gtex_full$SNP == lead_variant_id]
      pos <- gtex_full$position.exposure[gtex_full$SNP == lead_variant_id]
      
      window_start <- pos - 500000
      window_end <- pos + 500000
      
      gtex_window <- gtex_window %>%
        filter(chr %in% chromosome.exposure, 
               position.exposure >= window_start, 
               position.exposure <= window_end)
      
      gtex_window <- data.frame(gtex_window)
      
      if (nrow(gtex_window) == 0) next
      
      exposure_dat <- format_data(
        gtex_window,
        type = "exposure",
        snp_col = "SNP",
        beta_col = "beta.exposure",
        se_col = "se.exposure",
        effect_allele_col = "effect_allele.exposure",
        other_allele_col = "other_allele.exposure",
        pval_col = "pval.exposure",
        phenotype_col = "exposure",
        chr_col = "chromosome.exposure",
        pos_col = "position.exposure"
      )
      
      if (outcome_type == "overall") {
        # The overall GWAS has no allele frequencies, so take them from the colon GWAS
        colon_dat <- read_outcome_data(
          snps = gtex_window$rs_id_dbSNP151_GRCh38p7,
          filename = outcome_files[["colon"]],
          sep = " ",
          snp_col = "SNP",
          beta_col = "Effect",
          se_col = "StdErr",
          effect_allele_col = "Allele1",
          other_allele_col = "Allele2",
          eaf_col = "Freq1"
        )
        
        if (nrow(colon_dat) == 0) next
        
        dat_colon <- harmonise_data(exposure_dat, colon_dat, action = 3)
        dat_colon <- dat_colon[dat_colon$mr_keep == TRUE, ]
        
        if (nrow(dat_colon) == 0) next
        
        outcome_dat <- read_outcome_data(
          snps = gtex_window$rs_id_dbSNP151_GRCh38p7,
          filename = outcome_files[[outcome_type]],
          sep = "\t",
          snp_col = "variant_id",
          beta_col = "beta",
          se_col = "standard_error",
          effect_allele_col = "effect_allele",
          other_allele_col = "other_allele",
          pval_col = "p_value"
        )
        
        if (nrow(outcome_dat) == 0) next
        
        dat <- harmonise_data(exposure_dat, outcome_dat, action = 3)
        dat <- dat[dat$mr_keep == TRUE, ]
        
        if (nrow(dat) == 0) next
        
        dat <- dat %>%
          left_join(
            dat_colon %>% dplyr::select(SNP, eaf.outcome, effect_allele.exposure, other_allele.exposure),
            by = "SNP",
            suffix = c("", ".colon")
          ) %>%
          mutate(
            eaf_from_colon = case_when(
              effect_allele.exposure == effect_allele.exposure.colon ~ eaf.outcome.colon,
              effect_allele.exposure == other_allele.exposure.colon ~ 1 - eaf.outcome.colon,
              TRUE ~ NA_real_
            )
          )
        
      } else {
        outcome_dat <- read_outcome_data(
          snps = gtex_window$rs_id_dbSNP151_GRCh38p7,
          filename = outcome_files[[outcome_type]],
          sep = " ",
          snp_col = "SNP",
          beta_col = "Effect",
          se_col = "StdErr",
          effect_allele_col = "Allele1",
          other_allele_col = "Allele2",
          eaf_col = "Freq1"
        )
        
        if (nrow(outcome_dat) == 0) next
        
        dat <- harmonise_data(exposure_dat, outcome_dat, action = 3)
        dat <- dat[dat$mr_keep == TRUE, ]
        
        if (nrow(dat) == 0) next
      }
      
      if (outcome_type == "overall") {
        exposure_out <- dat %>%
          dplyr::select(
            SNP = SNP,
            A1 = effect_allele.exposure,
            A2 = other_allele.exposure,
            A1_freq = eaf_from_colon,
            beta = beta.exposure,
            se = se.exposure,
            p = pval.exposure
          )
      } else {
        exposure_out <- dat %>%
          dplyr::select(
            SNP = SNP,
            A1 = effect_allele.exposure,
            A2 = other_allele.exposure,
            A1_freq = eaf.outcome,
            beta = beta.exposure,
            se = se.exposure,
            p = pval.exposure
          )
      }
      
      exposure_out$n <- 580  # European-ancestry GTEx v8: transverse 306 + sigmoid 274
      
      if (outcome_type == "overall") {
        outcome_out <- dat %>%
          dplyr::select(
            SNP = SNP,
            A1 = effect_allele.outcome,
            A2 = other_allele.outcome,
            A1_freq = eaf_from_colon,
            beta = beta.outcome,
            se = se.outcome,
            p = pval.outcome
          )
      } else {
        outcome_out <- dat %>%
          dplyr::select(
            SNP = SNP,
            A1 = effect_allele.outcome,
            A2 = other_allele.outcome,
            A1_freq = eaf.outcome,
            beta = beta.outcome,
            se = se.outcome,
            p = pval.outcome
          )
      }
      outcome_out$n <- sample_sizes[[outcome_type]]$total
      outcome_out$case <- sample_sizes[[outcome_type]]$cases
      
      fwrite(exposure_out, paste0(gene_id, "_", outcome_type, "_GTExColon_exposure.csv"))
      fwrite(outcome_out, paste0(gene_id, "_", outcome_type, "_GTExColon_outcome.csv"))
      
    }, error = function(e) {
      print(paste("ERROR:", e$message))
    })
  }
}

sink()