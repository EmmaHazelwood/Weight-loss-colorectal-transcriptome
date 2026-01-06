library(dplyr)
library(data.table)
library(ieugwasr)
library(TwoSampleMR)
library(MRInstruments)
library(LDlinkR)
library(tidyr)

# Initialize data frames for collecting results across all analyses
all_fstats <- data.frame()
all_harmonised <- data.frame()

# Function to check if a SNP is palindromic (A/T or G/C)
is_palindromic <- function(a1, a2) {
  a1 <- toupper(as.character(a1))
  a2 <- toupper(as.character(a2))
  (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | 
    (a1 == "G" & a2 == "C") | (a1 == "C" & a2 == "G")
}

# Function to perform proxy SNP lookup
find_proxy_snps <- function(snps_for_proxy, exposure_dat, outcome_snps, token = "a563a6f57fcc") {
  proxySNPs <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(proxySNPs) <- c("missing_snp", "RS_Number", "REF", "ALT")
  
  for(j in snps_for_proxy) {
    print(paste("Looking for proxy for:", j))
    
    # Try to get proxy, handle errors (e.g., non-biallelic variants)
    proxy <- tryCatch({
      LDproxy(j, pop = "EUR", r2d = "r2", token = token)
    }, error = function(e) {
      print(paste("  Error getting proxy for", j, ":", e$message))
      return(NULL)
    })
    
    # Check if proxy query was successful
    if (is.null(proxy) || nrow(proxy) == 0) {
      print(paste("  No proxy data available for:", j))
      next
    }
    
    # Check for error messages in returned data
    if (ncol(proxy) == 1 && any(grepl("error", proxy[,1], ignore.case = TRUE))) {
      print(paste("  LDproxy error for", j, ":", proxy[1,1]))
      next
    }
    
    proxy2 <- proxy %>% 
      dplyr::filter(R2 > 0.8) %>%
      dplyr::filter(RS_Number != ".") %>%
      dplyr::filter(RS_Number %in% outcome_snps)
    
    if (nrow(proxy2) > 0) {
      # Try each proxy in order until we find a non-palindromic one
      proxy_found <- FALSE
      
      for (proxy_idx in 1:nrow(proxy2)) {
        proxy3 <- proxy2[proxy_idx, ]
        proxy3$missing_snp <- j
        
        # Extract effect and other alleles by column name, not position
        ea <- exposure_dat$effect_allele.exposure[which(exposure_dat$SNP == j)]
        oa <- exposure_dat$other_allele.exposure[which(exposure_dat$SNP == j)]
        
        # Take first value if multiple rows (shouldn't happen but just in case)
        ea <- ea[1]
        oa <- oa[1]
        
        Correlated_Alleles <- data.frame(proxy3$Correlated_Alleles)
        Correlated_Alleles$proxy3.Correlated_Alleles <- as.character(Correlated_Alleles$proxy3.Correlated_Alleles)
        Correlated_Alleles <- Correlated_Alleles %>% 
          separate(proxy3.Correlated_Alleles, c("a1exp", "a1out", "a2exp", "a2out"))
        Correlated_Alleles <- na.omit(Correlated_Alleles)
        
        if (nrow(Correlated_Alleles) != 0) {
          proxy4 <- cbind(proxy3, Correlated_Alleles)
          proxy4$ea <- ea
          proxy4$oa <- oa
          
          if (proxy4$a1exp == proxy4$ea) {
            proxy4$REF <- proxy4$a1out
            proxy4$ALT <- proxy4$a2out
          } else if(proxy4$a2exp == proxy4$ea) {
            proxy4$REF <- proxy4$a2out
            proxy4$ALT <- proxy4$a1out
          } else {
            proxy4$REF <- "NA"
            proxy4$ALT <- "NA"
          }
          
          # Check if the proxy SNP is palindromic
          if (proxy4$REF != "NA" && proxy4$ALT != "NA") {
            if (!is_palindromic(proxy4$REF, proxy4$ALT)) {
              # Found a non-palindromic proxy!
              proxy5 <- select(proxy4, missing_snp, RS_Number, REF, ALT)
              proxySNPs <- rbind(proxySNPs, proxy5)
              print(paste("  Found non-palindromic proxy:", proxy4$RS_Number, 
                          paste0(proxy4$REF, "/", proxy4$ALT)))
              proxy_found <- TRUE
              break  # Stop looking for this SNP
            } else {
              print(paste("  Skipping palindromic proxy:", proxy4$RS_Number, 
                          paste0(proxy4$REF, "/", proxy4$ALT)))
            }
          }
        }
      }
      
      if (!proxy_found) {
        print(paste("  No non-palindromic proxy found for:", j))
      }
    } else {
      print(paste("  No suitable proxies found for:", j))
    }
  }
  return(proxySNPs)
}

# Function to merge proxy SNPs with exposure data
merge_proxy_snps <- function(exposure_dat, proxySNPs, snps_for_proxy) {
  exposure_datproxies <- exposure_dat[which(exposure_dat$SNP %in% snps_for_proxy), ]
  
  if (nrow(proxySNPs) != 0) {
    exposure_dat2 <- merge(exposure_datproxies, proxySNPs, by.x = "SNP", by.y = "missing_snp")
    exposure_dat2$SNP <- exposure_dat2$RS_Number
    exposure_dat2$effect_allele.exposure <- exposure_dat2$REF
    exposure_dat2$other_allele.exposure <- exposure_dat2$ALT
    exposure_dat2 <- exposure_dat2[, 1:12]
    exposure_dat <- rbind(exposure_dat, exposure_dat2)
  }
  return(exposure_dat)
}

# Function to calculate r² and F-statistics
calculate_r2_and_f <- function(dat, sample_size) {
  # Check if EAF data is available
  if (all(is.na(dat$eaf.outcome))) {
    print("  No EAF data available - skipping r² and F-statistic calculation")
    return(list(r2 = NA, F_stat = NA, n_snps = nrow(dat)))
  }
  
  # Calculate proportion of variance explained (PVE) for each SNP
  dat$num <- 2 * (dat$beta.exposure^2) * dat$eaf.outcome * (1 - dat$eaf.outcome)
  dat$den <- 2 * (dat$beta.exposure^2) * dat$eaf.outcome * (1 - dat$eaf.outcome) + 
    ((dat$se.exposure^2) * 2 * sample_size * dat$eaf.outcome * (1 - dat$eaf.outcome))
  dat$pve <- dat$num / dat$den
  
  # Sum PVE across all SNPs (total r²)
  pve_total <- sum(dat$pve, na.rm = TRUE)
  
  # Calculate F-statistic
  F_stat <- (pve_total * (sample_size - 2)) / (1 - pve_total)
  
  return(list(r2 = pve_total, F_stat = F_stat, n_snps = nrow(dat)))
}

# Function to perform MR analysis with instrument strength calculation
perform_mr_analysis <- function(exposure_dat, outcome_dat, gene_id, eqtl_type, 
                                outcome_type, sample_size, outcome_n_cases = NA, outcome_n_controls = NA) {
  dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat, action = 3)
  dat <- dat[dat$mr_keep == TRUE, ]
  
  if (nrow(dat) == 0) {
    print("No rows after harmonization")
    return(NULL)
  }
  
  # Calculate instrument strength metrics
  print(paste("  Calculating instrument strength for", nrow(dat), "SNPs"))
  strength_metrics <- calculate_r2_and_f(dat, sample_size)
  
  # Only print if values are not NA
  if (!is.na(strength_metrics$r2) && !is.na(strength_metrics$F_stat)) {
    print(paste("  r² =", round(strength_metrics$r2, 4), 
                "| F-statistic =", round(strength_metrics$F_stat, 2)))
  }
  
  # Perform MR
  if (nrow(dat) == 1) {
    res <- mr(dat, method_list = "mr_wald_ratio")
  } else {
    res <- mr(dat, method_list = "mr_ivw")
  }
  
  # Add metadata and instrument strength metrics
  res$id.exposure <- gene_id
  res$type <- eqtl_type
  res$id.outcome <- outcome_type
  res$n_snps_used <- strength_metrics$n_snps
  res$r2 <- strength_metrics$r2
  res$F_statistic <- strength_metrics$F_stat
  res$sample_size <- sample_size
  
  # Store F statistics and r² in global collection (only once per gene-dataset pair)
  # Only store if we have valid r² and F values
  if (!is.na(strength_metrics$r2) && !is.na(strength_metrics$F_stat)) {
    existing <- all_fstats[all_fstats$gene == gene_id & all_fstats$dataset == eqtl_type, ]
    
    if (nrow(existing) == 0) {
      fstat_row <- data.frame(
        gene = gene_id,
        dataset = eqtl_type,
        n_snps = strength_metrics$n_snps,
        r2 = strength_metrics$r2,
        F_statistic = strength_metrics$F_stat
      )
      all_fstats <<- rbind(all_fstats, fstat_row)
    }
  }
  
  # Store harmonised data with metadata for supplementary table
  # Calculate total outcome sample size
  outcome_n_total <- NA
  if (!is.na(outcome_n_cases) && !is.na(outcome_n_controls)) {
    outcome_n_total <- outcome_n_cases + outcome_n_controls
  }
  
  harmonised_subset <- dat %>%
    select(SNP, effect_allele.exposure, other_allele.exposure, 
           beta.exposure, se.exposure, beta.outcome, se.outcome) %>%
    mutate(
      gene = gene_id,
      dataset = eqtl_type,
      outcome = outcome_type,
      N_exposure = sample_size,
      N_outcome = outcome_n_total
    ) %>%
    select(gene, dataset, outcome, SNP, effect_allele.exposure, 
           other_allele.exposure, beta.exposure, se.exposure, N_exposure,
           beta.outcome, se.outcome, N_outcome)
  
  all_harmonised <<- rbind(all_harmonised, harmonised_subset)
  
  return(res)
}

# Function to run GTEx/BARC analysis
run_tissue_eqtl_analysis <- function(pattern, tissue_name, outcome_file, outcome_snps, 
                                     outcome_type, ld_token, sample_size, 
                                     outcome_n_cases = NA, outcome_n_controls = NA,
                                     snp_col = "SNP", data_dir = ".") {
  print(paste("Processing:", tissue_name))
  
  # Get only clumped files matching the pattern in the specified directory
  x <- dir(path = data_dir, pattern = pattern)
  x <- substr(x, 1, 15)
  x <- unique(x)
  
  for (i in 1:length(x)) {
    # Only get files matching both the gene ID AND the clumped pattern
    y <- dir(path = data_dir, pattern = paste0(x[i], ".*", pattern))
    print(paste(i, x[i]))
    
    # Double check we have the clumped file
    clumped_file <- grep(pattern, y, value = TRUE)
    if (length(clumped_file) == 0) {
      print(paste("No clumped file found for", x[i]))
      next
    }
    
    # Construct full path to file
    full_path <- file.path(data_dir, clumped_file[1])
    
    exposure_dat <- read_exposure_data(
      filename = full_path,
      snp_col = snp_col,
      sep = ",",
      beta_col = "beta.exposure",
      se_col = "se.exposure",
      effect_allele_col = "effect_allele.exposure",
      other_allele_col = "other_allele.exposure",
      eaf_col = "eaf.exposure"
    )
    
    # Identify palindromic SNPs
    palindromic_snps <- exposure_dat$SNP[is_palindromic(exposure_dat$effect_allele.exposure, 
                                                        exposure_dat$other_allele.exposure)]
    
    if (length(palindromic_snps) > 0) {
      print(paste("  Found", length(palindromic_snps), "palindromic SNPs, searching for non-palindromic proxies"))
    }
    
    # Find SNPs that need proxies: either missing from outcome OR palindromic
    snps_missing_from_outcome <- setdiff(exposure_dat$SNP, outcome_snps)
    snps_needing_proxy <- unique(c(snps_missing_from_outcome, palindromic_snps))
    
    if (length(snps_needing_proxy) != 0) {
      print(paste("  Total SNPs needing proxies:", length(snps_needing_proxy), 
                  "(", length(snps_missing_from_outcome), "missing from outcome,",
                  length(palindromic_snps), "palindromic )"))
      proxySNPs <- find_proxy_snps(snps_needing_proxy, exposure_dat, outcome_snps, ld_token)
      
      if (nrow(proxySNPs) > 0) {
        # Remove original SNPs that we found proxies for
        exposure_dat <- exposure_dat[!exposure_dat$SNP %in% proxySNPs$missing_snp, ]
        
        # Merge in the proxy SNPs
        exposure_dat <- merge_proxy_snps(exposure_dat, proxySNPs, snps_needing_proxy)
        
        print(paste("  Replaced", nrow(proxySNPs), "SNPs with non-palindromic proxies"))
      }
    }
    
    if (length(intersect(exposure_dat$SNP, outcome_snps)) != 0) {
      outcome_dat <- read_outcome_data(
        snps = exposure_dat$SNP,
        filename = outcome_file,
        sep = ",",
        snp_col = "SNP",
        beta_col = "Effect",
        se_col = "StdErr",
        effect_allele_col = "Allele1",
        other_allele_col = "Allele2",
        eaf_col = "Freq1"
      )
      
      res <- perform_mr_analysis(exposure_dat, outcome_dat, x[i], tissue_name, 
                                 outcome_type, sample_size, outcome_n_cases, outcome_n_controls)
      
      if (!is.null(res)) {
        output_filename <- paste(x[i], tissue_name, outcome_type, "MR_Results.csv", sep = "_")
        fwrite(res, output_filename, quote = FALSE)
      }
    } else {
      print("No SNPs in common")
    }
  }
}

# Function to run eQTLGen analysis
run_eqtlgen_analysis <- function(outcome_file, outcome_snps, outcome_type, ld_token, 
                                 outcome_n_cases = NA, outcome_n_controls = NA, data_dir = ".") {
  print("Processing: eQTLGen")
  
  # Load pre-processed and filtered eQTLGen data
  eqtlgen_file <- file.path(data_dir, "eQTLGenGenes.csv")
  exposure_dat <- fread(eqtlgen_file)
  
  # Ensure gene column exists
  if (!"gene" %in% colnames(exposure_dat)) {
    exposure_dat$gene <- exposure_dat$ensembl_gene_id
  }
  
  # Identify palindromic SNPs
  palindromic_snps <- exposure_dat$SNP[is_palindromic(exposure_dat$effect_allele.exposure, 
                                                      exposure_dat$other_allele.exposure)]
  
  if (length(palindromic_snps) > 0) {
    print(paste("  Found", length(palindromic_snps), "palindromic SNPs, searching for non-palindromic proxies"))
  }
  
  # Find SNPs that need proxies
  snps_missing_from_outcome <- setdiff(exposure_dat$SNP, outcome_snps)
  snps_needing_proxy <- unique(c(snps_missing_from_outcome, palindromic_snps))
  
  if (length(snps_needing_proxy) != 0) {
    print(paste("  Total SNPs needing proxies:", length(snps_needing_proxy), 
                "(", length(snps_missing_from_outcome), "missing from outcome,",
                length(palindromic_snps), "palindromic )"))
    proxySNPs <- data.frame(matrix(ncol = 4, nrow = 0))
    
    for(j in snps_needing_proxy) {
      print(paste("Looking for proxy for:", j))
      
      proxy <- tryCatch({
        LDproxy(j, pop = "EUR", r2d = "r2", token = ld_token)
      }, error = function(e) {
        print(paste("  Error getting proxy for", j, ":", e$message))
        return(NULL)
      })
      
      if (is.null(proxy) || nrow(proxy) == 0) {
        print(paste("  No proxy data available for:", j))
        next
      }
      
      if (ncol(proxy) == 1 && any(grepl("error", proxy[,1], ignore.case = TRUE))) {
        print(paste("  LDproxy error for", j, ":", proxy[1,1]))
        next
      }
      
      proxy2 <- proxy %>% 
        dplyr::filter(R2 > 0.8) %>%
        dplyr::filter(RS_Number != ".") %>%
        dplyr::filter(RS_Number %in% outcome_snps)
      
      if (nrow(proxy2) > 0) {
        proxy_found <- FALSE
        
        for (proxy_idx in 1:nrow(proxy2)) {
          proxy3 <- proxy2[proxy_idx, ]
          proxy3$missing_snp <- j
          
          ea <- exposure_dat$effect_allele.exposure[which(exposure_dat$SNP == j)]
          oa <- exposure_dat$other_allele.exposure[which(exposure_dat$SNP == j)]
          ea <- unique(ea)[1]
          oa <- unique(oa)[1]
          
          Correlated_Alleles <- data.frame(proxy3$Correlated_Alleles)
          Correlated_Alleles$proxy3.Correlated_Alleles <- as.character(Correlated_Alleles$proxy3.Correlated_Alleles)
          Correlated_Alleles <- Correlated_Alleles %>% 
            separate(proxy3.Correlated_Alleles, c("a1exp", "a1out", "a2exp", "a2out"))
          Correlated_Alleles <- na.omit(Correlated_Alleles)
          
          if (nrow(Correlated_Alleles) != 0) {
            proxy4 <- cbind(proxy3, Correlated_Alleles)
            proxy4$ea <- ea
            proxy4$oa <- oa
            
            if (proxy4$a1exp == proxy4$ea) {
              proxy4$REF <- proxy4$a1out
              proxy4$ALT <- proxy4$a2out
            } else if(proxy4$a2exp == proxy4$ea) {
              proxy4$REF <- proxy4$a2out
              proxy4$ALT <- proxy4$a1out
            } else {
              proxy4$REF <- "NA"
              proxy4$ALT <- "NA"
            }
            
            if (proxy4$REF != "NA" && proxy4$ALT != "NA") {
              if (!is_palindromic(proxy4$REF, proxy4$ALT)) {
                proxy5 <- select(proxy4, missing_snp, RS_Number, REF, ALT)
                proxySNPs <- rbind(proxySNPs, proxy5)
                proxySNPs <- distinct(proxySNPs)
                print(paste("  Found non-palindromic proxy:", proxy4$RS_Number, 
                            paste0(proxy4$REF, "/", proxy4$ALT)))
                proxy_found <- TRUE
                break
              } else {
                print(paste("  Skipping palindromic proxy:", proxy4$RS_Number, 
                            paste0(proxy4$REF, "/", proxy4$ALT)))
              }
            }
          }
        }
        
        if (!proxy_found) {
          print(paste("  No non-palindromic proxy found for:", j))
        }
      } else {
        print(paste("  No suitable proxies found for:", j))
      }
    }
    
    if (nrow(proxySNPs) > 0) {
      exposure_datproxies <- exposure_dat[which(exposure_dat$SNP %in% snps_needing_proxy), ]
      
      if (nrow(proxySNPs) != 0) {
        exposure_dat2 <- merge(exposure_datproxies, proxySNPs, by.x = "SNP", by.y = "missing_snp")
        exposure_dat2$SNP <- exposure_dat2$RS_Number
        exposure_dat2$effect_allele.exposure <- exposure_dat2$REF
        exposure_dat2$other_allele.exposure <- exposure_dat2$ALT
        exposure_dat2 <- exposure_dat2[, 1:21]
        exposure_dat <- rbind(exposure_dat, exposure_dat2)
        
        print(paste("  Replaced", nrow(proxySNPs), "SNPs with non-palindromic proxies"))
      }
    }
  }
  
  outcome_dat <- read_outcome_data(
    snps = exposure_dat$SNP,
    filename = outcome_file,
    sep = ",",
    snp_col = "SNP",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1"
  )
  
  dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat, action = 3)
  dat <- dat[dat$mr_keep == TRUE, ]
  
  genes <- unique(dat$gene)
  
  for (i in 1:length(genes)) {
    print(paste(i, genes[i]))
    gene <- genes[i]
    dat2 <- dat[which(dat$gene == gene), ]
    
    if (nrow(dat2) == 0) {
      print("No rows")
      next
    }
    
    # Get sample size for this gene (use minimum if varies across SNPs)
    sample_size <- min(dat2$samplesize.exposure, na.rm = TRUE)
    
    # Calculate instrument strength
    print(paste("  Calculating instrument strength for", nrow(dat2), "SNPs"))
    strength_metrics <- calculate_r2_and_f(dat2, sample_size)
    
    # Only print if values are not NA
    if (!is.na(strength_metrics$r2) && !is.na(strength_metrics$F_stat)) {
      print(paste("  r² =", round(strength_metrics$r2, 4), 
                  "| F-statistic =", round(strength_metrics$F_stat, 2)))
    }
    
    # Perform MR
    if (nrow(dat2) == 1) {
      res <- mr(dat2, method_list = "mr_wald_ratio")
    } else {
      res <- mr(dat2, method_list = "mr_ivw")
    }
    
    # Add metadata and instrument strength metrics
    res$id.exposure <- genes[i]
    res$type <- "eQTLGen"
    res$id.outcome <- outcome_type
    res$n_snps_used <- strength_metrics$n_snps
    res$r2 <- strength_metrics$r2
    res$F_statistic <- strength_metrics$F_stat
    res$sample_size <- sample_size
    
    # Store F statistics (only once per gene-dataset pair and only if valid)
    if (!is.na(strength_metrics$r2) && !is.na(strength_metrics$F_stat)) {
      existing <- all_fstats[all_fstats$gene == genes[i] & all_fstats$dataset == "eQTLGen", ]
      
      if (nrow(existing) == 0) {
        fstat_row <- data.frame(
          gene = genes[i],
          dataset = "eQTLGen",
          n_snps = strength_metrics$n_snps,
          r2 = strength_metrics$r2,
          F_statistic = strength_metrics$F_stat
        )
        all_fstats <<- rbind(all_fstats, fstat_row)
      }
    }
    
    # Store harmonised data
    outcome_n_total <- NA
    if (!is.na(outcome_n_cases) && !is.na(outcome_n_controls)) {
      outcome_n_total <- outcome_n_cases + outcome_n_controls
    }
    
    harmonised_subset <- dat2 %>%
      select(SNP, effect_allele.exposure, other_allele.exposure, 
             beta.exposure, se.exposure, beta.outcome, se.outcome) %>%
      mutate(
        gene = genes[i],
        dataset = "eQTLGen",
        outcome = outcome_type,
        N_exposure = sample_size,
        N_outcome = outcome_n_total
      ) %>%
      select(gene, dataset, outcome, SNP, effect_allele.exposure, 
             other_allele.exposure, beta.exposure, se.exposure, N_exposure,
             beta.outcome, se.outcome, N_outcome)
    
    all_harmonised <<- rbind(all_harmonised, harmonised_subset)
    
    fwrite(res, paste(genes[i], "eQTLGen", outcome_type, "MR_Results.csv", sep = "_"), quote = FALSE)
  }
}


# Main analysis -----------------------------------------------------------

# Define outcome datasets
outcomes <- list(
  colon = list(
    file = "/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/colon_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
    name = "colon",
    temp_file = "/home/sw20203/MiniProject2/genexpression/Proxy/tempcolon.csv",
    n_cases = 28736,
    n_controls = 43099
  ),
  overall = list(
    file = "/projects/MRC_IEU/research/projects/icep2/wp2/049/working/data/Fernandez-Rozadilla/GCST90255675_buildGRCh37.tsv",
    name = "overall",
    snp_col = "variant_id",
    beta_col = "beta",
    se_col = "standard_error",
    ea_col = "effect_allele",
    oa_col = "other_allele",
    chr_col = "chromosome",
    pos_col = "base_pair_location",
    n_cases = 78473,
    n_controls = 107143
  ),
  distal = list(
    file = "/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/distal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
    name = "distal",
    temp_file = "/home/sw20203/MiniProject2/genexpression/Proxy/tempdistal.csv",
    n_cases = 12879,
    n_controls = 43099
  ),
  proximal = list(
    file = "/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/proximal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
    name = "proximal",
    temp_file = "/home/sw20203/MiniProject2/genexpression/Proxy/tempproximal.csv",
    n_cases = 14416,
    n_controls = 43099
  ),
  rectal = list(
    file = "/projects/MRC_IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/rectal_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
    name = "rectal",
    temp_file = "/home/sw20203/MiniProject2/genexpression/Proxy/temprectal.csv",
    n_cases = 14150,
    n_controls = 43099
  )
)

ld_token <- "a563a6f57fcc"

# Define sample sizes for each dataset
sample_sizes <- list(
  GTExColon = 686,
  BARC = 445
  # eQTLGen sample sizes vary by gene and are in the data file
)

# Process each outcome type
for (outcome_name in names(outcomes)) {
  outcome_info <- outcomes[[outcome_name]]
  
  sink(paste0("PerformingMR_", outcome_info$name, ".txt"))
  
  print(paste("Processing outcome:", outcome_info$name))
  
  # Load outcome data
  CRC <- fread(outcome_info$file)
  
  # For overall analysis, create temp file with standardized column names
  if (outcome_name == "overall") {
    temp_crc <- data.frame(
      SNP = CRC[[outcome_info$snp_col]],
      Effect = CRC[[outcome_info$beta_col]],
      StdErr = CRC[[outcome_info$se_col]],
      Allele1 = CRC[[outcome_info$ea_col]],
      Allele2 = CRC[[outcome_info$oa_col]],
      chr = CRC[[outcome_info$chr_col]],
      pos = CRC[[outcome_info$pos_col]]
    )
    
    temp_crc$Freq1 <- NA
    
    temp_file <- "/home/sw20203/MiniProject2/genexpression/Proxy/tempoverall.csv"
    fwrite(temp_crc, temp_file, quote = FALSE)
    outcome_file <- temp_file
    outcome_snps <- temp_crc$SNP
  } else {
    outcome_file <- outcome_info$temp_file
    outcome_snps <- CRC$SNP
  }
  
  # Define data directories
  gtex_dir <- "/projects/MRC_IEU/research/projects/icep2/wp2/034/working/data/ABHD11_paper/GTEx"
  barc_dir <- "/projects/MRC_IEU/research/projects/icep2/wp2/034/working/data/ABHD11_paper/BARCUVa-Seq"
  eqtlgen_dir <- "/projects/MRC_IEU/research/projects/icep2/wp2/034/working/data/ABHD11_paper/eQTLGen"
  
  # Run GTEx colon total eQTL analysis
  run_tissue_eqtl_analysis(
    pattern = "_colontotal_eqtl_clumped.csv",
    tissue_name = "GTExColon",
    outcome_file = outcome_file,
    outcome_snps = outcome_snps,
    outcome_type = outcome_info$name,
    ld_token = ld_token,
    sample_size = sample_sizes$GTExColon,
    outcome_n_cases = outcome_info$n_cases,
    outcome_n_controls = outcome_info$n_controls,
    snp_col = "SNP",
    data_dir = gtex_dir
  )
  
  # Run BARC eQTL analysis
  run_tissue_eqtl_analysis(
    pattern = "_BARC_eqtl_clumped.csv",
    tissue_name = "BARC",
    outcome_file = outcome_file,
    outcome_snps = outcome_snps,
    outcome_type = outcome_info$name,
    ld_token = ld_token,
    sample_size = sample_sizes$BARC,
    outcome_n_cases = outcome_info$n_cases,
    outcome_n_controls = outcome_info$n_controls,
    snp_col = "SNP",
    data_dir = barc_dir
  )
  
  # Run eQTLGen analysis (sample sizes are in the data file)
  run_eqtlgen_analysis(
    outcome_file = outcome_file,
    outcome_snps = outcome_snps,
    outcome_type = outcome_info$name,
    ld_token = ld_token,
    outcome_n_cases = outcome_info$n_cases,
    outcome_n_controls = outcome_info$n_controls,
    data_dir = eqtlgen_dir
  )
  
  sink()
  print(paste("Completed analysis for:", outcome_info$name))
}

print("All analyses complete!")
print("Results include r² and F-statistics calculated from harmonized SNPs used in each MR analysis.")

# Save summary files
print("Saving summary files...")

# Save F statistics summary
if (nrow(all_fstats) > 0) {
  fwrite(all_fstats, "instrument_strength_summary.csv", quote = FALSE, row.names = FALSE)
  print(paste("Saved F statistics for", nrow(all_fstats), "analyses to instrument_strength_summary.csv"))
}

# Save harmonised data for supplementary table
if (nrow(all_harmonised) > 0) {
  # Rename columns for clarity in supplementary table
  colnames(all_harmonised) <- c("Gene", "Dataset", "CRC_Subtype", "SNP", 
                                "Effect_Allele", "Other_Allele", 
                                "Beta_Exposure", "SE_Exposure", "N_Exposure",
                                "Beta_Outcome", "SE_Outcome", "N_Outcome")
  
  fwrite(all_harmonised, "supplementary_harmonised_data.csv", quote = FALSE, row.names = FALSE)
  print(paste("Saved harmonised data for", nrow(all_harmonised), "SNPs across all analyses to supplementary_harmonised_data.csv"))
}

print("fin")