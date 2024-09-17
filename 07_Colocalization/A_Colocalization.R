sink("Coloc.txt")

library(data.table)
library(ieugwasr)
library(dplyr)
library(TwoSampleMR)
library(coloc)
library(tidyverse)
overall_samplesize<-98715
proximal_samplesize<-57515
distal_samplesize<-55978
rectal_samplesize<-57249
left_samplesize<-70103
colon_samplesize<-71835

#Get a list of analyses needed to perform
df<-fread("/MiniProject2/MR_Results_Filtered.csv")
#FDR
df$bh_p =
  p.adjust(df$"pval",
           method = "BH")
df<-df[which(df$bh_p<0.05),]

#eQTLGen
eQTL<-df[which(df$type=="eQTLGen"),]

#overall
overall_eQTL<-eQTL[which(eQTL$id.outcome=="overall"),]
exposure_dat<-extract_instruments(outcomes=paste("eqtl-a-",overall_eQTL$id.exposure,sep=""), p1=1e-06)



outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)

dat<-unique(setDT(dat)[order(pval.exposure)], by = "id.exposure")

dat <- dat[match(unique(dat$gene), dat$gene),]

overall_results<-data.frame(matrix(ncol=5,nrow=0))

for (i in 1:nrow(dat)){
  print(i)
  x<-dat[i,]
  print(x$id.exposure)
  print(x$SNP)
  x$window_start<-x$pos.exposure-1000000
  x$window_end<-x$pos.exposure+1000000
  variants<-paste(x$chr.exposure,x$window_start,sep=":")
  variants<-paste(variants,x$window_end,sep="-")
  snps<-associations(id=x$id.exposure,variants=variants)
  outcome_dat <- read_outcome_data(
    snps = snps$rsid,
    filename = "overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
    sep = " ",
    snp_col = "SNP",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1"
  )
  fwrite(snps,"snps.csv",row.names = FALSE,quote=FALSE)
  snps <- read_exposure_data(
    filename = "snps.csv",
    snp_col = "rsid",
    sep=",",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "ea",
    other_allele_col = "nea",
    eaf_col = "eaf"
  )
  snps <- harmonise_data(exposure_dat = snps, outcome_dat = outcome_dat)
  snps$id.outcome<-"overall"
  write.csv(snps,paste("eQTLGen_Locus_Zoom_",x$id.exposure,"overall",x$SNP,sep=""),row.names = FALSE,quote=FALSE)
  
  dataset1 <- list(beta=snps$beta.exposure, varbeta = snps$se.exposure^2, MAF=snps$eaf.exposure, type = "quant", N = 31300)
  dataset2 <- list(beta=snps$beta.outcome, varbeta = snps$se.outcome^2, MAF=snps$eaf.outcome, type = "quant", N = overall_samplesize)
  result <- coloc.abf(dataset1, dataset2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  df2 <- data.frame(matrix(unlist(result$summary), nrow = 1, byrow = T))
  names(df2) <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
  df2$snp<-x$SNP
  df2$gene<-x$id.exposure
  df2$site<-"overall"
  overall_results<-rbind(overall_results,df2)
}
write.csv(overall_results,"Overall_eQTLGen_Coloc_Results.csv",row.names = FALSE,quote=FALSE)

#BARCeQTL
eQTL<-df[which(df$type=="BARCeQTL"),]
x<-dir("BARC-UVA",pattern="clumped.csv")
z<-substr(x,1,15)
colon<-fread("BARC-UVA/barcuvaseq.eqtls.allpairs.sumstats.txt")
colon$id.exposure<-substr(colon$gene_id,1,15)

#overall
CRC<-fread("overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
overall_eQTL<-eQTL[which(eQTL$id.outcome=="overall"),]
z1<-intersect(z,overall_eQTL$id.exposure)
overall_results<-data.frame(matrix(nrow=0,ncol=9))

for (i in 1:length(x)){
  y<-z1[i]
  print(i)
  print(y)
  
  exposure_dat <- read_exposure_data(
    filename = paste("BARC-UVA/",y,"_BARC_eqtl_clumped.csv",sep=""),
    snp_col = "SNP",
    sep=",",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure"
  )
  
  snps_for_proxy <- setdiff(exposure_dat$SNP,CRC$SNP)
  
  if (length(snps_for_proxy)!=0){
    proxySNPs=data.frame(matrix(ncol=4,nrow=0))
    for(j in snps_for_proxy) {
      print(j)
      proxy <- LDproxy(j, pop = "EUR", r2d = "r2", token = "a563a6f57fcc")
      proxy2 <- proxy %>% dplyr::filter(R2 > 0.8) %>% # minimum R2 0.8
        dplyr::filter(RS_Number != ".") %>% # remove ones without rsID
        dplyr::filter(RS_Number %in% CRC$SNP) # only keep SNPs that are in CRC data
      if (length(proxy2)!=0){
        proxy3 <- proxy2[1,] # keep rsID at the top with the highest LD
        proxy3$missing_snp <- j # label the SNP we were trying to proxy so we can merge back in when creating summary dat
        ea<-exposure_dat[which(exposure_dat$SNP==j),5]
        oa<-exposure_dat[which(exposure_dat$SNP==j),6]
        Correlated_Alleles<-data.frame(proxy3$Correlated_Alleles)
        Correlated_Alleles$proxy3.Correlated_Alleles<- as.character(Correlated_Alleles$proxy3.Correlated_Alleles)
        Correlated_Alleles<-Correlated_Alleles %>% separate(proxy3.Correlated_Alleles, c("a1exp", "a1out", "a2exp", "a2out"))
        Correlated_Alleles<-na.omit(Correlated_Alleles)
        if (nrow(Correlated_Alleles)!=0){
          proxy4<-cbind(proxy3,Correlated_Alleles)
          proxy4<-cbind(proxy4,ea)
          proxy4<-cbind(proxy4,oa)
          
          if (proxy4$a1exp==proxy4$ea){
            proxy4$REF<-proxy4$a1out
            proxy4$ALT<-proxy4$a2out
          } else if(proxy4$a2exp==proxy4$ea){
            proxy4$REF<-proxy4$a2out
            proxy4$ALT<-proxy4$a1out
          } else {
            proxy4$REF<-"NA"
            proxy4$ALT<-"NA"
          }
          proxy5<-select(proxy4,missing_snp,RS_Number,REF,ALT)
          proxySNPs<-rbind(proxySNPs,proxy5)
        }
      } else{
        proxySNPs<-data.frame(matrix(ncol=17,nrow=0))
        colnames(proxySNPs)<-"missing_snp" 
      }
      
      
      
      exposure_datproxies<-exposure_dat[which(exposure_dat$SNP %in% snps_for_proxy),]
      
      if (nrow(proxySNPs)!=0){
        
        exposure_dat2<-merge(exposure_datproxies,proxySNPs,by.x="SNP",by.y="missing_snp")
        exposure_dat2$SNP<-exposure_dat2$RS_Number
        exposure_dat2$effect_allele.exposure<-exposure_dat2$REF
        exposure_dat2$other_allele.exposure<-exposure_dat2$ALT
        exposure_dat2<-exposure_dat2[,1:12]
        exposure_dat<-rbind(exposure_dat,exposure_dat2)
      }}}
  
  
  outcome_dat <- read_outcome_data(
    snps = exposure_dat$rsid,
    filename = "overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
    sep = " ",
    snp_col = "SNP",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1"
  )
  
  
  dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)

  p<-unique(setDT(dat)[order(pval.exposure)])
  p$pos.exposure<-CRC$chrom_start[which(CRC$SNP==p$SNP)]
  p$window_start<-p$pos.exposure-1000000
  p$window_end<-p$pos.exposure+1000000
  p$chr.exposure<-CRC$chr_name[which(CRC$SNP==p$SNP)]
  snps<-colon[which(colon$id.exposure==y & colon$chr==p$chr.exposure & colon$pos>=p$window_start & colon$pos<=p$window_end),]
  p$id.exposure<-y
  outcome_dat <- read_outcome_data(
    snps = snps$variant_id,
    filename = "overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
    sep = " ",
    snp_col = "SNP",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1"
  )
  fwrite(snps,"snps.csv",row.names = FALSE,quote=FALSE)
  snps <- read_exposure_data(
    filename = "snps.csv",
    snp_col = "variant_id",
    sep=",",
    beta_col = "slope",
    se_col = "slope_se",
    effect_allele_col = "REF",
    other_allele_col = "ALT",
    eaf_col = "maf"
  )
  snps <- harmonise_data(exposure_dat = snps, outcome_dat = outcome_dat)
  p$id.outcome<-"overall"
  write.csv(snps,paste("BarcUVa-Seq_Locus_Zoom_",p$id.exposure,p$id.outcome,p$SNP,sep=""),row.names = FALSE,quote=FALSE)
  dataset1 <- list(beta=snps$beta.exposure, varbeta = snps$se.exposure^2, MAF=snps$eaf.exposure, type = "quant", N = 445)
  dataset2 <- list(beta=snps$beta.outcome, varbeta = snps$se.outcome^2, MAF=snps$eaf.outcome, type = "quant", N = overall_samplesize)
  result <- coloc.abf(dataset1, dataset2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  df2 <- data.frame(matrix(unlist(result$summary), nrow = 1, byrow = T))
  names(df2) <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
  df2$snp<-dat$SNP
  df2$gene<-y
  df2$site<-"overall"
  overall_results<-rbind(overall_results,df2)
}
write.csv(overall_results,"Overall_BARCeQTL_Coloc_Results.csv",row.names = FALSE,quote=FALSE)

#GTEx
eQTL<-df[which(df$type=="totaleQTL"),]
x<-dir("GTExMetaAnalysis",pattern=eQTL$id.exposure)
x<-x[grep("clumped",x)]
z<-substr(x,1,15)
colon<-fread("GTExMetaAnalysis/GTExMetaAnalysisResultsFixed.csv")
colon$id.exposure<-substr(colon$phenotype_id,1,15)

#overall
CRC<-fread("overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt")
overall_eQTL<-eQTL[which(eQTL$id.outcome=="overall"),]
z1<-intersect(z,overall_eQTL$id.exposure)
overall_results<-data.frame(matrix(nrow=0,ncol=9))

for (i in 1:length(x)){
  y<-z1[i]
  print(i)
  print(z1[i])
  
  exposure_dat <- read_exposure_data(
    filename = paste("GTExMetaAnalysis/",x[i],sep=""),
    snp_col = "SNP",
    sep=",",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure"
  )
  
  snps_for_proxy <- setdiff(exposure_dat$SNP,CRC$SNP)
  
  if (length(snps_for_proxy)!=0){
    proxySNPs=data.frame(matrix(ncol=4,nrow=0))
    for(j in snps_for_proxy) {
      print(j)
      proxy <- LDproxy(j, pop = "EUR", r2d = "r2", token = "a563a6f57fcc")
      proxy2 <- proxy %>% dplyr::filter(R2 > 0.8) %>% # minimum R2 0.8
        dplyr::filter(RS_Number != ".") %>% # remove ones without rsID
        dplyr::filter(RS_Number %in% CRC$SNP) # only keep SNPs that are in CRC data
      if (length(proxy2)!=0){
        proxy3 <- proxy2[1,] # keep rsID at the top with the highest LD
        proxy3$missing_snp <- j # label the SNP we were trying to proxy so we can merge back in when creating summary dat
        ea<-exposure_dat[which(exposure_dat$SNP==j),5]
        oa<-exposure_dat[which(exposure_dat$SNP==j),6]
        Correlated_Alleles<-data.frame(proxy3$Correlated_Alleles)
        Correlated_Alleles$proxy3.Correlated_Alleles<- as.character(Correlated_Alleles$proxy3.Correlated_Alleles)
        Correlated_Alleles<-Correlated_Alleles %>% separate(proxy3.Correlated_Alleles, c("a1exp", "a1out", "a2exp", "a2out"))
        Correlated_Alleles<-na.omit(Correlated_Alleles)
        if (nrow(Correlated_Alleles)!=0){
          proxy4<-cbind(proxy3,Correlated_Alleles)
          proxy4<-cbind(proxy4,ea)
          proxy4<-cbind(proxy4,oa)
          
          if (proxy4$a1exp==proxy4$ea){
            proxy4$REF<-proxy4$a1out
            proxy4$ALT<-proxy4$a2out
          } else if(proxy4$a2exp==proxy4$ea){
            proxy4$REF<-proxy4$a2out
            proxy4$ALT<-proxy4$a1out
          } else {
            proxy4$REF<-"NA"
            proxy4$ALT<-"NA"
          }
          proxy5<-select(proxy4,missing_snp,RS_Number,REF,ALT)
          proxySNPs<-rbind(proxySNPs,proxy5)
        }
      } else{
        proxySNPs<-data.frame(matrix(ncol=17,nrow=0))
        colnames(proxySNPs)<-"missing_snp" 
      }
      
      
      
      exposure_datproxies<-exposure_dat[which(exposure_dat$SNP %in% snps_for_proxy),]
      
      if (nrow(proxySNPs)!=0){
        
        exposure_dat2<-merge(exposure_datproxies,proxySNPs,by.x="SNP",by.y="missing_snp")
        exposure_dat2$SNP<-exposure_dat2$RS_Number
        exposure_dat2$effect_allele.exposure<-exposure_dat2$REF
        exposure_dat2$other_allele.exposure<-exposure_dat2$ALT
        exposure_dat2<-exposure_dat2[,1:12]
        exposure_dat<-rbind(exposure_dat,exposure_dat2)
      }}}
  
  
outcome_dat <- read_outcome_data(
    snps = exposure_dat$rsid,
    filename = "overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
    sep = " ",
    snp_col = "SNP",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1"
  )
  
  
  dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
  
  

  p<-unique(setDT(dat)[order(pval.exposure)], by = "id.exposure")
    p$pos.exposure<-CRC$chrom_start[which(CRC$SNP==p$SNP)]
    p$window_start<-p$pos.exposure-1000000
    p$window_end<-p$pos.exposure+1000000
    p$chr.exposure<-CRC$chr_name[which(CRC$SNP==p$SNP)]
    snps<-colon[which(colon$id.exposure==y),]
    colon_snp<-snps %>% separate(variant_id, c("CHR", "BP", "other_allele.exposure", "effect_allele.exposure"))
    colon_snp$chr<- as.numeric(gsub("chr", "", colon_snp$CHR)) #remove the "chr" letters
    snps<-colon_snp[which(colon_snp$id.exposure==y & colon_snp$chr==p$chr.exposure & colon_snp$BP>=p$window_start & colon_snp$BP<=p$window_end),]
    
    #Get rsIDs
    snp_info<-fread("GTEx-Colon/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt")
    snps$BP<-as.numeric(snps$BP)
    snps<-merge(snps,snp_info,by.x=c("CHR","BP"),by.y=c("chr","variant_pos"))
    
    p$id.exposure<-y
    outcome_dat <- read_outcome_data(
      snps = snps$rs_id_dbSNP151_GRCh38p7,
      filename = "overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt",
      sep = " ",
      snp_col = "SNP",
      beta_col = "Effect",
      se_col = "StdErr",
      effect_allele_col = "Allele1",
      other_allele_col = "Allele2",
      eaf_col = "Freq1"
    )
    fwrite(snps,"snps.csv",row.names = FALSE,quote=FALSE)
    snps <- read_exposure_data(
      filename = "snps.csv",
      snp_col = "rs_id_dbSNP151_GRCh38p7",
      sep=",",
      beta_col = "beta",
      se_col = "se",
      effect_allele_col = "ref",
      other_allele_col = "alt",
      eaf_col = "maf"
    )
    snps <- harmonise_data(exposure_dat = snps, outcome_dat = outcome_dat)
    p$id.outcome<-"overall"
    write.csv(snps,paste("GTEx_Locus_Zoom_",p$id.exposure,p$id.outcome,p$SNP,sep=""),row.names = FALSE,quote=FALSE)
    dataset1 <- list(beta=snps$beta.exposure, varbeta = snps$se.exposure^2, MAF=snps$eaf.outcome, type = "quant", N = 31684)
    dataset2 <- list(beta=snps$beta.outcome, varbeta = snps$se.outcome^2, MAF=snps$eaf.outcome, type = "quant", N = overall_samplesize)
    result <- coloc.abf(dataset1, dataset2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
    df2 <- data.frame(matrix(unlist(result$summary), nrow = 1, byrow = T))
    names(df2) <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
    df2$snp<-dat$SNP
    df2$gene<-y
    df2$site<-"overall"
    overall_results<-rbind(overall_results,df2)
}
write.csv(overall_results,"Overall_GTExeQTL_Coloc_Results.csv",row.names = FALSE,quote=FALSE)



sink()