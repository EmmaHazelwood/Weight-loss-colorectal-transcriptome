library(ieugwasr)
library(TwoSampleMR)
library(dplyr)

x<-dir(pattern="_colonBARC_eqtl.csv")
x<-substr(x,1,15)
x<-unique(x)

for (i in 1:length(x)){
  tryCatch({
  y<-dir(pattern=x[i])
  print(i)
  print(x[i])
  exposure_dat <- read_exposure_data(
    filename = grep("BARC_eqtl.csv",y,value=TRUE),
    sep = ",",
    snp_col = "variant_id",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    eaf_col = "eaf.exposure"
  )
  clump<- ld_clump(
    dplyr::tibble(rsid=exposure_dat$SNP, pval=exposure_dat$pval.exposure, id=exposure_dat$id.exposure),
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = "data/1000GenomesReferenceFiles/EUR"
  )
  exposure_dat<-exposure_dat[exposure_dat$SNP %in% clump$rsid ,]
  write.csv(exposure_dat, paste(x[i], "BARC_eqtl","clumped.csv", sep="_"))
  }, error=function(e){cat("ERROR", conditionMessage(e),"\n")})
}
