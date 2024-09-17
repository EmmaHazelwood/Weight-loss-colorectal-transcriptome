library(ieugwasr)
library(TwoSampleMR)
library(dplyr)
sink('ColonClumping.txt')

#We want to clump locally to avoid any server 503 problems


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
    effect_allele_col = "REF",
    other_allele_col = "ALT",
    eaf_col = "eaf.exposure"
  )
  clump<- ld_clump(
    dplyr::tibble(rsid=exposure_dat$SNP, pval=exposure_dat$pval.exposure, id=exposure_dat$id.exposure),
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = "EUR"
  )
  exposure_dat<-exposure_dat[exposure_dat$SNP %in% clump$rsid ,]
  write.csv(exposure_dat, paste(x[i], "BARC_eqtl","clumped.csv", sep="_"))
  }, error=function(e){cat("ERROR", conditionMessage(e),"\n")})
}

#Number 83 (ENSG00000139197) doesn't work - locally or using clump_data or ld_clump.
#Trying ENSG00000139197 in LDLink - downloaded R2 file - "r2_69956.txt"
#In R on Tom's server:
library(data.table)
library(dplyr)
x<-read_exposure_data(
  filename = "ENSG00000139197_colonBARC_eqtl.csv",
  sep = ",",
  snp_col = "variant_id",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "REF",
  other_allele_col = "ALT",
  eaf_col = "eaf.exposure"
)
y<-fread("r2_69956.txt")
x2<-arrange(x,pval.exposure)
x2<-x2[1,]
y2<-select(y,RS_number,rs1450962)
y3<-y2[which(y2$rs1450962<0.001),]
#no rows, so can only have this 1 SNP
write.csv(x2,"ENSG00000139197_BARC_eqtl_clumped.csv")
