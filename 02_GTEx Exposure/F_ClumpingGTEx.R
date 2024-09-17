sink('ColonClumping.txt')

library(ieugwasr)
library(TwoSampleMR)
library(dplyr)
library(data.table)


#We want to clump locally to avoid any server 503 problems

#Sigmoid
x<-dir(pattern="_colonsigmoid_eqtl_snps.csv")
x<-substr(x,1,15)
x<-unique(x)


for (i in 1:length(x)){
  tryCatch({
y<-dir(pattern=x[i])
print(i)
print(x[i])
exposure_dat <- read_exposure_data(
filename = grep("colonsigmoid_eqtl_snps.csv",y,value=TRUE),
sep = ",",
snp_col = "rs_id_dbSNP151_GRCh38p7",
beta_col = "beta.exposure",
se_col = "se.exposure",
effect_allele_col = "effect_allele.exposure",
other_allele_col = "other_allele.exposure",
eaf_col = "eaf.exposure",
phenotype_col = "phenotype_id"
)
rm(clump)
clump<- ld_clump(
  dplyr::tibble(rsid=exposure_dat$SNP, pval=exposure_dat$pval.exposure, id=exposure_dat$id.exposure),
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "1000GenomesReferenceFiles/EUR"
)
if (exists("clump")==TRUE){
exposure_dat<-exposure_dat[exposure_dat$SNP %in% clump$rsid ,]
fwrite(exposure_dat, paste(x[i], "colonsigmoid_eqtl","clumped.csv", sep="_"))
} else {
exposure_dat<-exposure_dat[order(exposure_dat$pval.exposure),]
exposure_dat<-exposure_dat[1,]
fwrite(exposure_dat, paste(x[i], "colonsigmoid_eqtl","clumped.csv", sep="_"))
}
}, error=function(e){cat("ERROR", conditionMessage(e),"\n")})
}
#Transverse
x<-dir(pattern="_colontransverse_eqtl_snps.csv")
x<-substr(x,1,15)
x<-unique(x)


for (i in 1:length(x)){
  tryCatch({
    y<-dir(pattern=x[i])
    print(i)
    print(x[i])
    exposure_dat <- read_exposure_data(
      filename = grep("colontransverse_eqtl_snps.csv",y,value=TRUE),
      sep = ",",
      snp_col = "rs_id_dbSNP151_GRCh38p7",
      beta_col = "beta.exposure",
      se_col = "se.exposure",
      effect_allele_col = "effect_allele.exposure",
      other_allele_col = "other_allele.exposure",
      eaf_col = "eaf.exposure",
      phenotype_col = "phenotype_id"
    )
    rm(clump)
    clump<- ld_clump(
      dplyr::tibble(rsid=exposure_dat$SNP, pval=exposure_dat$pval.exposure, id=exposure_dat$id.exposure),
      plink_bin = genetics.binaRies::get_plink_binary(),
      bfile = "1000GenomesReferenceFiles/EUR"
    )
    if (exists("clump")==TRUE){
      exposure_dat<-exposure_dat[exposure_dat$SNP %in% clump$rsid ,]
      fwrite(exposure_dat, paste(x[i], "colontransverse_eqtl","clumped.csv", sep="_"))
    } else {
      exposure_dat<-exposure_dat[order(exposure_dat$pval.exposure),]
      exposure_dat<-exposure_dat[1,]
      fwrite(exposure_dat, paste(x[i], "colontransverse_eqtl","clumped.csv", sep="_"))
    }
  }, error=function(e){cat("ERROR", conditionMessage(e),"\n")})
}

#Total
x<-dir(pattern="_colontotal_eqtl_snps.csv")
x<-substr(x,1,15)
x<-unique(x)


for (i in 1:length(x)){
  tryCatch({
    y<-dir(pattern=x[i])
    print(i)
    print(x[i])
    exposure_dat <- read_exposure_data(
      filename = grep("colontotal_eqtl_snps.csv",y,value=TRUE),
      sep = ",",
      snp_col = "rs_id_dbSNP151_GRCh38p7",
      beta_col = "beta.exposure",
      se_col = "se.exposure",
      effect_allele_col = "effect_allele.exposure",
      other_allele_col = "other_allele.exposure",
      eaf_col = "eaf.exposure",
      phenotype_col = "phenotype_id"
    )
    rm(clump)
    clump<- ld_clump(
      dplyr::tibble(rsid=exposure_dat$SNP, pval=exposure_dat$pval.exposure, id=exposure_dat$id.exposure),
      plink_bin = genetics.binaRies::get_plink_binary(),
      bfile = "1000GenomesReferenceFiles/EUR"
    )
    if (exists("clump")==TRUE){
      exposure_dat<-exposure_dat[exposure_dat$SNP %in% clump$rsid ,]
      fwrite(exposure_dat, paste(x[i], "colontotal_eqtl","clumped.csv", sep="_"))
    } else {
      exposure_dat<-exposure_dat[order(exposure_dat$pval.exposure),]
      exposure_dat<-exposure_dat[1,]
      fwrite(exposure_dat, paste(x[i], "colontotal_eqtl","clumped.csv", sep="_"))
    }
  }, error=function(e){cat("ERROR", conditionMessage(e),"\n")})
}

sink()