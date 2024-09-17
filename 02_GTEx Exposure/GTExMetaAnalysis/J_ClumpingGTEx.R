sink('ColonTotalClumping.txt')

library(ieugwasr)
library(TwoSampleMR)
library(dplyr)
library(data.table)


#We want to clump locally to avoid any server 503 problems

x<-dir(pattern="_colontotal_eqtl.csv")
x<-substr(x,1,15)
x<-unique(x)
maf<-fread("/home/sw20203/MiniProject2/genexpression/ColonSigmoid.csv")
maf$phenotype_id<-substr(maf$phenotype_id,1,15)
maf<-select(maf,phenotype_id,variant_id,maf)
SNP<-fread("biomartResults.txt")

#Add in rsID 
for (i in 1:length(x)){
  if(exists("file2")=="TRUE") {
rm(file2)
  }
y<-dir(pattern=x[i])
print(i)
print(x[i])
file<-fread(grep("colontotal_eqtl.csv",y,value=TRUE))
file<-merge(file,maf,by=c("variant_id","phenotype_id"),all.x=TRUE)
file<-merge(file,SNP,by.x=c("chr","chrom_start"),by.y=c("chrom_strand","chrom_start"),all.x=TRUE)
#For rsIDs that were found by biomart, need to confirm that compatible SNP with the effect and other alleles - already true for myvariant SNPs as 
file2<-file[which(file$allele!="NA"),]
file<-file[which(file$allele=="NA"),]
library(stringr)
a<-str_split_fixed(file2$allele,pattern="/",n=Inf)
file2<-cbind(file2,a)
if(nrow(file2)>=1){
b<-ncol(file2)-22
} else {
b=0
 }
if(b==0){
"Not needed"

} else if(b==2){
file2<- filter(file2,effect_allele.exposure==V1 | effect_allele.exposure==V2)
file2<- filter(file2,other_allele.exposure==V1 | other_allele.exposure==V2)

} else if(b==3) {
file2<- filter(file2,effect_allele.exposure==V1 | effect_allele.exposure==V2 | effect_allele.exposure==V3)
file2<- filter(file2,other_allele.exposure==V1 | other_allele.exposure==V2 | other_allele.exposure==V3)

} else if(b==4) {
file2<- filter(file2,effect_allele.exposure==V1 | effect_allele.exposure==V2 | effect_allele.exposure==V3 | effect_allele.exposure==V4)
file2<- filter(file2,other_allele.exposure==V1 | other_allele.exposure==V2 | other_allele.exposure==V3 | other_allele.exposure==V4)
  
} else if(b==5) {
file2<- filter(file2,effect_allele.exposure==V1 | effect_allele.exposure==V2 | effect_allele.exposure==V3 | effect_allele.exposure==V4 | effect_allele.exposure==V5)
file2<- filter(file2,other_allele.exposure==V1 | other_allele.exposure==V2 | other_allele.exposure==V3 | other_allele.exposure==V4 | other_allele.exposure==V5)

} else if(b==6) {
file2<- filter(file2,effect_allele.exposure==V1 | effect_allele.exposure==V2 | effect_allele.exposure==V3 | effect_allele.exposure==V4 | effect_allele.exposure==V5 | effect_allele.exposure==V6)
file2<- filter(file2,other_allele.exposure==V1 | other_allele.exposure==V2 | other_allele.exposure==V3 | other_allele.exposure==V4 | other_allele.exposure==V5 | effect_allele.exposure==V6)

} else if(b==7) {
file2<- filter(file2,effect_allele.exposure==V1 | effect_allele.exposure==V2 | effect_allele.exposure==V3 | effect_allele.exposure==V4 | effect_allele.exposure==V5 | effect_allele.exposure==V6 | effect_allele.exposure==V7)
file2<- filter(file2,other_allele.exposure==V1 | other_allele.exposure==V2 | other_allele.exposure==V3 | other_allele.exposure==V4 | other_allele.exposure==V5 | effect_allele.exposure==V6 | effect_allele.exposure==V7)

} else if(b==8) {
file2<- filter(file2,effect_allele.exposure==V1 | effect_allele.exposure==V2 | effect_allele.exposure==V3 | effect_allele.exposure==V4 | effect_allele.exposure==V5 | effect_allele.exposure==V6 | effect_allele.exposure==V7 | effect_allele.exposure==V8)
file2<- filter(file2,other_allele.exposure==V1 | other_allele.exposure==V2 | other_allele.exposure==V3 | other_allele.exposure==V4 | other_allele.exposure==V5 | effect_allele.exposure==V6 | effect_allele.exposure==V7 | effect_allele.exposure==V8)

} else if(b==9) {
file2<- filter(file2,effect_allele.exposure==V1 | effect_allele.exposure==V2 | effect_allele.exposure==V3 | effect_allele.exposure==V4 | effect_allele.exposure==V5 | effect_allele.exposure==V6 | effect_allele.exposure==V7 | effect_allele.exposure==V8 | effect_allele.exposure==V9)
file2<- filter(file2,other_allele.exposure==V1 | other_allele.exposure==V2 | other_allele.exposure==V3 | other_allele.exposure==V4 | other_allele.exposure==V5 | effect_allele.exposure==V6 | effect_allele.exposure==V7 | effect_allele.exposure==V8 | effect_allele.exposure==V9)

} else if(b==9) {
file2<- filter(file2,effect_allele.exposure==V1 | effect_allele.exposure==V2 | effect_allele.exposure==V3 | effect_allele.exposure==V4 | effect_allele.exposure==V5 | effect_allele.exposure==V6 | effect_allele.exposure==V7 | effect_allele.exposure==V8 | effect_allele.exposure==V9 | effect_allele.exposure==V10)
file2<- filter(file2,other_allele.exposure==V1 | other_allele.exposure==V2 | other_allele.exposure==V3 | other_allele.exposure==V4 | other_allele.exposure==V5 | effect_allele.exposure==V6 | effect_allele.exposure==V7 | effect_allele.exposure==V8 | effect_allele.exposure==V9 | effect_allele.exposure==V10)
} else {
  print("not working")
}


file2<-select(file2,phenotype_id,refsnp_id,maf,chr,chrom_start,beta.exposure,se.exposure,pval,effect_allele.exposure,other_allele.exposure)
file<-select(file,phenotype_id,refsnp_id,maf,chr,chrom_start,beta.exposure,se.exposure,pval,effect_allele.exposure,other_allele.exposure)
file<-rbind(file,file2)


fwrite(file,paste(x[i], "colontotal_eqtl.csv", sep="_"))}

x<-x[-94]

for (i in 1:length(x)){
y<-dir(pattern=x[i])
print(i)
print(x[i])
exposure_dat <- read_exposure_data(
filename = grep("colontotal_eqtl.csv",y,value=TRUE),
sep = ",",
snp_col = "refsnp_id",
beta_col = "beta.exposure",
se_col = "se.exposure",
effect_allele_col = "effect_allele.exposure",
other_allele_col = "other_allele.exposure",
eaf_col = "maf",
phenotype_col = "phenotype_id"
)
exposure_dat<-exposure_dat[which(exposure_dat$pval.exposure<=0.000005),]
clump<- ld_clump(
dplyr::tibble(rsid=exposure_dat$SNP, pval=exposure_dat$pval.exposure, id=exposure_dat$id.exposure),
plink_bin = genetics.binaRies::get_plink_binary(),
bfile = "/projects/MRC_IEU/research/projects/icep2/wp2/034/working/data/1000GenomesReferenceFiles/EUR"
)
exposure_dat<-exposure_dat[exposure_dat$SNP %in% clump$rsid ,]
fwrite(exposure_dat, paste(x[i], "colontotal_eqtl","clumped.csv", sep="_"))
}

sink()

#For 94 (ENSG00000139197) no clumping results? Leaving out and trying again