sink("r2andFstats.txt")

library(dplyr)
library(data.table)
library(ieugwasr)
library(TwoSampleMR)
library(MRInstruments)
library(LDlinkR)
library(tidyr)


setwd("r2_and_F_stats")

#GTEx
#Total
print("GTEx")
y1<-dir("exposuredat/p0.000005",pattern="_colontotal_eqtl_clumped.csv")
x<-substr(y1,1,15)
x<-unique(x)

results<-data.frame(matrix(ncol=4))
colnames(results)<-c("Gene","Dataset","r2","F stat")
for (i in 1:length(x)){
  y<-y1[grep(x[i],y1)]
  print(i)
  print(x[i])
  
  exposure_dat <- read_exposure_data(
    filename = paste("exposuredat/p0.000005/",y,sep=""),
    snp_col = "SNP",
    sep=",",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure"
  )
  exposure_dat$samplesize.exposure<-686
  
  exposure_dat$num <- 2*(exposure_dat$beta.exposure^2)*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
  exposure_dat$den <- 2*(exposure_dat$beta.exposure^2)*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure) + ((exposure_dat$se.exposure^2)*2*exposure_dat$samplesize.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure))
  exposure_dat$pve <- exposure_dat$num/exposure_dat$den              
  pve_exposure_dat = sum(exposure_dat$pve)
  r <- pve_exposure_dat
  F=((pve_exposure_dat)*(exposure_dat$samplesize.exposure-2))/(1-pve_exposure_dat)
  res<-c(x[i],"ColonTotal",r,F)
  results<-rbind(results,res)
}

#Sigmoid
print("GTEx sigmoid")
y1<-dir("exposuredat/p0.000005",pattern="_colonsigmoid_eqtl_clumped.csv")
x<-substr(y1,1,15)
x<-unique(x)

for (i in 1:length(x)){
  y<-y1[grep(x[i],y1)]
  print(i)
  print(x[i])
  
  exposure_dat <- read_exposure_data(
    filename = paste("exposuredat/p0.000005/",y,sep=""),
    snp_col = "SNP",
    sep=",",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure"
  )
  exposure_dat$samplesize.exposure<-318
  
  exposure_dat$num <- 2*(exposure_dat$beta.exposure^2)*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
  exposure_dat$den <- 2*(exposure_dat$beta.exposure^2)*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure) + ((exposure_dat$se.exposure^2)*2*exposure_dat$samplesize.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure))
  exposure_dat$pve <- exposure_dat$num/exposure_dat$den              
  pve_exposure_dat = sum(exposure_dat$pve)
  r <- pve_exposure_dat
  F=((pve_exposure_dat)*(exposure_dat$samplesize.exposure-2))/(1-pve_exposure_dat)
  res<-c(x[i],"Colonsigmoid",r,F)
  results<-rbind(results,res)
}

#transverse
print("GTEx transverse")
y1<-dir("exposuredat/p0.000005",pattern="_colontransverse_eqtl_clumped.csv")
x<-substr(y1,1,15)
x<-unique(x)

for (i in 1:length(x)){
  y<-y1[grep(x[i],y1)]
  print(i)
  print(x[i])
  
  exposure_dat <- read_exposure_data(
    filename = paste("exposuredat/p0.000005/",y,sep=""),
    snp_col = "SNP",
    sep=",",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure"
  )
  exposure_dat$samplesize.exposure<-368
  
  exposure_dat$num <- 2*(exposure_dat$beta.exposure^2)*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
  exposure_dat$den <- 2*(exposure_dat$beta.exposure^2)*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure) + ((exposure_dat$se.exposure^2)*2*exposure_dat$samplesize.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure))
  exposure_dat$pve <- exposure_dat$num/exposure_dat$den              
  pve_exposure_dat = sum(exposure_dat$pve)
  r <- pve_exposure_dat
  F=((pve_exposure_dat)*(exposure_dat$samplesize.exposure-2))/(1-pve_exposure_dat)
  res<-c(x[i],"Colontransverse",r,F)
  results<-rbind(results,res)
}

#BARC
print("BARC")
y1<-dir("BARC-UVA",pattern="_BARC_eqtl_clumped.csv")
x<-substr(y1,1,15)
x<-unique(x)

for (i in 1:length(x)){
  y<-y1[grep(x[i],y1)]
  print(i)
  print(x[i])
  
  exposure_dat <- read_exposure_data(
    filename = paste("BARC-UVA/",y,sep=""),
    snp_col = "SNP",
    sep=",",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure"
  )
  exposure_dat$samplesize.exposure<-445
  
  exposure_dat$num <- 2*(exposure_dat$beta.exposure^2)*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure)
  exposure_dat$den <- 2*(exposure_dat$beta.exposure^2)*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure) + ((exposure_dat$se.exposure^2)*2*exposure_dat$samplesize.exposure*exposure_dat$eaf.exposure*(1-exposure_dat$eaf.exposure))
  exposure_dat$pve <- exposure_dat$num/exposure_dat$den              
  pve_exposure_dat = sum(exposure_dat$pve)
  r <- pve_exposure_dat
  F=((pve_exposure_dat)*(exposure_dat$samplesize.exposure-2))/(1-pve_exposure_dat)
  res<-c(x[i],"BARC",r,F)
  results<-rbind(results,res)
}

#eQTLGen
exposure_dat<-fread("eQTLGen/eQTLGenGenes.csv")

#Filter for cis SNPs
exposure_dat$gene<-exposure_dat$ensembl_gene_id
exposure_dat$chromosome_name<-as.numeric(exposure_dat$chromosome_name)
exposure_dat<-exposure_dat[which((exposure_dat$chr.exposure=exposure_dat$chromosome_name) & (exposure_dat$pos.exposure <=exposure_dat$cisend)& (exposure_dat$pos.exposure >=exposure_dat$cisstart)),]


for (i in 1:length(unique(exposure_dat$gene))){
  exposure_dat2 <- exposure_dat[which(exposure_dat$gene==unique(exposure_dat$gene)[i]),]  
  exposure_dat2$num <- 2*(exposure_dat2$beta.exposure^2)*exposure_dat2$eaf.exposure*(1-exposure_dat2$eaf.exposure)
  exposure_dat2$den <- 2*(exposure_dat2$beta.exposure^2)*exposure_dat2$eaf.exposure*(1-exposure_dat2$eaf.exposure) + ((exposure_dat2$se.exposure^2)*2*exposure_dat2$samplesize.exposure*exposure_dat2$eaf.exposure*(1-exposure_dat2$eaf.exposure))
  exposure_dat2$pve <- exposure_dat2$num/exposure_dat2$den              
  pve_exposure_dat2 = sum(exposure_dat2$pve)
  r <- pve_exposure_dat2
  F=((pve_exposure_dat2)*(min(exposure_dat2$samplesize.exposure)-2))/(1-pve_exposure_dat2)
  res<-c(unique(exposure_dat$gene)[i],"eQTLGen",r,F)
  results<-rbind(results,res)
}

fwrite(results,"r2_and_F_stats_results.csv",quote=FALSE,row.names = FALSE)

sink()