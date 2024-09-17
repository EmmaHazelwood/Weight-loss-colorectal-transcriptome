library(dplyr)
library(data.table)
library(ieugwasr)
library(TwoSampleMR)
library(MRInstruments)
library(LDlinkR)
library(tidyr)

sink("PerformingMRoverall.txt")

CRC<-fread("overall_CRC_GWAS_noUKBio_summary_stats_annotated.txt")

#GTEx
print("GTEx")
x<-dir(pattern="_colontotal_eqtl_clumped.csv")
x<-substr(x,1,15)
x<-unique(x)

for (i in 1:length(x)){
  y<-dir(pattern=x[i])
  print(i)
  print(x[i])
  
  exposure_dat <- read_exposure_data(
    filename = grep("colontotal_eqtl_clumped.csv",y,value=TRUE),
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
  
  if (length(intersect(exposure_dat$SNP,CRC$SNP)!=0)){
    outcome_dat <- read_outcome_data(
      snps = exposure_dat$SNP,
      filename = "Proxy/tempoverall.csv",
      sep = ",",
      snp_col = "SNP",
      beta_col = "Effect",
      se_col = "StdErr",
      effect_allele_col = "Allele1",
      other_allele_col = "Allele2",
      eaf_col = "Freq1"
    )
    
    
    dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
    dat<- subset(dat, (as.character(effect_allele.exposure) == as.character(effect_allele.outcome)) & (as.character(other_allele.exposure) == as.character(other_allele.outcome)))
    dat<- subset(dat, (as.character(palindromic) == "FALSE"))
    
    if (nrow(dat)==0){
      print("No rows")
    } else if (nrow(dat)==1){
      res <- mr(dat, method_list="mr_wald_ratio")
      res$id.exposure <- paste(x[i])
      res$type <- "totaleQTL"
      res$id.outcome <- "overall"
      write.csv(res, paste(x[i], "total","overall_MR_Results.csv", sep="_"),quote = FALSE)
    } else {
      res <- mr(dat, method_list="mr_ivw")
      res$id.exposure <- paste(x[i])
      res$type <- "totaleQTL"
      res$id.outcome <- "overall"
      write.csv(res, paste(x[i], "total","overall_MR_Results.csv", sep="_"),quote = FALSE)
    }
  } else {
    print("No SNPs in common")
  }
}

#Sigmoid
print("Sigmoid")
x<-dir(pattern="_colonsigmoid_eqtl_clumped.csv")
x<-substr(x,1,15)
x<-unique(x)

for (i in 1:length(x)){
  y<-dir(pattern=x[i])
  print(i)
  print(x[i])
  
  exposure_dat <- read_exposure_data(
    filename = grep("colonsigmoid_eqtl_clumped.csv",y,value=TRUE),
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
  
  if (length(intersect(exposure_dat$SNP,CRC$SNP)!=0)){
    outcome_dat <- read_outcome_data(
      snps = exposure_dat$SNP,
      filename = "Proxy/tempoverall.csv",
      sep = ",",
      snp_col = "SNP",
      beta_col = "Effect",
      se_col = "StdErr",
      effect_allele_col = "Allele1",
      other_allele_col = "Allele2",
      eaf_col = "Freq1"
    )
    
    
    dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
    dat<- subset(dat, (as.character(effect_allele.exposure) == as.character(effect_allele.outcome)) & (as.character(other_allele.exposure) == as.character(other_allele.outcome)))
    dat<- subset(dat, (as.character(palindromic) == "FALSE"))
    
    if (nrow(dat)==0){
      print("No rows")
    } else if (nrow(dat)==1){
      res <- mr(dat, method_list="mr_wald_ratio")
      res$id.exposure <- paste(x[i])
      res$type <- "sigmoideQTL"
      res$id.outcome <- "overall"
      write.csv(res, paste(x[i], "sigmoid","overall_MR_Results.csv", sep="_"),quote = FALSE)
    } else {
      res <- mr(dat, method_list="mr_ivw")
      res$id.exposure <- paste(x[i])
      res$type <- "sigmoideQTL"
      res$id.outcome <- "overall"
      write.csv(res, paste(x[i], "sigmoid","overall_MR_Results.csv", sep="_"),quote = FALSE)
    }
  } else {
    print("No SNPs in common")
  }
}

#Transverse
print("Transverse")
x<-dir(pattern="_colontransverse_eqtl_clumped.csv")
x<-substr(x,1,15)
x<-unique(x)
x<-x[-72]

for (i in 1:length(x)){
  y<-dir(pattern=x[i])
  print(i)
  print(x[i])
  
  exposure_dat <- read_exposure_data(
    filename = grep("colontransverse_eqtl_clumped.csv",y,value=TRUE),
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
  
  if (length(intersect(exposure_dat$SNP,CRC$SNP)!=0)){
    outcome_dat <- read_outcome_data(
      snps = exposure_dat$SNP,
      filename = "Proxy/tempoverall.csv",
      sep = ",",
      snp_col = "SNP",
      beta_col = "Effect",
      se_col = "StdErr",
      effect_allele_col = "Allele1",
      other_allele_col = "Allele2",
      eaf_col = "Freq1"
    )
    
    
    dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
    dat<- subset(dat, (as.character(effect_allele.exposure) == as.character(effect_allele.outcome)) & (as.character(other_allele.exposure) == as.character(other_allele.outcome)))
    dat<- subset(dat, (as.character(palindromic) == "FALSE"))
    
    if (nrow(dat)==0){
      print("No rows")
    } else if (nrow(dat)==1){
      res <- mr(dat, method_list="mr_wald_ratio")
      res$id.exposure <- paste(x[i])
      res$type <- "transverseeQTL"
      res$id.outcome <- "overall"
      write.csv(res, paste(x[i], "transverse","overall_MR_Results.csv", sep="_"),quote = FALSE)
    } else {
      res <- mr(dat, method_list="mr_ivw")
      res$id.exposure <- paste(x[i])
      res$type <- "transverseeQTL"
      res$id.outcome <- "overall"
      write.csv(res, paste(x[i], "transverse","overall_MR_Results.csv", sep="_"),quote = FALSE)
    }
  } else {
    print("No SNPs in common")
  }
}

#BarcUVa-Seq
print("BARC")
x<-dir(pattern="BARC_eqtl_clumped.csv")
x<-substr(x,1,15)
x<-unique(x)

for (i in 1:length(x)){
  y<-dir(pattern=x[i])
  print(i)
  print(x[i])
  
  exposure_dat <- read_exposure_data(
    filename = grep("BARC_eqtl_clumped.csv",y,value=TRUE),
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
    snps = exposure_dat$SNP,
    filename = "Proxy/tempoverall.csv",
    sep = ",",
    snp_col = "SNP",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1"
  )
  
  
  dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
  dat<- subset(dat, (as.character(effect_allele.exposure) == as.character(effect_allele.outcome)) & (as.character(other_allele.exposure) == as.character(other_allele.outcome)))
  dat<- subset(dat, (as.character(palindromic) == "FALSE"))
  
  if (nrow(dat)==0){
    print("No rows")
  } else if (nrow(dat)==1){
    res <- mr(dat, method_list="mr_wald_ratio")
    res$id.exposure <- paste(x[i])
    res$type <- "BARCeQTL"
    res$id.outcome <- "overall"
    write.csv(res, paste(x[i], "BARC","overall_MR_Results.csv", sep="_"),quote = FALSE)
  } else {
    res <- mr(dat, method_list="mr_ivw")
    res$id.exposure <- paste(x[i])
    res$type <- "BARCeQTL"
    res$id.outcome <- "overall"
    write.csv(res, paste(x[i], "BARC","overall_MR_Results.csv", sep="_"),quote = FALSE)
  }
}


#eQTLGen
print("eQTLGen")
exposure_dat<-fread("eQTLGenGenes.csv")

#Filter for cis SNPs
exposure_dat$gene<-exposure_dat$ensembl_gene_id
exposure_dat$chromosome_name<-as.numeric(exposure_dat$chromosome_name)
exposure_dat<-exposure_dat[which((exposure_dat$chr.exposure=exposure_dat$chromosome_name) & (exposure_dat$pos.exposure <=exposure_dat$cisend)& (exposure_dat$pos.exposure >=exposure_dat$cisstart)),]


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
      ea<-distinct(exposure_dat[which(exposure_dat$SNP==j),9])
      oa<-distinct(exposure_dat[which(exposure_dat$SNP==j),10])
      Correlated_Alleles<-data.frame(proxy3$Correlated_Alleles)
      Correlated_Alleles$proxy3.Correlated_Alleles<- as.character(Correlated_Alleles$proxy3.Correlated_Alleles)
      Correlated_Alleles<-Correlated_Alleles %>% separate(proxy3.Correlated_Alleles, c("a1exp", "a1out", "a2exp", "a2out"))
      Correlated_Alleles<-na.omit(Correlated_Alleles)
      if (nrow(Correlated_Alleles)!=0){
        proxy4<-cbind(proxy3,Correlated_Alleles)
        proxy4<-cbind(proxy4,ea)
        proxy4<-cbind(proxy4,oa)
        
        if (proxy4$a1exp==proxy4$effect_allele.exposure){
          proxy4$REF<-proxy4$a1out
          proxy4$ALT<-proxy4$a2out
        } else if(proxy4$a2exp==proxy4$effect_allele.exposure){
          proxy4$REF<-proxy4$a2out
          proxy4$ALT<-proxy4$a1out
        } else {
          proxy4$REF<-"NA"
          proxy4$ALT<-"NA"
        }
        proxy5<-select(proxy4,missing_snp,RS_Number,REF,ALT)
        proxySNPs<-rbind(proxySNPs,proxy5)
        proxySNPs<-distinct(proxySNPs)
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
      exposure_dat2<-exposure_dat2[,1:21]
      exposure_dat<-rbind(exposure_dat,exposure_dat2)
    }}}


outcome_dat <- read_outcome_data(
  snps = exposure_dat$SNP,
  filename = "Proxy/tempoverall.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1"
)

dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
dat<- subset(dat, (as.character(effect_allele.exposure) == as.character(effect_allele.outcome)) & (as.character(other_allele.exposure) == as.character(other_allele.outcome)))
dat<- subset(dat, (as.character(palindromic) == "FALSE"))
dat<-distinct(dat)

genes<-unique(dat$gene)

for (i in 1:length(genes)){
  print(i)
  gene<-genes[i]
dat2<-dat[which(dat$gene==gene),]
if (nrow(dat2)==0){
  print("No rows")
} else if (nrow(dat2)==1){
  res <- mr(dat2, method_list="mr_wald_ratio")
  res$id.exposure <- paste(genes[i])
  res$type <- "eQTLGen"
  res$id.outcome <- "overall"
  write.csv(res, paste(genes[i], "eQTLGen","overall_MR_Results.csv", sep="_"),quote = FALSE)
} else {
  res <- mr(dat2, method_list="mr_ivw")
  res$id.exposure <- paste(genes[i])
  res$type <- "eQTLGen"
  res$id.outcome <- "overall"
  write.csv(res, paste(genes[i], "eQTLGen","overall_MR_Results.csv", sep="_"),quote = FALSE)
}}

sink()