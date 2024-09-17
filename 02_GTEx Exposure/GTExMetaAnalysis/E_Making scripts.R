for(k in 1:100){
  code <- '#Meta-analysis
library(meta)
library(data.table)
library(dplyr)

data<-fread(paste0("genes",j,".csv"))

results=data.frame(matrix(ncol=10,nrow=0))
ls<-unique(data$phenotype_id)

for (a in 1:length(ls)){
  print(ls[a])
  
data_2<-data[which(data$phenotype_id %in% ls[a]),]

#Split into SNPs
list<-unique(data_2$SNP)

for (i in 1:length(list)){
  print(list[i])

data_1<-data_2[which(data_2$SNP %in% list[i]),]
  
meta_res<-metagen(data_1$b, data_1$se, sm="OR", studlab=data_1$dataset)

#view forest plot, remove if not interested
#forest(metagen(data_1$b, data_1$se, sm="OR", studlab=data_1$dataset))

#save fixed-effects meta-analysis results as a dataframe
fixed<-data.frame( meta_res$TE.fixed, meta_res$seTE.fixed, meta_res$lower.fixed, meta_res$upper.fixed, meta_res$pval.fixed, meta_res$zval.fixed, meta_res$I2)

fixed$test<-"fixed"

fixed<-rename(fixed, c("beta"="meta_res.TE.fixed", "se"="meta_res.seTE.fixed", "LCI"="meta_res.lower.fixed",
                       
                       "UCI"="meta_res.upper.fixed", "pval"="meta_res.pval.fixed", "zval"="meta_res.zval.fixed",
                       
                       "I2"="meta_res.I2"))

#save random effects meta-analysis as a data.frame

random<-data.frame( meta_res$TE.random, meta_res$seTE.random, meta_res$lower.random, meta_res$upper.random, meta_res$pval.random, meta_res$zval.random, meta_res$I2)

random$test<-"random"

random<-rename(random, c("beta"="meta_res.TE.random", "se"="meta_res.seTE.random", "LCI"="meta_res.lower.random",
                         
                         "UCI"="meta_res.upper.random", "pval"="meta_res.pval.random", "zval"="meta_res.zval.random",
                         
                         "I2"="meta_res.I2"))

#combine fixed and random effects meta-analysis results into 1 dataframe

both<-rbind(fixed, random)
both$phenotype_id<-ls[a]
both$variant_id<-list[i]
results=rbind(results,both)
}
}
write.csv(results,paste0("GTExOverall",j,".csv"),row.names=FALSE,quote=FALSE,sep=",")'
  code <- gsub("j", k, code)
  file <- paste0("GTExmeta-analysis", k, ".R")
  writeLines(code, file)
  }
