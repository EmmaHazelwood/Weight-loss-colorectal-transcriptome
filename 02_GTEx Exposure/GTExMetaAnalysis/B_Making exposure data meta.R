#module add languages/R-4.0.3-gcc9.1.0
library(data.table)
library(dplyr)

CRC<- fread("/projects/MRC-IEU/research/projects/icep2/wp2/005/working/data/gecco/annotated/colon_CRC_GWAS_noUKBio_summary_stats_annotated.txt")

#colon eQTLs
x<-dir(pattern="_colonoverall_eqtl.csv")
x<-substr(x,1,15)
x<-unique(x)

for (i in 1:length(x)){
  y<-dir(pattern=x[i])
  print(x[i])
  sig<-read.csv(grep("overall_eqtl.csv",y,value=TRUE),sep=",")
  SNPs <- subset(CRC,chr_pos %in% sig$pos)
  write.csv(SNPs, paste(x[i], "overall","colon_outcome.csv", sep="_"))
  SNPsID<- select(SNPs,chr_pos,SNP)
  sigID<-merge(sig,SNPsID,by.x="pos",by.y="chr_pos")
  write.csv(sigID, paste(x[i], "overall","colon_exposure.csv", sep="_"))
}
