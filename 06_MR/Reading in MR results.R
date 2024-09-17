library(data.table)
library(tidyverse)

tbl_fread <- list.files(pattern = "MR_Results.csv") %>% map_df(~fread(.))
tbl_fread<-as.data.frame(tbl_fread)

tbl_fread<-merge(genes,tbl_fread,by.x="ensembl_gene_id",by.y="id.exposure")

tbl_fread$bh_p <- p.adjust(tbl_fread$"pval",method = "BH")
tbl_fread<-tbl_fread[order(tbl_fread$bh_p),]

#Need to swap BARCUVa directions
tbl_fread$b[which(tbl_fread$type=="BARCeQTL")]<-(-1*tbl_fread$b[which(tbl_fread$type=="BARCeQTL")])

write.csv(tbl_fread,"MR_Results.csv",quote = FALSE,row.names = FALSE)
write.csv(tbl_fread,"MR_Results.csv",quote = FALSE,row.names = FALSE)
