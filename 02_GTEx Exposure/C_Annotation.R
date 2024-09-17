#Illumina® HumanHT-12 v4 Expression BeadChip Assay
library(dplyr)
library(data.table)

res<-fread("Genesforloop.csv")

ref <- read.delim("Ensembl/GPL10558-50081.txt", quote="")

res$ID <- res$Probe_Id
res$FC<-exp(res$logFC)

mergeall <- merge(ref,res,by="ID")

merge <- dplyr::select(mergeall,ID,ILMN_Gene,Entrez_Gene_ID,Chromosome,Synonyms,FC,P.Value)

write.csv(merge,"Merge.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)

#Get Ensembl ID
library("biomaRt")                                                                                                                   
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
entrezgene = data.frame(merge$Entrez_Gene_ID)
entrezgene <- na.omit(entrezgene)
genes <- getBM(filters="entrezgene_id", attributes=c("ensembl_gene_id","entrezgene_id","start_position"), values=entrezgene, mart=ensembl)
genes$cisend <- genes$start_position+1000000
genes$cisstart <- genes$start_position-1000000
genes$Entrez_Gene_ID <- genes$entrezgene_id
x <- merge(merge,genes,by="Entrez_Gene_ID")
genes2 <- x[ -c(1:3,5:7,9) ]
genes2 <- genes2[, c(2,1,3,4,5)]
genes2 <- rename(genes2, c("ensembl_gene_id"="ensembl_gene_id","chr"="Chromosome","gene_start"="start_position","cis_start"="cisstart","cis_end"="cisend"))
write.csv(genes2,"Genesforloop.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(genes2,"Genesforloop.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl",host="uswest.ensembl.org")                                                                         
entrezgene = data.frame(merge$Entrez_Gene_ID)
entrezgene <- na.omit(entrezgene)
genes <- getBM(filters="entrezgene_id", attributes=c("ensembl_gene_id","entrezgene_id","start_position"), values=entrezgene, mart=ensembl)
genes$cisstart <- genes$start_position-1000000
genes$cisend <- genes$start_position+1000000
genes$Entrez_Gene_ID <- genes$entrezgene_id
x <- merge(merge,genes,by="Entrez_Gene_ID")
genes2 <- x[ -c(1:3,5:7,9) ]
genes2 <- genes2[, c(2,1,3,4,5)]
genes2 <- rename(genes2, c("ensembl_gene_id"="ensembl_gene_id","chr"="Chromosome","gene_start"="start_position","cis_start"="cisstart","cis_end"="cisend"))
write.csv(genes2,"Genesforloop.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(genes2,"Genesforloop.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)


#Can stop running here 


#Entrezgene IDs without results
genes$Entrez_Gene_ID <- genes$entrezgene_id
x <- merge(mergeall,genes,by="Entrez_Gene_ID")
found <- x$ID
found<- unique(found) #Some probe IDs have multiple entrez gene IDs
missing <- mergeall[ ! mergeall$ID %in% found, ]

#Unigenes
unigenes <- mergeall[mergeall$Unigene_ID != "", ]
unigenesid <- unigenes$Unigene_ID
#Tried: clusterProfiler, DAVID, biomart

#Other genes
unigenesprobeID <- unigenes$ID
missing <- missing[ ! missing$ID %in% unigenesprobeID, ]
#Still 87 missing
#Tried looking one up on DAVID and it says no conversion

#Looking into duplicates
#There are duplicates with the same EntrezID

#(on linux) get gencode ID from entrez id
library(data.table)
library(dplyr)
metadata <- fread("gencode.v38.metadata.EntrezGene")
genes <- fread("Merge.csv")
metadata$Entrez_Gene_ID <- metadata$V2
merge <- merge(genes,metadata,by="Entrez_Gene_ID")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("biomaRt")
library(biomaRt)
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")                                                                         
entrezgene = genes$Entrez_Gene_ID 
entrezgene <- na.omit(entrezgene)
genes <- getBM(filters="entrezgene_id", attributes=c("ensembl_gene_id","entrezgene_id","start_position","chromosome_name"), values=entrezgene, mart=ensembl,useCache=FALSE)                                                                                                                 
genes$cisstart <- genes$start_position-1000000
genes$cisend <- genes$start_position+1000000
genes$Entrez_Gene_ID <- genes$entrezgene_id
x <- merge(metadata,genes,by="Entrez_Gene_ID")
#We just want gencodeid	chr	gene_start	cis_start	cis_end
genes2 <- x[, -c(1,3,4,5) ]
genes2 <- genes2[, c(1,3,2,4,5)]
genes2 <- rename(genes2, c("id"="V1","chr"="chromosome_name","gene_start"="start_position","cis_start"="cisstart","cis_end"="cisend"))
#have transcript ID from metdata - need to match this to gene ID
#BiocManager::install("rtracklayer")
library(rtracklayer)
gtf <- rtracklayer::import('gencode.v38.annotation.gtf')
gtf_df=as.data.frame(gtf)
gtf_df$ensembl_gene_id <- substr(gtf_df$gene_id, 0, 15)
geneid<-gtf_df[,c(10,12,28)]
genes2$transcript_id <- genes2$id
x <- merge(geneid,genes2,by="transcript_id")
x$ensembl_gene_id <- x$gene_id
x <- x[, -c(1:2) ]
write.csv(x,"Genesforloop.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
