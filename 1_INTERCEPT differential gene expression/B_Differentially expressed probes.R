library(data.table)
Results <- read.csv("Results.csv")
n<-round(0.02*nrow(Results),digits=0)

res<-Results[1:n,]
res$FC<-2^res$logFC

# limma logFC is on the log2 scale: keep fold change > 1.3 in either direction
res<-res[abs(res$logFC)>log2(1.3),]

fwrite(res,"Differentially_expressed_probes.csv")
 