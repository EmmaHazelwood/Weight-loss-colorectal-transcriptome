library(data.table)
Results <- read.csv("Results.csv")
n<-round(0.02*nrow(Results),digits=0)

res<-Results[1:n,]
res$FC<-exp(res$logFC)

res1<-res[res$FC>1.4,]
res2<-res[res$FC<0.714,]
res<-rbind(res1,res2)

fwrite(res,"Genesforloop.csv")
