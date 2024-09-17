library(data.table)
library(dplyr)

x<-dir(pattern="Overall")

results=data.frame(matrix(ncol=10,nrow=0))
for (i in 1:length(x)){
  print(i)
  y<-dir(pattern=x[i])
  x2<-read.table(y,sep=",")
  results=rbind(results,x2)
}

write.table(results, file = "GTExMetaAnalysisResults.txt", sep = "\t",
            row.names = F)
