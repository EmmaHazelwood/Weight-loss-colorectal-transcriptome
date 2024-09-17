#install.packages("data.table")
#1

library(data.table)

chr1 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr1.csv",pattern = "part-00007",full.names = TRUE)
chr1 <- fread(chr1)
chr1 <- data.frame(chr1)

chr2 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr2.csv",pattern = "part-00007",full.names = TRUE)
chr2 <- fread(chr2)
chr2 <- data.frame(chr2)

chr3 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr3.csv",pattern = "part-00007",full.names = TRUE)
chr3 <- fread(chr3)
chr3 <- data.frame(chr3)

chr4 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr4.csv",pattern = "part-00007",full.names = TRUE)
chr4 <- fread(chr4)
chr4 <- data.frame(chr4)

chr5 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr5.csv",pattern = "part-00007",full.names = TRUE)
chr5 <- fread(chr5)
chr5 <- data.frame(chr5)

chr6 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr6.csv",pattern = "part-00007",full.names = TRUE)
chr6 <- fread(chr6)
chr6 <- data.frame(chr6)

chr7 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr7.csv",pattern = "part-00007",full.names = TRUE)
chr7 <- fread(chr7)
chr7 <- data.frame(chr7)

chr8 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr8.csv",pattern = "part-00007",full.names = TRUE)
chr8 <- fread(chr8)
chr8 <- data.frame(chr8)

chr9 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr9.csv",pattern = "part-00007",full.names = TRUE)
chr9 <- fread(chr9)
chr9 <- data.frame(chr9)

chr10 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr10.csv",pattern = "part-00007",full.names = TRUE)
chr10 <- fread(chr10)
chr10 <- data.frame(chr10)

chr11 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr11.csv",pattern = "part-00007",full.names = TRUE)
chr11 <- fread(chr11)
chr11 <- data.frame(chr11)

chr12 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr12.csv",pattern = "part-00007",full.names = TRUE)
chr12 <- fread(chr12)
chr12 <- data.frame(chr12)

chr13 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr13.csv",pattern = "part-00007",full.names = TRUE)
chr13 <- fread(chr13)
chr13 <- data.frame(chr13)

chr14 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr14.csv",pattern = "part-00007",full.names = TRUE)
chr14 <- fread(chr14)
chr14 <- data.frame(chr14)

chr15 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr15.csv",pattern = "part-00007",full.names = TRUE)
chr15 <- fread(chr15)
chr15 <- data.frame(chr15)

chr16 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr16.csv",pattern = "part-00007",full.names = TRUE)
chr16 <- fread(chr16)
chr16 <- data.frame(chr16)

chr17 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr17.csv",pattern = "part-00007",full.names = TRUE)
chr17 <- fread(chr17)
chr17 <- data.frame(chr17)

chr18 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr18.csv",pattern = "part-00007",full.names = TRUE)
chr18 <- fread(chr18)
chr18 <- data.frame(chr18)

chr19 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr19.csv",pattern = "part-00007",full.names = TRUE)
chr19 <- fread(chr19)
chr19 <- data.frame(chr19)

chr20 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr20.csv",pattern = "part-00007",full.names = TRUE)
chr20 <- fread(chr20)
chr20 <- data.frame(chr20)

chr21 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr21.csv",pattern = "part-00007",full.names = TRUE)
chr21 <- fread(chr21)
chr21 <- data.frame(chr21)

chr22 <- dir(path="Colon_Transverse.v8.EUR.allpairs.chr22.csv",pattern = "part-00007",full.names = TRUE)
chr22 <- fread(chr22)
chr22 <- data.frame(chr22)

chrX <- dir(path="Colon_Transverse.v8.EUR.allpairs.chrX.csv",pattern = "part-00007",full.names = TRUE)
chrX <- fread(chrX)
chrX <- data.frame(chrX)

chr <- rbind(chr1,chr2)
chr <- rbind(chr,chr3)
chr <- rbind(chr,chr3)
chr <- rbind(chr,chr4)
chr <- rbind(chr,chr5)
chr <- rbind(chr,chr6)
chr <- rbind(chr,chr7)
chr <- rbind(chr,chr8)
chr <- rbind(chr,chr9)
chr <- rbind(chr,chr10)
chr <- rbind(chr,chr11)
chr <- rbind(chr,chr12)
chr <- rbind(chr,chr13)
chr <- rbind(chr,chr14)
chr <- rbind(chr,chr15)
chr <- rbind(chr,chr16)
chr <- rbind(chr,chr17)
chr <- rbind(chr,chr18)
chr <- rbind(chr,chr19)
chr <- rbind(chr,chr20)
chr <- rbind(chr,chr21)
chr <- rbind(chr,chr22)
chr <- rbind(chr,chrX)

write.csv(chr,"ColonTransverse.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
