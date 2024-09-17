library(sparklyr)
spark_install(version="3.0.0")
sc<-spark_connect(master="local")

chr1 <- spark_read_parquet(sc, "chr1", "Colon_Transverse.v8.EUR.allpairs.chr1.parquet")
spark_write_csv(chr1, "Colon_Transverse.v8.EUR.allpairs.chr1.csv")

chr2 <- spark_read_parquet(sc, "chr2", "Colon_Transverse.v8.EUR.allpairs.chr2.parquet")
spark_write_csv(chr2, "Colon_Transverse.v8.EUR.allpairs.chr2.csv")

chr3 <- spark_read_parquet(sc, "chr3", "Colon_Transverse.v8.EUR.allpairs.chr3.parquet")
spark_write_csv(chr3, "Colon_Transverse.v8.EUR.allpairs.chr3.csv")

chr4 <- spark_read_parquet(sc, "chr4", "Colon_Transverse.v8.EUR.allpairs.chr4.parquet")
spark_write_csv(chr4, "Colon_Transverse.v8.EUR.allpairs.chr4.csv")

chr5 <- spark_read_parquet(sc, "chr5", "Colon_Transverse.v8.EUR.allpairs.chr5.parquet")
spark_write_csv(chr5, "Colon_Transverse.v8.EUR.allpairs.chr5.csv")

chr6 <- spark_read_parquet(sc, "chr6", "Colon_Transverse.v8.EUR.allpairs.chr6.parquet")
spark_write_csv(chr6, "Colon_Transverse.v8.EUR.allpairs.chr6.csv")

chr7 <- spark_read_parquet(sc, "chr7", "Colon_Transverse.v8.EUR.allpairs.chr7.parquet")
spark_write_csv(chr7, "Colon_Transverse.v8.EUR.allpairs.chr7.csv")

chr8 <- spark_read_parquet(sc, "chr8", "Colon_Transverse.v8.EUR.allpairs.chr8.parquet")
spark_write_csv(chr8, "Colon_Transverse.v8.EUR.allpairs.chr8.csv")

chr9 <- spark_read_parquet(sc, "chr9", "Colon_Transverse.v8.EUR.allpairs.chr9.parquet")
spark_write_csv(chr9, "Colon_Transverse.v8.EUR.allpairs.chr9.csv")

chr10 <- spark_read_parquet(sc, "chr10", "Colon_Transverse.v8.EUR.allpairs.chr10.parquet")
spark_write_csv(chr10, "Colon_Transverse.v8.EUR.allpairs.chr10.csv")

chr11 <- spark_read_parquet(sc, "chr11", "Colon_Transverse.v8.EUR.allpairs.chr11.parquet")
spark_write_csv(chr11, "Colon_Transverse.v8.EUR.allpairs.chr11.csv")

chr12 <- spark_read_parquet(sc, "chr12", "Colon_Transverse.v8.EUR.allpairs.chr12.parquet")
spark_write_csv(chr12, "Colon_Transverse.v8.EUR.allpairs.chr12.csv")

chr13 <- spark_read_parquet(sc, "chr13", "Colon_Transverse.v8.EUR.allpairs.chr13.parquet")
spark_write_csv(chr13, "Colon_Transverse.v8.EUR.allpairs.chr13.csv")

chr14 <- spark_read_parquet(sc, "chr14", "Colon_Transverse.v8.EUR.allpairs.chr14.parquet")
spark_write_csv(chr14, "Colon_Transverse.v8.EUR.allpairs.chr14.csv")

chr15 <- spark_read_parquet(sc, "chr15", "Colon_Transverse.v8.EUR.allpairs.chr15.parquet")
spark_write_csv(chr15, "Colon_Transverse.v8.EUR.allpairs.chr15.csv")

chr16 <- spark_read_parquet(sc, "chr16", "Colon_Transverse.v8.EUR.allpairs.chr16.parquet")
spark_write_csv(chr16, "Colon_Transverse.v8.EUR.allpairs.chr16.csv")

chr17 <- spark_read_parquet(sc, "chr17", "Colon_Transverse.v8.EUR.allpairs.chr17.parquet")
spark_write_csv(chr17, "Colon_Transverse.v8.EUR.allpairs.chr17.csv")

chr18 <- spark_read_parquet(sc, "chr18", "Colon_Transverse.v8.EUR.allpairs.chr18.parquet")
spark_write_csv(chr18, "Colon_Transverse.v8.EUR.allpairs.chr18.csv")

chr19 <- spark_read_parquet(sc, "chr19", "Colon_Transverse.v8.EUR.allpairs.chr19.parquet")
spark_write_csv(chr19, "Colon_Transverse.v8.EUR.allpairs.chr19.csv")

chr20 <- spark_read_parquet(sc, "chr20", "Colon_Transverse.v8.EUR.allpairs.chr20.parquet")
spark_write_csv(chr20, "Colon_Transverse.v8.EUR.allpairs.chr20.csv")

chr21 <- spark_read_parquet(sc, "chr21", "Colon_Transverse.v8.EUR.allpairs.chr21.parquet")
spark_write_csv(chr21, "Colon_Transverse.v8.EUR.allpairs.chr21.csv")

chr22 <- spark_read_parquet(sc, "chr22", "Colon_Transverse.v8.EUR.allpairs.chr22.parquet")
spark_write_csv(chr22, "Colon_Transverse.v8.EUR.allpairs.chr22.csv")

chrX <- spark_read_parquet(sc, "chrX", "Colon_Transverse.v8.EUR.allpairs.chrX.parquet")
spark_write_csv(chrX, "Colon_Transverse.v8.EUR.allpairs.chrX.csv")