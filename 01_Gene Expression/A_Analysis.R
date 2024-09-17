
library(statmod)
library(limma)

idatfiles <- dir(pattern="idat")
bgxfile <- dir(pattern="bgx") 
x <- read.idat(idatfiles, bgxfile,dateinfo = TRUE)
mapping <- read.delim("mapping.txt", header = TRUE, sep = "\t")
x$targets$Participant <- mapping$participant.id
x$targets$Time <- mapping$time.point
x$targets$BMI <- mapping$bmi
x$targets$Age <- mapping$age
x$targets$Sex <- mapping$sex
x$targets$waist <- mapping$waist
x$targets$hip <- mapping$hip
x$targets$igf1 <- mapping$igf1
x$targets$insulin <- mapping$insulin
x$targets$crp <- mapping$crp
x$targets$glucose <- mapping$glucose
x$targets$HOMA.IR <- mapping$HOMA.IR
x$targets$hdl <- mapping$hdl
x$targets$tc <- mapping$tc
x$targets$tg<-mapping$tg
x$targets$mTOR <- mapping$mTOR
x$targets$M30 <- mapping$M30
x$targets$Ki67 <- mapping$Ki67
x$targets$Apoptosis <- mapping$Apoptosis
x$targets$Date<-mapping$date

#Gives information on the probes
table(x$genes$Status)

#Expression value for each probe - intensities on log scale
options(digits=3)
head(x$E)

#p-values for how different each probe is (within each sample) from the negative control probes.
head(x$targets,1)

#Manually calculate and add in p values
x$other$Detection <- detectionPValues(x,negctrl="negative")

#Plot MA-plots
plotMD(x)

#Box plot of intensity per sample
P <- boxplot(log2(x$E),range=0,ylab="log2 intensity")

#Estimating the overall proportion of the regular probes that correspond to expressed transcript
pe <- propexpr(x)
dim(pe) <- c(1,24)
dimnames(pe) <- list()
pe

#Background correction and normalize:
#can't use normalizeWithinArrays() as don't have R and G components. backgroundCorrect(method="subtract") ran without errors but did not alter data at all. However, eith method="normexp" this does work.
y <- neqc(x)
#neqc() performs normexp background correction using negative controls, then quantile normalizes and finally log2 transforms. It also automatically removes the control probes, leaving only the regular probes in y
dim(x)
dim(y)
#Gets rid of 887

#Compare densities before and after normalization
plotDensities(x,legend=FALSE)
plotDensities(y,legend=FALSE)

#Filter out probes that are not expressed. Keep probes that are expressed in at least three arrays according to a detection p-value of 5%
expressed <- rowSums(y$other$Detection < 0.05) >=3
y <- y[expressed,]
dim(y) #Excludes 20,770 genes, so 26553 left

#Multi-dimensional scaling plot - should hopefully group individuals together.
plotMDS(y,labels=y$targets$Participant)

#Export raw data for PCA later
write.csv(y$E,"Data for PCA.csv",row.names = FALSE,quote=FALSE)

#Differential expression between before and after weight loss
#Moderated paired t-test - simple, not adjusting for anything
Time <- factor(y$targets$Time,levels = c("Baseline","Followup"))
Participant <- factor(y$targets$Participant)
design <- model.matrix(~Participant+Time)
fit <- lmFit(y,design)
fit <- eBayes(fit) #automatically adjusts B statistic by assuming that 1% of the genes are differentially expressed - can be altered
topTable(fit,coef="TimeFollowup")

#Weighting by array quality
arrayw <- arrayWeights(y)
barplot(arrayw,xlab="Array",ylab="Weight",col="white",las=2)
fitw <- lmFit(y,design,weights=arrayw)
fitw <- eBayes(fitw) 
topTable(fitw,coef="TimeFollowup")

#Adjusting for within patient correlation
#Within-patient correlation
design <- model.matrix(~0+Time)
colnames(design)<- levels(Time)
dupcor <- duplicateCorrelation(y,design,block=y$targets$Participant)
dupcor$consensus.correlation
#within donor is 0.0865 - small, but positive, as expected

design <- model.matrix(~Participant+Time, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
topTable(fit, coef="TimeFollowup")

#Adjusting for other factors
Age <- as.numeric(y$targets$Age)
Sex <- factor(y$targets$Sex,levels = c("F","M"))

#Including Age and Sex
design <- model.matrix(~Participant+Time+Age+Sex, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
topTable(fit, coef="TimeFollowup")
#Don't want to adjust for BMI as this will remove effect of weight loss

#Save as table
res <- topTable(fit, coef="TimeFollowup",number=26553)
write.csv(res,"Results.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)

#Effect of date
Date <- as.numeric(y$targets$date)
df<-y$targets
Participant.tbl<-dplyr::select(df,Participant,Time,Date)
Participant.tbl<-Participant.tbl[order(Participant.tbl$Time),]
Participant.tbl<-Participant.tbl[match(unique(Participant.tbl$Participant),Participant.tbl$Participant),]
df<-merge(df,Participant.tbl,by=c("Participant","Time"),all.x = T)
df<-df[match(y$targets$Participant,df$Participant),]
y$targets$Date<-df$Date.x
Date<-as.POSIXct(y$targets$Date,format="%d/%m/%Y")

Date<-as.numeric(Date)


design <- model.matrix(~Participant+Time+Age+Sex+Date, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
resDate <- topTable(fit, coef="TimeFollowup",number=26553)

comp<-merge(res,resDate,by="Probe_Id")
plot(comp$logFC.x,comp$logFC.y)
plot(comp$P.Value.x,comp$P.Value.y)


plotMDS(y,labels=y$targets$Date)

#Trying to explain the changes
BMI <- as.numeric(y$targets$BMI)
WC <- as.numeric(y$targets$waist)
HC <- as.numeric(y$targets$hip)
WHR <- as.numeric(WC/HC)
IGF1 <- as.numeric(y$targets$igf1)
Insulin <- as.numeric(y$targets$insulin)
CRP <- as.numeric(y$targets$crp)
Glu <- as.numeric(y$targets$glucose)
HOMAIR <- as.numeric(y$targets$HOMA.IR)
HDL <- as.numeric(y$targets$hdl)
TC <- as.numeric(y$targets$tc)
TG <- as.numeric(y$targets$tg)
mTOR <- as.numeric(y$targets$mTOR)
M30 <- as.numeric(y$targets$M30)
Ki67 <- as.numeric(y$targets$Ki67)
Apo <- as.numeric(y$targets$Apoptosis)

#BMI
design <- model.matrix(~Participant+Time+Age+Sex+BMI, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res2 <- topTable(fit, coef="TimeFollowup",number=26553)

#WC
design <- model.matrix(~Participant+Time+Age+Sex+WC, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res3 <- topTable(fit, coef="TimeFollowup",number=26553)

#HC
design <- model.matrix(~Participant+Time+Age+Sex+HC, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res4 <- topTable(fit, coef="TimeFollowup",number=26553)

#WHR
design <- model.matrix(~Participant+Time+Age+Sex+WHR, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res5 <- topTable(fit, coef="TimeFollowup",number=26553)

#Apo
design <- model.matrix(~Participant+Time+Age+Sex+Apo, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res16 <- topTable(fit, coef="TimeFollowup",number=26553)

#Removing missing data
idatfiles <- dir(pattern="idat")
idatfiles <- idatfiles[idatfiles!="3998479018_J_Grn.idat"]
idatfiles <- idatfiles[idatfiles!="3998479018_L_Grn.idat"]
bgxfile <- dir(pattern="bgx") 
x <- read.idat(idatfiles, bgxfile,dateinfo = TRUE)
mapping <- read.delim("mapping.txt", header = TRUE, sep = "\t")
mapping <- mapping[mapping$participant.id!="INT-22",]
x$targets$Participant <- mapping$participant.id
x$targets$Time <- mapping$time.point
x$targets$BMI <- mapping$bmi
x$targets$Age <- mapping$age
x$targets$Sex <- mapping$sex
x$targets$waist <- mapping$waist
x$targets$hip <- mapping$hip
x$targets$igf1 <- mapping$igf1
x$targets$insulin <- mapping$insulin
x$targets$crp <- mapping$crp
x$targets$glucose <- mapping$glucose
x$targets$HOMA.IR <- mapping$HOMA.IR
x$targets$hdl <- mapping$hdl
x$targets$tc <- mapping$tc
x$targets$tg<-mapping$tg
x$targets$mTOR <- mapping$mTOR
x$targets$M30 <- mapping$M30
x$targets$Ki67 <- mapping$Ki67
x$targets$Apoptosis <- mapping$Apoptosis

x$other$Detection <- detectionPValues(x,negctrl="negative")

y <- neqc(x)

expressed <- rowSums(y$other$Detection < 0.05) >=3
y <- y[expressed,]

Time <- factor(y$targets$Time,levels = c("Baseline","Followup"))
Participant <- factor(y$targets$Participant)
arrayw <- arrayWeights(y)
design <- model.matrix(~0+Time)
colnames(design)<- levels(Time)
dupcor <- duplicateCorrelation(y,design,block=y$targets$Participant)
dupcor$consensus.correlation
Age <- as.numeric(y$targets$Age)
Sex <- factor(y$targets$Sex,levels = c("F","M"))

BMI <- as.numeric(y$targets$BMI)
WC <- as.numeric(y$targets$waist)
HC <- as.numeric(y$targets$hip)
WHR <- as.numeric(WC/HC)
IGF1 <- as.numeric(y$targets$igf1)
Insulin <- as.numeric(y$targets$insulin)
CRP <- as.numeric(y$targets$crp)
Glu <- as.numeric(y$targets$glucose)
HOMAIR <- as.numeric(y$targets$HOMA.IR)
HDL <- as.numeric(y$targets$hdl)
TC <- as.numeric(y$targets$tc)
TG <- as.numeric(y$targets$tg)
mTOR <- as.numeric(y$targets$mTOR)
M30 <- as.numeric(y$targets$M30)
Ki67 <- as.numeric(y$targets$Ki67)
Apo <- as.numeric(y$targets$Apoptosis)

#IGF1
design <- model.matrix(~Participant+Time+Age+Sex+IGF1, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res6 <- topTable(fit, coef="TimeFollowup",number=26553)

#Insulin
design <- model.matrix(~Participant+Time+Age+Sex+Insulin, correlation=dupcor$consensus.correlation)


fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res7 <- topTable(fit, coef="TimeFollowup",number=26553)

#CRP
design <- model.matrix(~Participant+Time+Age+Sex+CRP, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res8 <- topTable(fit, coef="TimeFollowup",number=26553)

#Glu
design <- model.matrix(~Participant+Time+Age+Sex+Glu, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res9 <- topTable(fit, coef="TimeFollowup",number=26553)

#HOMAIR
design <- model.matrix(~Participant+Time+Age+Sex+HOMAIR, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res10 <- topTable(fit, coef="TimeFollowup",number=26553)

#HDL
design <- model.matrix(~Participant+Time+Age+Sex+HDL, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res11 <- topTable(fit, coef="TimeFollowup",number=26553)

#TC
design <- model.matrix(~Participant+Time+Age+Sex+TC, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res12 <- topTable(fit, coef="TimeFollowup",number=26553)

#TG
design <- model.matrix(~Participant+Time+Age+Sex+TG, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res17 <- topTable(fit, coef="TimeFollowup",number=26553)

#Remove missing data again
idatfiles <- dir(pattern="idat")
idatfiles <- idatfiles[idatfiles!="3998479018_J_Grn.idat"]
idatfiles <- idatfiles[idatfiles!="3998479018_L_Grn.idat"]
idatfiles <- idatfiles[idatfiles!="3998479018_A_Grn.idat"]
idatfiles <- idatfiles[idatfiles!="3998479018_C_Grn.idat"]
idatfiles <- idatfiles[idatfiles!="3998479022_E_Grn.idat"]
idatfiles <- idatfiles[idatfiles!="3998479022_I_Grn.idat"]
bgxfile <- dir(pattern="bgx") 
x <- read.idat(idatfiles, bgxfile,dateinfo = TRUE)
mapping <- read.delim("mapping.txt", header = TRUE, sep = "\t")
mapping <- mapping[mapping$participant.id!="INT-22",]
mapping <- mapping[mapping$participant.id!="INT-19",]
x$targets$Participant <- mapping$participant.id
x$targets$Time <- mapping$time.point
x$targets$BMI <- mapping$bmi
x$targets$Age <- mapping$age
x$targets$Sex <- mapping$sex
x$targets$waist <- mapping$waist
x$targets$hip <- mapping$hip
x$targets$igf1 <- mapping$igf1
x$targets$insulin <- mapping$insulin
x$targets$crp <- mapping$crp
x$targets$glucose <- mapping$glucose
x$targets$HOMA.IR <- mapping$HOMA.IR
x$targets$hdl <- mapping$hdl
x$targets$tc <- mapping$tc
x$targets$mTOR <- mapping$mTOR
x$targets$M30 <- mapping$M30
x$targets$Ki67 <- mapping$Ki67
x$targets$Apoptosis <- mapping$Apoptosis

x$other$Detection <- detectionPValues(x,negctrl="negative")

y <- neqc(x)

expressed <- rowSums(y$other$Detection < 0.05) >=3
y <- y[expressed,]

Time <- factor(y$targets$Time,levels = c("Baseline","Followup"))
Participant <- factor(y$targets$Participant)
arrayw <- arrayWeights(y)
design <- model.matrix(~0+Time)
colnames(design)<- levels(Time)
dupcor <- duplicateCorrelation(y,design,block=y$targets$Participant)
dupcor$consensus.correlation
Age <- as.numeric(y$targets$Age)
Sex <- factor(y$targets$Sex,levels = c("F","M"))

BMI <- as.numeric(y$targets$BMI)
WC <- as.numeric(y$targets$waist)
HC <- as.numeric(y$targets$hip)
WHR <- as.numeric(WC/HC)
IGF1 <- as.numeric(y$targets$igf1)
Insulin <- as.numeric(y$targets$insulin)
CRP <- as.numeric(y$targets$crp)
Glu <- as.numeric(y$targets$glucose)
HOMAIR <- as.numeric(y$targets$HOMA.IR)
HDL <- as.numeric(y$targets$hdl)
TC <- as.numeric(y$targets$tc)
mTOR <- as.numeric(y$targets$mTOR)
M30 <- as.numeric(y$targets$M30)
Ki67 <- as.numeric(y$targets$Ki67)
Apo <- as.numeric(y$targets$Apoptosis)

#Ki67
design <- model.matrix(~Participant+Time+Age+Sex+Ki67, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res15 <- topTable(fit, coef="TimeFollowup",number=26553)

#Remove missing data again again
idatfiles <- dir(pattern="idat")
idatfiles <- idatfiles[idatfiles!="3998479018_J_Grn.idat"]
idatfiles <- idatfiles[idatfiles!="3998479018_L_Grn.idat"]
idatfiles <- idatfiles[idatfiles!="3998479018_A_Grn.idat"]
idatfiles <- idatfiles[idatfiles!="3998479018_C_Grn.idat"]
idatfiles <- idatfiles[idatfiles!="3998479022_E_Grn.idat"]
idatfiles <- idatfiles[idatfiles!="3998479022_I_Grn.idat"]
idatfiles <- idatfiles[idatfiles!="3998479018_K_Grn.idat"]
idatfiles <- idatfiles[idatfiles!="3998479018_D_Grn.idat"]
idatfiles <- idatfiles[idatfiles!="3998479022_G_Grn.idat"]
bgxfile <- dir(pattern="bgx") 
x <- read.idat(idatfiles, bgxfile,dateinfo = TRUE)
mapping <- read.delim("mapping.txt", header = TRUE, sep = "\t")
mapping <- mapping[mapping$IDAT.file!="3998479018_J",]
mapping <- mapping[mapping$IDAT.file!="3998479018_L",]
mapping <- mapping[mapping$IDAT.file!="3998479018_A",]
mapping <- mapping[mapping$IDAT.file!="3998479018_C",]
mapping <- mapping[mapping$IDAT.file!="3998479022_E",]
mapping <- mapping[mapping$IDAT.file!="3998479022_I",]
mapping <- mapping[mapping$IDAT.file!="3998479018_K",]
mapping <- mapping[mapping$IDAT.file!="3998479018_D",]
mapping <- mapping[mapping$IDAT.file!="3998479022_G",]
x$targets$Participant <- mapping$participant.id
x$targets$Time <- mapping$time.point
x$targets$BMI <- mapping$bmi
x$targets$Age <- mapping$age
x$targets$Sex <- mapping$sex
x$targets$waist <- mapping$waist
x$targets$hip <- mapping$hip
x$targets$igf1 <- mapping$igf1
x$targets$insulin <- mapping$insulin
x$targets$crp <- mapping$crp
x$targets$glucose <- mapping$glucose
x$targets$HOMA.IR <- mapping$HOMA.IR
x$targets$hdl <- mapping$hdl
x$targets$tc <- mapping$tc
x$targets$mTOR <- mapping$mTOR
x$targets$M30 <- mapping$M30
x$targets$Ki67 <- mapping$Ki67
x$targets$Apoptosis <- mapping$Apoptosis

x$other$Detection <- detectionPValues(x,negctrl="negative")

y <- neqc(x)

expressed <- rowSums(y$other$Detection < 0.05) >=3
y <- y[expressed,]

Time <- factor(y$targets$Time,levels = c("Baseline","Followup"))
Participant <- factor(y$targets$Participant)
arrayw <- arrayWeights(y)
design <- model.matrix(~0+Time)
colnames(design)<- levels(Time)
dupcor <- duplicateCorrelation(y,design,block=y$targets$Participant)
dupcor$consensus.correlation
Age <- as.numeric(y$targets$Age)
Sex <- factor(y$targets$Sex,levels = c("F","M"))

BMI <- as.numeric(y$targets$BMI)
WC <- as.numeric(y$targets$waist)
HC <- as.numeric(y$targets$hip)
WHR <- as.numeric(WC/HC)
IGF1 <- as.numeric(y$targets$igf1)
Insulin <- as.numeric(y$targets$insulin)
CRP <- as.numeric(y$targets$crp)
Glu <- as.numeric(y$targets$glucose)
HOMAIR <- as.numeric(y$targets$HOMA.IR)
HDL <- as.numeric(y$targets$hdl)
TC <- as.numeric(y$targets$tc)
mTOR <- as.numeric(y$targets$mTOR)
M30 <- as.numeric(y$targets$M30)
Ki67 <- as.numeric(y$targets$Ki67)
Apo <- as.numeric(y$targets$Apoptosis)

#mTOR
design <- model.matrix(~Participant+Time+Age+Sex+mTOR, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res13 <- topTable(fit, coef="TimeFollowup",number=26553)

#m30
design <- model.matrix(~Participant+Time+Age+Sex+M30, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
res14 <- topTable(fit, coef="TimeFollowup",number=26553)

#Do any other variables explain the initial results?
resa <- res[1:26553,]
genelist <- resa$Probe_Id
a<-res2
res2 <- a[a$Probe_Id %in% genelist ,]
a<-res3
res3 <- a[a$Probe_Id %in% genelist ,]
a<-res4
res4 <- a[a$Probe_Id %in% genelist ,]
a<-res5
res5 <- a[a$Probe_Id %in% genelist ,]
a<-res6
res6 <- a[a$Probe_Id %in% genelist ,]
a<-res7
res7 <- a[a$Probe_Id %in% genelist ,]
a<-res8
res8 <- a[a$Probe_Id %in% genelist ,]
a<-res9
res9 <- a[a$Probe_Id %in% genelist ,]
a<-res10
res10 <- a[a$Probe_Id %in% genelist ,]
a<-res11
res11 <- a[a$Probe_Id %in% genelist ,]
a<-res12
res12 <- a[a$Probe_Id %in% genelist ,]
a<-res13
res13 <- a[a$Probe_Id %in% genelist ,]
a<-res14
res14 <- a[a$Probe_Id %in% genelist ,]
a<-res15
res15 <- a[a$Probe_Id %in% genelist ,]
a<-res16
res16 <- a[a$Probe_Id %in% genelist ,]
a<-res17
res17 <- a[a$Probe_Id %in% genelist ,]

write.csv(res2,"Results BMI.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res3,"Results WC.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res4,"Results HC.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res5,"Results WHR.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res6,"Results IGF1.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res7,"Results Insulin.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res8,"Results CRP.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res9,"Results Glucose.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res10,"Results HOMAIR.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res11,"Results HDL.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res12,"Results TC.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res13,"Results Ki67.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res14,"Results mTOR.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res15,"Results m30.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res16,"Results Apo.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)
write.csv(res17,"Results TG.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)

#Effect of date