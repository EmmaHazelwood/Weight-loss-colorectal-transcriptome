
library(statmod)
library(limma)
library(illuminaio)

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
write.csv(res,"~/OneDrive - University of Cambridge/Documents/Bristol work/ABHD11 paper/Analysis/Results.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)

