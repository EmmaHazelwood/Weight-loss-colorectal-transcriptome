
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

#Sample information
head(x$targets,1)

#Detection p-values relative to the negative control probes
x$other$Detection <- detectionPValues(x,negctrl="negative")

#Plot MA-plots
plotMD(x)

#Box plot of intensity per sample
P <- boxplot(log2(x$E),range=0,ylab="log2 intensity")

#Estimated proportion of probes corresponding to expressed transcripts
pe <- propexpr(x)
dim(pe) <- c(1,24)
dimnames(pe) <- list()
pe

#normexp background correction using negative controls, quantile normalisation and log2 transformation; control probes are removed
y <- neqc(x)
dim(x)
dim(y)

#Densities before and after normalisation
plotDensities(x,legend=FALSE)
plotDensities(y,legend=FALSE)

#Keep probes expressed (detection P < 0.05) in at least three arrays
expressed <- rowSums(y$other$Detection < 0.05) >=3
y <- y[expressed,]
dim(y)

#Multi-dimensional scaling plot
plotMDS(y,labels=y$targets$Participant)

#Normalised expression for PCA
write.csv(y$E,"Data for PCA.csv",row.names = FALSE,quote=FALSE)

#Differential expression between before and after weight loss
#Moderated paired t-test, unadjusted
Time <- factor(y$targets$Time,levels = c("Baseline","Followup"))
Participant <- factor(y$targets$Participant)
design <- model.matrix(~Participant+Time)
fit <- lmFit(y,design)
fit <- eBayes(fit)
topTable(fit,coef="TimeFollowup")

#Weighting by array quality
arrayw <- arrayWeights(y)
barplot(arrayw,xlab="Array",ylab="Weight",col="white",las=2)
fitw <- lmFit(y,design,weights=arrayw)
fitw <- eBayes(fitw) 
topTable(fitw,coef="TimeFollowup")

#Within-participant correlation
design <- model.matrix(~0+Time)
colnames(design)<- levels(Time)
dupcor <- duplicateCorrelation(y,design,block=y$targets$Participant)
dupcor$consensus.correlation

design <- model.matrix(~Participant+Time, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
topTable(fit, coef="TimeFollowup")

#Additional covariates
Age <- as.numeric(y$targets$Age)
Sex <- factor(y$targets$Sex,levels = c("F","M"))

#Including Age and Sex
design <- model.matrix(~Participant+Time+Age+Sex, correlation=dupcor$consensus.correlation)
fit <- lmFit(y,design,weights=arrayw,block=y$targets$Participant,correlation=dupcor$consensus.correlation)
fit <-eBayes(fit)
topTable(fit, coef="TimeFollowup")
#Not adjusted for BMI, as weight loss is the exposure

#Save as table
res <- topTable(fit, coef="TimeFollowup",number=26553)
write.csv(res,"Results.csv", row.names = FALSE,col.names = TRUE,quote = FALSE)

