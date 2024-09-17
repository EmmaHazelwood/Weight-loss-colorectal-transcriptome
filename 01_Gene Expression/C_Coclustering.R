library(statmod)
library(limma)
library(dplyr)
library(data.table)
library(ComplexHeatmap)

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

#Manually calculate and add in p values
x$other$Detection <- detectionPValues(x,negctrl="negative")

#Background correction and normalize:
#can't use normalizeWithinArrays() as don't have R and G components. backgroundCorrect(method="subtract") ran without errors but did not alter data at all. However, eith method="normexp" this does work.
y <- neqc(x)

#Filter out probes that are not expressed. Keep probes that are expressed in at least three arrays according to a detection p-value of 5%
expressed <- rowSums(y$other$Detection < 0.05) >=3
y <- y[expressed,]


#Standardise - divide all continuous values by SD
genes<-t(apply(y$E,1, function(x){x/sd(x)}))
genes<-t(apply(y$E,1, function(x){x/sd(x)}))
bmi<-as.data.frame(y$targets$BMI/(sd(y$targets$BMI)))
bmi<-cbind(bmi,y$targets$IDATfile)
bmi$IDATfile<-substr(bmi$`y$targets$IDATfile`,1,16)
ins<-as.data.frame(as.numeric(y$targets$insulin)/(sd(na.omit(as.numeric(y$targets$insulin)))))
ins<-cbind(ins,y$targets$IDATfile)
ins$IDATfile<-substr(ins$`y$targets$IDATfile`,1,16)
y$targets$WHR<-y$targets$waist/y$targets$hip
WHR<-as.data.frame(as.numeric(y$targets$WHR)/(sd(na.omit(as.numeric(y$targets$WHR)))))
WHR<-cbind(WHR,y$targets$IDATfile)
WHR$IDATfile<-substr(WHR$`y$targets$IDATfile`,1,16)
igf1<-as.data.frame(as.numeric(y$targets$igf1)/(sd(na.omit(as.numeric(y$targets$igf1)))))
igf1<-cbind(igf1,y$targets$IDATfile)
igf1$IDATfile<-substr(igf1$`y$targets$IDATfile`,1,16)
crp<-as.data.frame(as.numeric(y$targets$crp)/(sd(na.omit(as.numeric(y$targets$crp)))))
crp<-cbind(crp,y$targets$IDATfile)
crp$IDATfile<-substr(crp$`y$targets$IDATfile`,1,16)
glucose<-as.data.frame(as.numeric(y$targets$glucose)/(sd(na.omit(as.numeric(y$targets$glucose)))))
glucose<-cbind(glucose,y$targets$IDATfile)
glucose$IDATfile<-substr(glucose$`y$targets$IDATfile`,1,16)
HOMA.IR<-as.data.frame(as.numeric(y$targets$HOMA.IR)/(sd(na.omit(as.numeric(y$targets$HOMA.IR)))))
HOMA.IR<-cbind(HOMA.IR,y$targets$IDATfile)
HOMA.IR$IDATfile<-substr(HOMA.IR$`y$targets$IDATfile`,1,16)
hdl<-as.data.frame(as.numeric(y$targets$hdl)/(sd(na.omit(as.numeric(y$targets$hdl)))))
hdl<-cbind(hdl,y$targets$IDATfile)
hdl$IDATfile<-substr(hdl$`y$targets$IDATfile`,1,16)
tc<-as.data.frame(as.numeric(y$targets$tc)/(sd(na.omit(as.numeric(y$targets$tc)))))
tc<-cbind(tc,y$targets$IDATfile)
tc$IDATfile<-substr(tc$`y$targets$IDATfile`,1,16)
tg<-as.data.frame(as.numeric(y$targets$tg)/(sd(na.omit(as.numeric(y$targets$tg)))))
tg<-cbind(tg,y$targets$IDATfile)
tg$IDATfile<-substr(tg$`y$targets$IDATfile`,1,16)
mTOR<-as.data.frame(as.numeric(y$targets$mTOR)/(sd(na.omit(as.numeric(y$targets$mTOR)))))
mTOR<-cbind(mTOR,y$targets$IDATfile)
mTOR$IDATfile<-substr(mTOR$`y$targets$IDATfile`,1,16)
M30<-as.data.frame(as.numeric(y$targets$M30)/(sd(na.omit(as.numeric(y$targets$M30)))))
M30<-cbind(M30,y$targets$IDATfile)
M30$IDATfile<-substr(M30$`y$targets$IDATfile`,1,16)
Ki67<-as.data.frame(as.numeric(y$targets$Ki67)/(sd(na.omit(as.numeric(y$targets$Ki67)))))
Ki67<-cbind(Ki67,y$targets$IDATfile)
Ki67$IDATfile<-substr(Ki67$`y$targets$IDATfile`,1,16)
Apoptosis<-as.data.frame(as.numeric(y$targets$Apoptosis)/(sd(na.omit(as.numeric(y$targets$Apoptosis)))))
Apoptosis<-cbind(Apoptosis,y$targets$IDATfile)
Apoptosis$IDATfile<-substr(Apoptosis$`y$targets$IDATfile`,1,16)

df<-data.frame(matrix(ncol=26568,nrow=8))
colnames(df)<-"Participant"
df$Participant<-unique(y$targets$Participant)
colnames(df)[2:26554]<-rownames(y$E)
colnames(df)[26555]<-"BMI"
colnames(df)[26556]<-"Insulin"
colnames(df)[26557]<-"WHR"
colnames(df)[26558]<-"igf1"
colnames(df)[26559]<-"crp"
colnames(df)[26560]<-"glucose"
colnames(df)[26561]<-"HOMA.IR"
colnames(df)[26562]<-"hdl"
colnames(df)[26563]<-"tc"
colnames(df)[26564]<-"tg"
colnames(df)[26565]<-"mTOR"
colnames(df)[26566]<-"M30"
colnames(df)[26567]<-"Ki67"
colnames(df)[26568]<-"Apoptosis"

for (i in unique(y$targets$Participant)){
samples<-as.data.frame(cbind(substr(y$targets$IDATfile,1,16),y$targets$Participant,y$targets$Time))
colnames(samples)<-c("IDATfile","Participant","Time")
samplesb<-samples[samples$Participant==i & samples$Time=="Baseline",]
samplesa<-samples[samples$Participant==i & samples$Time=="Followup",]

#Change in genes
genesb<-as.data.frame(genes[,which(colnames(genes) %in% samplesb$IDATfile)])
genesb$bav<-rowMeans(genesb)
genesa<-as.data.frame(genes[,which(colnames(genes) %in% samplesa$IDATfile)])
genesa$aav<-rowMeans(genesa)
ba<-as.data.frame(cbind(genesb$bav,genesa$aav))
rownames(ba)<-rownames(genesa)
ba$change<-ba$V2-ba$V1
df[df$Participant==i,2:26554]<-ba$change

#Change in bmi
bmib<-as.data.frame(bmi$`y$targets$BMI/(sd(y$targets$BMI))`[bmi$IDATfile %in% samplesb$IDATfile])
bav<-colMeans(bmib)
bmia<-as.data.frame(bmi$`y$targets$BMI/(sd(y$targets$BMI))`[bmi$IDATfile %in% samplesa$IDATfile])
aav<-colMeans(bmia)
change<-aav-bav
df[df$Participant==i,26555]<-change

#Change in insulin
insb<-as.data.frame(ins$`as.numeric(y$targets$insulin)/(sd(na.omit(as.numeric(y$targets$insulin))))`[ins$IDATfile %in% samplesb$IDATfile])
bav<-colMeans(insb)
insa<-as.data.frame(ins$`as.numeric(y$targets$insulin)/(sd(na.omit(as.numeric(y$targets$insulin))))`[ins$IDATfile %in% samplesa$IDATfile])
aav<-colMeans(insa)
change<-aav-bav
df[df$Participant==i,26556]<-change

#WHR
WHRb<-as.data.frame(WHR$`as.numeric(y$targets$WHR)/(sd(na.omit(as.numeric(y$targets$WHR))))`[WHR$IDATfile %in% samplesb$IDATfile])
bav<-colMeans(WHRb)
WHRa<-as.data.frame(WHR$`as.numeric(y$targets$WHR)/(sd(na.omit(as.numeric(y$targets$WHR))))`[WHR$IDATfile %in% samplesa$IDATfile])
aav<-colMeans(WHRa)
change<-aav-bav
df[df$Participant==i,26557]<-change

#igf1
igf1b<-as.data.frame(igf1$`as.numeric(y$targets$igf1)/(sd(na.omit(as.numeric(y$targets$igf1))))`[igf1$IDATfile %in% samplesb$IDATfile])
bav<-colMeans(igf1b)
igf1a<-as.data.frame(igf1$`as.numeric(y$targets$igf1)/(sd(na.omit(as.numeric(y$targets$igf1))))`[igf1$IDATfile %in% samplesa$IDATfile])
aav<-colMeans(igf1a)
change<-aav-bav
df[df$Participant==i,26558]<-change

#crp
crpb<-as.data.frame(crp$`as.numeric(y$targets$crp)/(sd(na.omit(as.numeric(y$targets$crp))))`[crp$IDATfile %in% samplesb$IDATfile])
bav<-colMeans(crpb)
crpa<-as.data.frame(crp$`as.numeric(y$targets$crp)/(sd(na.omit(as.numeric(y$targets$crp))))`[crp$IDATfile %in% samplesa$IDATfile])
aav<-colMeans(crpa)
change<-aav-bav
df[df$Participant==i,26559]<-change

#glucose
glucoseb<-as.data.frame(glucose$`as.numeric(y$targets$glucose)/(sd(na.omit(as.numeric(y$targets$glucose))))`[glucose$IDATfile %in% samplesb$IDATfile])
bav<-colMeans(glucoseb)
glucosea<-as.data.frame(glucose$`as.numeric(y$targets$glucose)/(sd(na.omit(as.numeric(y$targets$glucose))))`[glucose$IDATfile %in% samplesa$IDATfile])
aav<-colMeans(glucosea)
change<-aav-bav
df[df$Participant==i,26560]<-change

#HOMA.IR
HOMA.IRb<-as.data.frame(HOMA.IR$`as.numeric(y$targets$HOMA.IR)/(sd(na.omit(as.numeric(y$targets$HOMA.IR))))`[HOMA.IR$IDATfile %in% samplesb$IDATfile])
bav<-colMeans(HOMA.IRb)
HOMA.IRa<-as.data.frame(HOMA.IR$`as.numeric(y$targets$HOMA.IR)/(sd(na.omit(as.numeric(y$targets$HOMA.IR))))`[HOMA.IR$IDATfile %in% samplesa$IDATfile])
aav<-colMeans(HOMA.IRa)
change<-aav-bav
df[df$Participant==i,26561]<-change

#hdl
hdlb<-as.data.frame(hdl$`as.numeric(y$targets$hdl)/(sd(na.omit(as.numeric(y$targets$hdl))))`[hdl$IDATfile %in% samplesb$IDATfile])
bav<-colMeans(hdlb)
hdla<-as.data.frame(hdl$`as.numeric(y$targets$hdl)/(sd(na.omit(as.numeric(y$targets$hdl))))`[hdl$IDATfile %in% samplesa$IDATfile])
aav<-colMeans(hdla)
change<-aav-bav
df[df$Participant==i,26562]<-change

#tc
tcb<-as.data.frame(tc$`as.numeric(y$targets$tc)/(sd(na.omit(as.numeric(y$targets$tc))))`[tc$IDATfile %in% samplesb$IDATfile])
bav<-colMeans(tcb)
tca<-as.data.frame(tc$`as.numeric(y$targets$tc)/(sd(na.omit(as.numeric(y$targets$tc))))`[tc$IDATfile %in% samplesa$IDATfile])
aav<-colMeans(tca)
change<-aav-bav
df[df$Participant==i,26563]<-change

#tg
tgb<-as.data.frame(tg$`as.numeric(y$targets$tg)/(sd(na.omit(as.numeric(y$targets$tg))))`[tg$IDATfile %in% samplesb$IDATfile])
bav<-colMeans(tgb)
tga<-as.data.frame(tg$`as.numeric(y$targets$tg)/(sd(na.omit(as.numeric(y$targets$tg))))`[tg$IDATfile %in% samplesa$IDATfile])
aav<-colMeans(tga)
change<-aav-bav
df[df$Participant==i,26564]<-change

#mTOR
mTORb<-as.data.frame(mTOR$`as.numeric(y$targets$mTOR)/(sd(na.omit(as.numeric(y$targets$mTOR))))`[mTOR$IDATfile %in% samplesb$IDATfile])
bav<-colMeans(mTORb)
mTORa<-as.data.frame(mTOR$`as.numeric(y$targets$mTOR)/(sd(na.omit(as.numeric(y$targets$mTOR))))`[mTOR$IDATfile %in% samplesa$IDATfile])
aav<-colMeans(mTORa)
change<-aav-bav
df[df$Participant==i,26565]<-change

#M30
M30b<-as.data.frame(M30$`as.numeric(y$targets$M30)/(sd(na.omit(as.numeric(y$targets$M30))))`[M30$IDATfile %in% samplesb$IDATfile])
bav<-colMeans(M30b)
M30a<-as.data.frame(M30$`as.numeric(y$targets$M30)/(sd(na.omit(as.numeric(y$targets$M30))))`[M30$IDATfile %in% samplesa$IDATfile])
aav<-colMeans(M30a)
change<-aav-bav
df[df$Participant==i,26566]<-change

#Ki67
Ki67b<-as.data.frame(Ki67$`as.numeric(y$targets$Ki67)/(sd(na.omit(as.numeric(y$targets$Ki67))))`[Ki67$IDATfile %in% samplesb$IDATfile])
bav<-colMeans(Ki67b)
Ki67a<-as.data.frame(Ki67$`as.numeric(y$targets$Ki67)/(sd(na.omit(as.numeric(y$targets$Ki67))))`[Ki67$IDATfile %in% samplesa$IDATfile])
aav<-colMeans(Ki67a)
change<-aav-bav
df[df$Participant==i,26567]<-change

#Apoptosis
Apoptosisb<-as.data.frame(Apoptosis$`as.numeric(y$targets$Apoptosis)/(sd(na.omit(as.numeric(y$targets$Apoptosis))))`[Apoptosis$IDATfile %in% samplesb$IDATfile])
bav<-colMeans(Apoptosisb)
Apoptosisa<-as.data.frame(Apoptosis$`as.numeric(y$targets$Apoptosis)/(sd(na.omit(as.numeric(y$targets$Apoptosis))))`[Apoptosis$IDATfile %in% samplesa$IDATfile])
aav<-colMeans(Apoptosisa)
change<-aav-bav
df[df$Participant==i,26568]<-change

}

#Limit to differentially expressed genes
all<-fread("Results.csv")
n<-round(0.02*nrow(all),digits=0)
de<-all[1:n,]
de$FC<-exp(de$logFC)
de<-de[de$FC>1.4 | de$FC<0.714,]
de<-de$Array_Address_Id
df2<-df[,c(1,which(colnames(df) %in% de))]
df<-cbind(df2,df[,26555:26568])

#Linear regression
res<-data.frame(matrix(nrow=14,ncol=314))
rownames(res)<-c("BMI","Insulin","WHR","igf1","crp","glucose","HOMA.IR","hdl","tc","tg","mTOR","M30","Ki67","Apoptosis")
colnames(res)<-colnames(df)[2:(ncol(df)-14)]

res2<-apply(df[2:(ncol(df)-14)],2,function(x){lm(x~df$BMI)})
res3<-apply(df[2:(ncol(df)-14)],2,function(x){lm(x~df$Insulin)})
res4<-apply(df[2:(ncol(df)-14)],2,function(x){lm(x~df$WHR)})
res5<-apply(df[2:(ncol(df)-14)],2,function(x){lm(x~df$igf1)})
res6<-apply(df[2:(ncol(df)-14)],2,function(x){lm(x~df$crp)})
res7<-apply(df[2:(ncol(df)-14)],2,function(x){lm(x~df$glucose)})
res8<-apply(df[2:(ncol(df)-14)],2,function(x){lm(x~df$HOMA.IR)})
res9<-apply(df[2:(ncol(df)-14)],2,function(x){lm(x~df$hdl)})
res10<-apply(df[2:(ncol(df)-14)],2,function(x){lm(x~df$tc)})
res11<-apply(df[2:(ncol(df)-14)],2,function(x){lm(x~df$tg)})
res12<-apply(df[2:(ncol(df)-14)],2,function(x){lm(x~df$mTOR)})
res13<-apply(df[2:(ncol(df)-14)],2,function(x){lm(x~df$M30)})
res14<-apply(df[2:(ncol(df)-14)],2,function(x){lm(x~df$Ki67)})
res15<-apply(df[2:(ncol(df)-14)],2,function(x){lm(x~df$Apoptosis)})

#Get coefficients
for (i in colnames(df)[2:(ncol(df)-14)]){
  print(i)
  res[1,which(colnames(res)==i)]<-res2[[i]]$coefficients[2]
  res[2,which(colnames(res)==i)]<-res3[[i]]$coefficients[2]
  res[3,which(colnames(res)==i)]<-res4[[i]]$coefficients[2]
  res[4,which(colnames(res)==i)]<-res5[[i]]$coefficients[2]
  res[5,which(colnames(res)==i)]<-res6[[i]]$coefficients[2]
  res[6,which(colnames(res)==i)]<-res7[[i]]$coefficients[2]
  res[7,which(colnames(res)==i)]<-res8[[i]]$coefficients[2]
  res[8,which(colnames(res)==i)]<-res9[[i]]$coefficients[2]
  res[9,which(colnames(res)==i)]<-res10[[i]]$coefficients[2]
  res[10,which(colnames(res)==i)]<-res11[[i]]$coefficients[2]
  res[11,which(colnames(res)==i)]<-res12[[i]]$coefficients[2]
  res[12,which(colnames(res)==i)]<-res13[[i]]$coefficients[2]
  res[13,which(colnames(res)==i)]<-res14[[i]]$coefficients[2]
  res[14,which(colnames(res)==i)]<-res15[[i]]$coefficients[2]
  }

row.names(res)<-c("BMI","Insulin","WHR","IGF-1","CRP","Glucose","HOMA-IR","HDL cholesterol","Total cholesterol","Triglyceride","mTOR","M30","Ki67","Apoptosis")

res<-data.matrix(res)

jpeg("Coclustering Heatmap.jpeg",width = 250, height = 200, units = "mm",res=300)
Heatmap(res,show_column_name=FALSE,row_dend_reorder=FALSE,heatmap_legend_param = list(title = "Coefficient"))
dev.off()

#Looking at ABHD11 in heatmap
df<-res[,colnames(res)=="6110592"]
jpeg("ABHD11 Heatmap.jpeg",width = 250, height = 200, units = "mm",res=300)
Heatmap(df,show_column_name=FALSE,row_dend_reorder=FALSE,heatmap_legend_param = list(title = "Coefficient"))
dev.off()

#Find genes which are the same and different between BMI and WHR
dat<-res[c(1,2),]
same1<-dat[,which(dat[1,]>0)]
same1<-same1[,which(same1[2,]>0)]
same2<-dat[,which(dat[1,]<0)]
same2<-same2[,which(same2[2,]<0)]
same<-cbind(same1,same2)

dif1<-dat[,which(dat[1,]>0)]
dif1<-dif1[,which(dif1[2,]<0)]
dif2<-dat[,which(dat[1,]<0)]
dif2<-dif2[,which(dif2[2,]>0)]
dif<-cbind(dif1,dif2)

samegenes<-data.frame(colnames(same))
difgenes<-data.frame(colnames(dif))

#Get symbol
genes<-dplyr::select(all,Symbol,Array_Address_Id)
samegenes<-merge(samegenes,genes,by.x="colnames.same.",by.y="Array_Address_Id",all.x=TRUE)
write.csv(samegenes,"Same.csv",quote=FALSE,row.names = FALSE)

difgenes<-merge(difgenes,genes,by.x="colnames.dif.",by.y="Array_Address_Id",all.x=TRUE)
write.csv(difgenes,"Different.csv",quote=FALSE,row.names = FALSE)

#Mean r2 for each factor
res2<-apply(df[2:(ncol(df)-14)],2,function(x){anova(lm(x~df$BMI))})
res3<-apply(df[2:(ncol(df)-14)],2,function(x){anova(lm(x~df$Insulin))})
res4<-apply(df[2:(ncol(df)-14)],2,function(x){anova(lm(x~df$WHR))})
res5<-apply(df[2:(ncol(df)-14)],2,function(x){anova(lm(x~df$igf1))})
res6<-apply(df[2:(ncol(df)-14)],2,function(x){anova(lm(x~df$crp))})
res7<-apply(df[2:(ncol(df)-14)],2,function(x){anova(lm(x~df$glucose))})
res8<-apply(df[2:(ncol(df)-14)],2,function(x){anova(lm(x~df$HOMA.IR))})
res9<-apply(df[2:(ncol(df)-14)],2,function(x){anova(lm(x~df$hdl))})
res10<-apply(df[2:(ncol(df)-14)],2,function(x){anova(lm(x~df$tc))})
res11<-apply(df[2:(ncol(df)-14)],2,function(x){anova(lm(x~df$tg))})
res12<-apply(df[2:(ncol(df)-14)],2,function(x){anova(lm(x~df$mTOR))})
res13<-apply(df[2:(ncol(df)-14)],2,function(x){anova(lm(x~df$M30))})
res14<-apply(df[2:(ncol(df)-14)],2,function(x){anova(lm(x~df$Ki67))})
res15<-apply(df[2:(ncol(df)-14)],2,function(x){anova(lm(x~df$Apoptosis))})


res<-data.frame(matrix(nrow=14,ncol=314))
rownames(res)<-c("BMI","Insulin","WHR","igf1","crp","glucose","HOMA.IR","hdl","tc","tg","mTOR","M30","Ki67","Apoptosis")
colnames(res)<-colnames(df)[2:(ncol(df)-14)]


for (i in colnames(df)[2:(ncol(df)-14)]){
  print(i)
  res[1,which(colnames(res)==i)]<-res2[[i]]$`Sum Sq`[1]/sum(res2[[i]]$`Sum Sq`,na.rm=T)
  res[2,which(colnames(res)==i)]<-res3[[i]]$`Sum Sq`[1]/sum(res3[[i]]$`Sum Sq`,na.rm=T)
  res[3,which(colnames(res)==i)]<-res4[[i]]$`Sum Sq`[1]/sum(res4[[i]]$`Sum Sq`,na.rm=T)
  res[4,which(colnames(res)==i)]<-res5[[i]]$`Sum Sq`[1]/sum(res5[[i]]$`Sum Sq`,na.rm=T)
  res[5,which(colnames(res)==i)]<-res6[[i]]$`Sum Sq`[1]/sum(res6[[i]]$`Sum Sq`,na.rm=T)
  res[6,which(colnames(res)==i)]<-res7[[i]]$`Sum Sq`[1]/sum(res7[[i]]$`Sum Sq`,na.rm=T)
  res[7,which(colnames(res)==i)]<-res8[[i]]$`Sum Sq`[1]/sum(res8[[i]]$`Sum Sq`,na.rm=T)
  res[8,which(colnames(res)==i)]<-res9[[i]]$`Sum Sq`[1]/sum(res9[[i]]$`Sum Sq`,na.rm=T)
  res[9,which(colnames(res)==i)]<-res10[[i]]$`Sum Sq`[1]/sum(res10[[i]]$`Sum Sq`,na.rm=T)
  res[10,which(colnames(res)==i)]<-res11[[i]]$`Sum Sq`[1]/sum(res11[[i]]$`Sum Sq`,na.rm=T)
  res[11,which(colnames(res)==i)]<-res12[[i]]$`Sum Sq`[1]/sum(res12[[i]]$`Sum Sq`,na.rm=T)
  res[12,which(colnames(res)==i)]<-res13[[i]]$`Sum Sq`[1]/sum(res13[[i]]$`Sum Sq`,na.rm=T)
  res[13,which(colnames(res)==i)]<-res14[[i]]$`Sum Sq`[1]/sum(res14[[i]]$`Sum Sq`,na.rm=T)
  res[14,which(colnames(res)==i)]<-res15[[i]]$`Sum Sq`[1]/sum(res15[[i]]$`Sum Sq`,na.rm=T)
  }

results<-data.frame(matrix(nrow=14,ncol=2))
rownames(results)<-rownames(res)
colnames(results)<-c("Mean r2","Median r2")

for (i in rownames(results)){
  a<-res[which(rownames(res)==i),]
  results[which(rownames(results)==i),1]<-mean(as.numeric(a[1,]))
  results[which(rownames(results)==i),2]<-median(as.numeric(a[1,]))
}
