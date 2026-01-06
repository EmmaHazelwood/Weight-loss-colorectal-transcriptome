library(dplyr)
library(ggplot2)
library(data.table)
library(ggpubr)

dat <- fread("Figures/TCGA.csv")
dat$Upper_whisker <- dat$`Upper whisker`
dat$Gene <- factor(dat$Gene, levels = c("ABHD11","ATP5MC2","CHMP2A"))

rectal<-dat[dat$Site=="Rectum",]
colon<-dat[dat$Site=="Colon",]

# Make a violin plot
p1 <- ggplot(rectal,aes(group=Gene,x=as.factor(Tissue),ymin=Min,lower=Q1,middle=Med,upper=Q3,ymax=Upper_whisker,fill=Tissue)) + 
  geom_boxplot(stat="identity",aes(group = interaction(Gene,Tissue))) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Gene expression value") +
  facet_wrap(~Gene,scales = "free",nrow=1) + 
  theme(text = element_text(size=20))+
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  scale_fill_manual(values = c("#317EC2","#C03830"))+
  theme(legend.position = "bottom")
p1

p2 <- ggplot(colon,aes(group=Gene,x=as.factor(Tissue),ymin=Min,lower=Q1,middle=Med,upper=Q3,ymax=Upper_whisker,fill=Tissue)) + 
  geom_boxplot(stat="identity",aes(group = interaction(Gene,Tissue))) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="Gene expression value") +
  facet_wrap(~Gene,scales = "free",nrow=1) + 
  theme(text = element_text(size=20))+
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  scale_fill_manual(values = c("#317EC2","#C03830")) +
  theme(legend.position = "none")
p2

p<-ggarrange(p2,NULL,p1,heights=c(5,0.3,5),ncol=1)
p


ggsave(filename="TCGA Plot.png", plot=last_plot(),width = 0.35*440, height = 0.5*270, units = "mm",bg="white")





