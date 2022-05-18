
###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 09/09/20                  Last Modified: 02/03/22 ###
###########################################################################
###########################################################################
###                   hetperkb.R         		                     	###
###########################################################################

#### PREREQUISITES #####
# load packages
library(ggplot2)
library(ggsignif)
library(RcppCNPy)
library(tidyr)
library(reshape2)
library(ggpubr)
library(cowplot)
library(plotrix)
library(wesanderson)
library(plyr)
library(dplyr)
library(DataCombine)
library(DescTools)
library(rcartocolor)
library(car)
library(detectRUNS)
library(tidyverse)
library(randomcoloR)
library(IdeoViz)
library(gridExtra)
library(reshape2)


# Load colors

col_az <- "#CC3A8E" #carto_pal(12,"Vivid")[11]
col_wtx <- "#9f82ce" #carto_pal(7,"Purp")[5]
col_etx <-"#f2855d" #carto_pal(7,"Peach")[5]
col_mx <- "#3598FF"
cols <- c(col_az,col_wtx,col_etx, col_mx)
colss <-c(rep(col_az,28),rep(col_mx,4),rep(col_wtx,31),rep(col_etx,3))
cols_ll <- adjustcolor(cols, alpha.f = 0.2)
col_azwtx <- "#90CCC2"
cols <- c(col_az,col_wtx,col_etx, col_mx,col_azwtx)
cols_l <- adjustcolor(cols, alpha.f = 0.3)
setwd("~/Documents/Thesis_Research/Final Results/Ch3/revise/slim/final_output/csv/")

#R2: with modern contraction of TX pop

# get files

for (i in 1:100)
{
  f <- read.csv(paste("final/run",i,".postbot.csv",sep=""), header=F)
  colnames(f) <- c("gen","pop","fit","size","het","froh","vsdel","sdel", "model","wdel","vsben","sben","moben","wben")
  f <- f[-c(1,2),]
  f <- na.omit(f)
  name <- paste("run",i,".file",sep="")
  assign(name,f)
}

dataframe <- NULL
# Parse data by pop
 for (i in 1:50)
  {
    name <- paste("run",i,".df",sep="")
    f <- get(paste("run",i,".file",sep=""))
    init.fit.p1 <- f$fit[1]
    init.fit.p2 <- f$fit[2]
    df.p1 <- f[which(f$pop == "p1"),]
    df.p2 <- f[which(f$pop == "p2"),]
    
    p1 <- data.frame(Run=i,Pop="P1",Gen=df.p1$gen, Relfit=log2(df.p1$fit/init.fit.p1), Size=df.p1$size, Het=df.p1$het,Froh=df.p1$froh,Del=(df.p1$vsdel+df.p1$sdel)*10,
                     Ben=(df.p1$vsben+df.p1$sben),WDel=(df.p1$model+df.p1$wdel)*10,WBen=(df.p1$moben+df.p1$wben))
    p2 <- data.frame(Run=i,Pop="P2",Gen=df.p2$gen, Relfit=log2(df.p2$fit/init.fit.p2), Size=df.p2$size, Het=df.p2$het,Froh=df.p2$froh,Del=(df.p2$vsdel+df.p2$sdel)*10,
                     Ben=(df.p2$vsben+df.p2$sben),WDel=(df.p2$model+df.p2$wdel)*10,WBen=(df.p2$moben+df.p2$wben))
    
    df <- rbind(p1,p2)
    dataframe <- rbind(dataframe,df)
  }

df.p1 <- dataframe[which(dataframe$Pop == "P1"),]
df.p2 <- dataframe[which(dataframe$Pop == "P2"),]

unique(df.p1$Gen)

# Get mean
meanframe <- NULL
for (i in unique(df.p1$Gen))
{
  a <- df.p1[which(df.p1$Gen == i),]
  b <- df.p2[which(df.p2$Gen == i),]
  
  p1 <- data.frame(Pop="P1",Gen=i, Relfit=mean(a$Relfit), Size=mean(a$Size), Het=mean(a$Het),Froh=mean(a$Froh),Del=mean(a$Del),
                   Ben=mean(a$Ben),WDel=mean(a$WDel),WBen=mean(a$WBen))
  p2 <- data.frame(Pop="P2",Gen=i, Relfit=mean(b$Relfit), Size=mean(b$Size), Het=mean(b$Het),Froh=mean(b$Froh),Del=mean(b$Del),
                   Ben=mean(b$Ben),WDel=mean(b$WDel),WBen=mean(b$WBen))
  df <- rbind(p1,p2)
  meanframe <- rbind(meanframe,df)
}

# Plot data
dataframe$Gen <- dataframe$Gen-125003
meanframe$Gen <- meanframe$Gen-125003


p1.mean <- meanframe[which(meanframe$Pop == "P1"),]
p2.mean <- meanframe[which(meanframe$Pop == "P2"),]
df.p1 <- dataframe[which(dataframe$Pop == "P1"),]
df.p2 <- dataframe[which(dataframe$Pop == "P2"),]

#Population size
p0 <- ggplot(dataframe, aes(x=Gen, col=Pop)) +
  xlab("Generation Post-Bottleneck") +
  ylab("Population size")+
  theme_classic(base_size = 15)+
  ylim(0,8000)+
  scale_fill_manual(values=cols[1:2])+
  scale_x_continuous(limits=c(1247, 1500), breaks = c(1300,1500))+
  scale_color_discrete(name = "Population", labels = c("Pop1", "Pop2"))+
  theme(legend.title=element_blank(),
        axis.text.x=element_text(size=15),
        legend.position = "none")

p0 <- p0 + geom_line(aes(y=Size), color="black", data = p1.mean, size=3) +
  geom_line(aes(y=Size), color="black",data = p2.mean, size=3) +
  geom_line(aes(y=Size), color=col_az, data = p1.mean, size=2) +
  geom_line(aes(y=Size), color=col_wtx,data = p2.mean, size=2) 

p0

# Heterozygosity
p1 <- ggplot(dataframe, aes(x=Gen, col=Pop)) +
  xlab("Generation Post-Bottleneck") +
  ylab("Heterozygosity")+
  theme_classic(base_size = 22)+
  ylim(0.00055,0.00075)+
  scale_fill_manual(values=cols[1:2])+
  scale_x_continuous(limits=c(0, 1500), breaks = c(0,500,1000,1300,1500))+
  scale_color_discrete(name = "Population", labels = c("Pop1", "Pop2"))+
  theme(legend.title=element_blank(),
        axis.text.x=element_text(size=15),
        legend.position = "none")

for (i in 1:100)
{
  a <- df.p1[which(df.p1$Run == i),]
  b <- df.p2[which(df.p2$Run == i),]
  p1 <- p1 + geom_line(aes(x=Gen,y=Het),data = a, col=cols_l[1])+
    geom_line(aes(x=Gen,y=Het),data = b, col=cols_l[2])
}

p1 <- p1 + geom_line(aes(y=Het), color="black", data = p1.mean, size=3) +
  geom_line(aes(y=Het), color="black",data = p2.mean, size=3) +
  geom_line(aes(y=Het), color=col_az, data = p1.mean, size=2) +
  geom_line(aes(y=Het), color=col_wtx,data = p2.mean, size=2) 

p1 

# FROH
p2 <- ggplot(dataframe, aes(x=Gen, col=Pop)) +
  xlab("Generation Post-Bottleneck") +
  ylab(expression(F[ROH]))+
  theme_classic(base_size = 22)+
  ylim(0,0.15)+
  scale_fill_manual(values=cols[1:2])+
  scale_x_continuous(limits=c(0, 1500), breaks = c(0,500,1000,1300,1500))+
  scale_color_discrete(name = "Population", labels = c("Pop1", "Pop2"))+
  theme(legend.title=element_blank(),
        axis.text.x=element_text(size=15),
        legend.position = "none")

for (i in 1:100)
{
  a <- df.p1[which(df.p1$Run == i),]
  b <- df.p2[which(df.p2$Run == i),]
  p2 <- p2 + geom_line(aes(x=Gen,y=Froh),data = a, col=cols_l[1])+
    geom_line(aes(x=Gen,y=Froh),data = b, col=cols_l[2])
}

p2 <- p2 + geom_line(aes(y=Froh), color="black", data = p1.mean, size=3) +
  geom_line(aes(y=Froh), color="black",data = p2.mean, size=3) +
  geom_line(aes(y=Froh), color=col_az, data = p1.mean, size=2) +
  geom_line(aes(y=Froh), color=col_wtx,data = p2.mean, size=2) 

p2 

# Relative fitness
p3 <- ggplot(dataframe, aes(x=Gen, col=Pop)) +
  xlab("Generation Post-Bottleneck") +
  ylab(expression(paste("Relative Fitness = ",log[2],"(",frac(bar(F)[t(ind)], bar(F)[t0(ind)]),")")))+
  theme_classic(base_size = 22)+
  ylim(-1,2)+
  scale_fill_manual(values=cols[1:2])+
  scale_x_continuous(limits=c(0, 1500), breaks = c(0,500,1000,1300,1500))+
  scale_color_discrete(name = "Population", labels = c("Pop1", "Pop2"))+
  theme(legend.title=element_blank(),
        axis.text.x=element_text(size=15),
        legend.position = "none")

for (i in 1:100)
{
  a <- df.p1[which(df.p1$Run == i),]
  b <- df.p2[which(df.p2$Run == i),]
  p3 <- p3 + geom_line(aes(x=Gen,y=Relfit),data = a, col=cols_l[1])+
    geom_line(aes(x=Gen,y=Relfit),data = b, col=cols_l[2])
}

p3 <- p3 + geom_line(aes(y=Relfit), color="black", data = p1.mean, size=3) +
  geom_line(aes(y=Relfit), color="black",data = p2.mean, size=3) +
  geom_line(aes(y=Relfit), color=col_az, data = p1.mean, size=2) +
  geom_line(aes(y=Relfit), color=col_wtx,data = p2.mean, size=2) 

p3

# Strongly Deleterious mutations
p4 <- ggplot(dataframe, aes(x=Gen, col=Pop)) +
  xlab("Generation Post-Bottleneck") +
  ylab(expression(paste("Strongly deleterious")))+
  theme_classic(base_size = 20)+
  ylim(0,300)+
  scale_fill_manual(values=cols[1:2])+
  scale_x_continuous(limits=c(0, 1500), breaks = c(0,500,1000,1300,1500))+
  scale_color_discrete(name = "Population", labels = c("Pop1", "Pop2"))+
  theme(legend.title=element_blank(),
        axis.text.x = element_text(size=15),
        legend.position = "none")

for (i in 1:100)
{
  a <- df.p1[which(df.p1$Run == i),]
  b <- df.p2[which(df.p2$Run == i),]
  p4 <- p4 + geom_line(aes(x=Gen,y=Del),data = a, col=cols_l[1])+
    geom_line(aes(x=Gen,y=Del),data = b, col=cols_l[2])
}

p4 <- p4 + geom_line(aes(y=Del), color="black", data = p1.mean, size=3) +
  geom_line(aes(y=Del), color="black",data = p2.mean, size=3) +
  geom_line(aes(y=Del), color=col_az, data = p1.mean, size=2) +
  geom_line(aes(y=Del), color=col_wtx,data = p2.mean, size=2)

p4

# Weakly Deleterious mutations

p5 <- ggplot(dataframe, aes(x=Gen, col=Pop)) +
  xlab("Generation Post-Bottleneck")+
  ylab(expression(paste("Weakly deleterious")))+
  theme_classic(base_size = 20)+
  ylim(1900,2500)+
  scale_fill_manual(values=cols[1:2])+
  scale_x_continuous(limits=c(0, 1500), breaks = c(0,500,1000,1300,1500))+
  scale_color_discrete(name = "Population", labels = c("Pop1", "Pop2"))+
  theme(legend.title=element_blank(),
        axis.text.x = element_text(size=15),
        legend.position = "none")

for (i in 1:100)
{
  a <- df.p1[which(df.p1$Run == i),]
  b <- df.p2[which(df.p2$Run == i),]
  p5 <- p5 + geom_line(aes(x=Gen,y=WDel),data = a, col=cols_l[1])+
    geom_line(aes(x=Gen,y=WDel),data = b, col=cols_l[2])
}

p5 <- p5 + geom_line(aes(y=WDel), color="black", data = p1.mean, size=3) +
  geom_line(aes(y=WDel), color="black",data = p2.mean, size=3) +
  geom_line(aes(y=WDel), color=col_az, data = p1.mean, size=2) +
  geom_line(aes(y=WDel), color=col_wtx,data = p2.mean, size=2)

p5

# Strongly Beneficial mutations
p6 <- ggplot(dataframe, aes(x=Gen, col=Pop)) +
  xlab("Generation Post-Bottleneck") +
  ylab(expression(paste("Strongly beneficial")))+
  theme_classic(base_size = 20)+
  scale_fill_manual(values=cols[1:2])+
  scale_x_continuous(limits=c(0, 1500), breaks = c(0,500,1000,1300,1500))+
  scale_y_continuous(limits=c(164, 170), breaks = c(164,165,166,167,168,169,170))+
  scale_color_discrete(name = "Population", labels = c("Pop1", "Pop2"))+
  theme(legend.title=element_blank(),
        axis.text.x = element_text(size=15),
        legend.position = "none")

for (i in 1:100)
{
  a <- df.p1[which(df.p1$Run == i),]
  b <- df.p2[which(df.p2$Run == i),]
  p6 <- p6 + geom_line(aes(x=Gen,y=Ben),data = a, col=cols_l[1])+
    geom_line(aes(x=Gen,y=Ben),data = b, col=cols_l[2])
}

p6 <- p6 + geom_line(aes(y=Ben), color="black", data = p1.mean, size=3) +
  geom_line(aes(y=Ben), color="black",data = p2.mean, size=3) +
  geom_line(aes(y=Ben), color=col_az, data = p1.mean, size=2) +
  geom_line(aes(y=Ben), color=col_wtx,data = p2.mean, size=2)

p6

# Weakly Beneficial mutations

p7 <- ggplot(dataframe, aes(x=Gen, col=Pop)) +
  xlab("Generation Post-Bottleneck")+
  ylab(expression(paste("Weakly beneficial")))+
  theme_classic(base_size = 20)+
  ylim(140,190)+
  scale_fill_manual(values=cols[1:2])+
  scale_x_continuous(limits=c(0, 1500), breaks = c(0,500,1000,1300,1500))+
  scale_color_discrete(name = "Population", labels = c("Pop1", "Pop2"))+
  theme(legend.title=element_blank(),
        axis.text.x = element_text(size=15),
        legend.position = "none")

for (i in 1:100)
{
  a <- df.p1[which(df.p1$Run == i),]
  b <- df.p2[which(df.p2$Run == i),]
  p7 <- p7 + geom_line(aes(x=Gen,y=WBen),data = a, col=cols_l[1])+
    geom_line(aes(x=Gen,y=WBen),data = b, col=cols_l[2])
}

p7 <- p7 + geom_line(aes(y=WBen), color="black", data = p1.mean, size=3) +
  geom_line(aes(y=WBen), color="black",data = p2.mean, size=3) +
  geom_line(aes(y=WBen), color=col_az, data = p1.mean, size=2) +
  geom_line(aes(y=WBen), color=col_wtx,data = p2.mean, size=2)

p7

#### FST analysis ####

fst_all <- read.table("../fst_allmuts.txt",header = T)
fst_ben <- read.table("../fst_ben.txt",header = T)
fst_del <- read.table("../fst_del.txt",header = T)
fst_neut <- read.table("../fst_neut.txt",header = T)

fst_all_sites <- read.table("../fst_allmuts.sites.txt",header = T)
fst_ben_sites <- read.table("../fst_ben.sites.txt",header = T)
fst_del_sites <- read.table("../fst_del.sites.txt",header = T)
fst_neut_sites <- read.table("../fst_neut.sites.txt",header = T)


gens <- c(125002,125252,125502,125752,126002,126302,126502)
colnames(fst_all) <- gens
colnames(fst_ben) <- gens
colnames(fst_del) <- gens
colnames(fst_neut) <- gens
colnames(fst_all_sites) <- gens
colnames(fst_ben_sites) <- gens
colnames(fst_del_sites) <- gens
colnames(fst_neut_sites) <- gens

df.fst <- NULL
for (i in 1:7)
{
  df <- data.frame(Gen=as.character(gens[i]-125002),All=fst_all[,i],Ben=fst_ben[,i],Del=fst_del[,i],Neut=fst_neut[,i],
                   AllSites=fst_all_sites[,i],BenSites=fst_ben_sites[,i],DelSites=fst_del_sites[,i],NeutSites=fst_neut_sites[,i])
  df.fst <- rbind(df.fst,df)
}

df.fst$Gen <- factor(df.fst$Gen,levels=as.character(gens-125002))
fst.mean <- ddply(df.fst, "Gen", summarise, Mean_All=mean(All),Mean_Ben=mean(Ben),Mean_Del=mean(Del),Mean_Neut=mean(Neut),
                  SD_All=sd(All),SD_Ben=sd(Ben),SD_Del=sd(Del),SD_Neut=sd(Neut),
                  Mean_AllSites=mean(AllSites),Mean_BenSites=mean(BenSites),Mean_DelSites=mean(DelSites),Mean_NeutSites=mean(NeutSites),
                  SD_AllSites=sd(AllSites),SD_BenSites=sd(BenSites),SD_DelSites=sd(DelSites),SD_NeutSites=sd(NeutSites))
fst.mean$Gen <- as.character(fst.mean$Gen)
fst.mean$Gen <- factor(fst.mean$Gen,levels=as.character(gens-125002))

# FST
ggplot(fst.mean,aes(x=Gen))+
  geom_line(aes(y=Mean_All,group=1),col="black",size=2)+
  geom_line(aes(y=Mean_Ben,group=1),col="#00BFC4",size=2)+
  geom_line(aes(y=Mean_Del,group=1),col="#F8766D",size=2)+
  geom_line(aes(y=Mean_Neut,group=1),col="gray50",size=2)+
  
  geom_errorbar(aes(ymin=Mean_All-SD_All, ymax=Mean_All+SD_All), col="black",width=.2,size=1.5,position=position_dodge(.9))+
  geom_errorbar(aes(ymin=Mean_Ben-SD_Ben, ymax=Mean_Ben+SD_Ben), col="#00BFC4",width=.2,size=1.5,position=position_dodge(.9))+
  geom_errorbar(aes(ymin=Mean_Del-SD_Del, ymax=Mean_Del+SD_Del), col="#F8766D",width=.2,size=1.5,position=position_dodge(.9))+
  geom_errorbar(aes(ymin=Mean_Neut-SD_Neut, ymax=Mean_Neut+SD_Neut), col="gray50",width=.2,size=1.5,position=position_dodge(.9))+
  
  geom_point(aes(y=Mean_All),col="black",pch=21, size=4, fill="black")+
  geom_point(aes(y=Mean_Ben),col="black",pch=21, size=4, fill="#00BFC4")+
  geom_point(aes(y=Mean_Del),col="black",pch=21, size=4, fill="#F8766D")+
  geom_point(aes(y=Mean_Neut),col="black",pch=21, size=4, fill="gray50")+
  
  theme_classic(base_size = 22)+
  labs(x="Generation Post-Bottleneck", 
       y=expression(paste(F[ST])))+
  theme(axis.title.x = element_text(size=15),
    legend.title=element_blank(),
        legend.position = "none")


kruskal.test(All ~ Gen, data = df.fst) # Kruskal-Wallis chi-squared = 616.1, df = 6, p-value < 2.2e-16
pairwise.wilcox.test(df.fst$Ben, df.fst$Gen,
                     p.adjust.method = "fdr")

#Sites

ggplot(fst.mean,aes(x=Gen))+
  geom_line(aes(y=Mean_AllSites,group=1),col="black",size=2)+
  geom_line(aes(y=Mean_BenSites,group=1),col="#00BFC4",size=2)+
  geom_line(aes(y=Mean_DelSites,group=1),col="#F8766D",size=2)+
  geom_line(aes(y=Mean_NeutSites,group=1),col="gray50",size=2)+
  
  geom_errorbar(aes(ymin=Mean_AllSites-SD_AllSites, ymax=Mean_AllSites+SD_AllSites), col="black",width=.2,size=1.5,position=position_dodge(.9))+
  geom_errorbar(aes(ymin=Mean_BenSites-SD_BenSites, ymax=Mean_BenSites+SD_BenSites), col="#00BFC4",width=.2,size=1.5,position=position_dodge(.9))+
  geom_errorbar(aes(ymin=Mean_DelSites-SD_DelSites, ymax=Mean_DelSites+SD_DelSites), col="#F8766D",width=.2,size=1.5,position=position_dodge(.9))+
  geom_errorbar(aes(ymin=Mean_NeutSites-SD_NeutSites, ymax=Mean_NeutSites+SD_NeutSites), col="gray50",width=.2,size=1.5,position=position_dodge(.9))+
  
  geom_point(aes(y=Mean_AllSites),col="black",pch=21, size=4, fill="black")+
  geom_point(aes(y=Mean_BenSites),col="black",pch=21, size=4, fill="#00BFC4")+
  geom_point(aes(y=Mean_DelSites),col="black",pch=21, size=4, fill="#F8766D")+
  geom_point(aes(y=Mean_NeutSites),col="black",pch=21, size=4, fill="gray50")+
  
  theme_classic(base_size = 22)+
  labs(x="Generation Post-Bottleneck", 
       y=expression(paste("No. of SNPs")))+
  theme(axis.title.x = element_text(size=15),
        legend.title=element_blank(),
        legend.position = "none")



df.fst <- NULL
for (i in 1:7)
{
  df <- data.frame(Gen=gens[i]-107502,Fst=fst_all[,i])
  df.fst <- rbind(df.fst,df)
}
fst.mean <- ddply(df.fst, "Gen", summarise, Mean_Fst=mean(Fst))
pairwise.wilcox.test(df.fst$Fst, df.fst$Gen,
                     p.adjust.method = "fdr")

dft <-diff(fst.mean$Mean_Fst)/diff(fst.mean$Gen)
ggplot()+
  geom_point(aes(x=Gen,y=Mean_Fst),data=fst.mean, col="black",pch=22, size=4, fill="red")

plot(1:5,dft,type="b")


#### Rxy analysis ####

gens <- c(125002,125252,125502,125752,126002,126302,126502)
rxy <- read.table("../Ryx.sims.txt",header = T)
colnames(rxy)[c(3,4)] <- c("Rxy_Del","Rxy_Ben")
rxy$Gen <- rxy$Gen-125002
rxy$Gen <- as.character(rxy$Gen)

rxy$Gen <- factor(rxy$Gen,levels=as.character(gens-125002))

rxy.df <- rbind(data.frame(Gen=rxy$Gen,Run=rxy$Run,Mutation="Deleterious",Rxy=rxy$Rxy_Del),
                data.frame(Gen=rxy$Gen,Run=rxy$Run,Mutation="Beneficial",Rxy=rxy$Rxy_Ben))

rxy_long <- melt(rxy.df, id.vars=c("Gen","Run","Mutation"))


rxy_long$Mutation <- factor(rxy_long$Mutation,levels=c("Deleterious","Beneficial"))
rxy_long$Gen <- factor(rxy_long$Gen,levels=as.character(gens-125002))
ggplot(rxy_long, aes(x=factor(Gen),y=1/value,fill=factor(Mutation)))+
  geom_point(position=position_jitterdodge(), size=2,aes(color=Gen),color="black",pch=21)+
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  geom_hline(yintercept = 1, linetype=2)+
  theme_classic(base_size = 22)+
  scale_fill_discrete(name = "Type")+
  labs(x="Generation Post-Bottleneck", 
       y=expression(paste(R[XY])))+
  theme(legend.title=element_blank(),
        legend.position = "none")
  

pairwise.wilcox.test(rxy$Rxy_Del, rxy$Gen,
                     p.adjust.method = "fdr")

pairwise.wilcox.test(rxy$Rxy_Ben, rxy$Gen,
                     p.adjust.method = "fdr")

#mean and sd
rxy.mean <- NULL
for (gen in c(0,250,500,750,1000,1300,1500))
{
  df <- rbind(data.frame(Gen=gen,Mutation="Deleterious", meanRxy=mean(1/rxy$Rxy_Del[which(rxy$Gen==gen)]),
                   sdRxy=sd(1/rxy$Rxy_Del[which(rxy$Gen==gen)])),
              data.frame(Gen=gen,Mutation="Beneficial", meanRxy=mean(1/rxy$Rxy_Ben[which(rxy$Gen==gen)]),
                         sdRxy=sd(1/rxy$Rxy_Ben[which(rxy$Gen==gen)])))
  rxy.mean <- rbind(rxy.mean,df)
}

rxy.mean

# Mutation age analysis

age <- read.table("../Age.sims.txt",header = T)
age_long <- melt(age, id.vars=c("Run","Type","Site"))

age$Age <- abs(age$Age-126502)

Snps <- NULL
for (i in 1:100)
{
  b_s <- nrow(age[which(age$Run ==i & age$Type=="Beneficial" & age$Site =="Shared"),])
  b_paz <- nrow(age[which(age$Run ==i & age$Type=="Beneficial" & age$Site =="Private_AZ"),])
  b_ptx <- nrow(age[which(age$Run ==i & age$Type=="Beneficial" & age$Site =="Private_TX"),])
  d_s <- nrow(age[which(age$Run ==i & age$Type=="Deleterious" & age$Site =="Shared"),])
  d_paz <- nrow(age[which(age$Run ==i & age$Type=="Deleterious" & age$Site =="Private_AZ"),])
  d_ptx <- nrow(age[which(age$Run ==i & age$Type=="Deleterious" & age$Site =="Private_TX"),])
  snps <- c(i,d_s,d_paz,d_ptx,b_s,b_paz,b_ptx)
  Snps <- as.data.frame(rbind(Snps,snps))
}

meanSnps <- rbind(data.frame(Site="AZ",Type="Deleterious",Mean=mean(c(Snps$V2,Snps$V3)),SD=sd(c(Snps$V2,Snps$V3))),
                  data.frame(Site="WTX",Type="Deleterious",Mean=mean(c(Snps$V2,Snps$V4)),SD=sd(c(Snps$V2,Snps$V4))),
                  data.frame(Site="AZ",Type="Beneficial",Mean=mean(c(Snps$V5,Snps$V6)),SD=sd(c(Snps$V5,Snps$V6))),
                  data.frame(Site="WTX",Type="Beneficial",Mean=mean(c(Snps$V5,Snps$V7)),SD=sd(c(Snps$V5,Snps$V7))))

wilcox.test(c(Snps$V2,Snps$V3), c(Snps$V2,Snps$V4)) # W = 25326, p-value = 4.09e-06
wilcox.test(c(Snps$V5,Snps$V6), c(Snps$V5,Snps$V7))

meanSnps$Type <- factor(meanSnps$Type,levels=c("Deleterious","Beneficial"))
ggplot(meanSnps,aes(x=Type, y=Mean, fill=Site)) +
  geom_bar(position=position_dodge(), stat="identity", alpha=0.7, color="black") +
  geom_errorbar( aes(x=Type, ymin=Mean-SD, ymax=Mean+SD), position=position_dodge(.9),colour="black", alpha=0.9, width=0.3)+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:2])+
  labs(x="Type of mutation", 
       y=expression("No. of SNPs"))+
  theme(legend.title=element_blank(),
        legend.position = "none")
  

age1 <- age[which(age$Run ==1),]
del <- age1[which(age1$Type=="Deleterious"),]
ben <- age1[which(age1$Type=="Beneficial"),]

azDel <- del[which(del$MAF1 !=0),]
txDel <- del[which(del$MAF2 !=0),]
azBen <- ben[which(ben$MAF1 !=0),]
txBen <- ben[which(ben$MAF2 !=0),]

df <- rbind(data.frame(Site="AZ",Type="Deleterious", Age=azDel$Age,MAF=azDel$MAF1),
            data.frame(Site="AZ",Type="Beneficial", Age=azBen$Age,MAF=azBen$MAF1),
            data.frame(Site="TX",Type="Deleterious", Age=txDel$Age,MAF=txDel$MAF2),
            data.frame(Site="TX",Type="Beneficial", Age=txBen$Age,MAF=txBen$MAF2))

df$Type <- factor(df$Type,levels=c("Deleterious","Beneficial"))

df2 <- df[which(df$Type=="Deleterious"),]

ggplot(df, aes(x=Type,y=Age,fill=Site))+
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  #geom_point(position=position_jitterdodge(), size=2,aes(color=Gen),color="black",pch=21)+
  ylim(0,150000)+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:2])+
  labs(x="Type of mutation", 
       y=expression("Age of mutation"))+
  theme(legend.title=element_blank(),
        legend.position = "none")

wilcox.test(df$Age[which(df$Site=="AZ",df$Type=="Deleterious")],
            df$Age[which(df$Site=="TX",df$Type=="Deleterious")]) # W = 28142764, p-value = 1.454e-08

wilcox.test(df$Age[which(df$Site=="AZ",df$Type=="Beneficial")],
            df$Age[which(df$Site=="TX",df$Type=="Beneficial")]) # W = 28142764, p-value = 1.454e-08

ggplot(df, aes(x=Type,y=MAF,fill=Site))+
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  #geom_point(position=position_jitterdodge(), size=2,aes(color=Gen),color="black",pch=21)+
  ylim(0,1)+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:2])+
  labs(x="Type of mutation", 
       y=expression("MAF of mutation"))+
  theme(legend.title=element_blank(),
        legend.position = "none")

wilcox.test(df$MAF[which(df$Site=="AZ",df$Type=="Deleterious")],
            df$MAF[which(df$Site=="TX",df$Type=="Deleterious")]) # W = W = 23921820, p-value < 2.2e-16

wilcox.test(df$MAF[which(df$Site=="AZ",df$Type=="Beneficial")],
            df$MAF[which(df$Site=="TX",df$Type=="Beneficial")]) # W = 28142764, p-value = 1.454e-08

#Age plot
del_age <- del[which(del$variable=="Age"),]
del_age$value <- abs(del_age$value-126502)

p1 <- ggplot(del_age) +
  xlab(expression("Age of mutation (Ya)"))+
  ylab("No. of deleterious mutations")+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[c(1,2,5)])+
  ylim(0,500)+
  scale_x_continuous(limits=c(0, 4500), breaks = c(0,200,750,1500,3000,4500))+
  theme(legend.title=element_blank(),
        legend.position = "none")

for (i in 1:100)
{
  df <- del_age[which(del_age$Run == i),]
  df_privaz <- df[which(df$Site == "Private_AZ"),]
  df_privaz <- df_privaz[which(df_privaz$value <= 4500),]
  df_privtx <- df[which(df$Site == "Private_TX"),]
  df_privtx <- df_privtx[which(df_privtx$value <= 4500),]
  p1 <- p1 + geom_histogram(alpha=0.2, bins = 100, color=cols_l[1],fill=cols[1],aes(x=value),data=df_privaz)
  p1 <- p1 + geom_histogram(alpha=0.4, bins = 100, color=cols_l[2],fill=cols[2],aes(x=value),data=df_privtx)
}

p1 <- p1+ geom_vline(xintercept = 200, linetype=2)+ geom_vline(xintercept = 1500, linetype=2)
p1
  

# Shared MAF

shared.del <- age[which(age$Site=="Shared" & age$Type=="Deleterious"),]

az.del <- data.frame(Site="AZ",Age=shared.del$Age,MAF=shared.del$MAF1)
tx.del <- data.frame(Site="TX",Age=shared.del$Age,MAF=shared.del$MAF2)
df.del <- rbind(az.del,tx.del)
az.del$Age <- abs(az.del$Age-126502)
tx.del$Age <- abs(tx.del$Age-126502)
df.del$Age <- abs(df.del$Age-126502)

p1 <- ggplot(df.del,aes(fill=Site, y=MAF, x=Age))+
  ylim(0,1)+
  geom_point(color=adjustcolor(col_az,alpha.f = 0.2), size=0.5, shape=16)+ 
  geom_point(color=adjustcolor(col_wtx,alpha.f = 0.2),  size=0.5, shape=16)+ 
  scale_x_reverse(limits=c(2000, 0), breaks = c(0,200,750,1500,2000))+
  geom_vline(xintercept = 200, linetype=2)+ geom_vline(xintercept = 1500, linetype=2)+
  theme_classic(base_size = 15)+
  labs(y="Derived Allele Frequency", x="Years before present")+
  theme(legend.title=element_blank(),
        legend.position = "none")

p2 <- p1 +geom_smooth(method=loess, span=0.80, aes(x=Age, y=MAF), data = az.del, color=col_az, fill="#69b3a2", se=F, size=2.5)+
  geom_smooth(method=loess, span=0.80, aes(x=Age, y=MAF), data = tx.del, color=col_wtx, fill="#69b3a2", se=F, size=2.5)+
  geom_smooth(method=loess, span=0.80, aes(x=Age, y=MAF), data = az.del, color="black", fill="#69b3a2", se=F, size=0.5)+
  geom_smooth(method=loess, span=0.80, aes(x=Age, y=MAF), data = tx.del, color="black", fill="#69b3a2", se=F, size=0.5)

p2

# Private MAF
privaz.del <- age[which(age$Site=="Private_AZ" & age$Type=="Deleterious"),]
privtx.del <- age[which(age$Site=="Private_TX" & age$Type=="Deleterious"),]

az.del <- data.frame(Site="AZ",Age=privaz.del$Age,MAF=privaz.del$MAF1)
tx.del <- data.frame(Site="TX",Age=privtx.del$Age,MAF=privtx.del$MAF2)
df.del <- rbind(az.del,tx.del)
az.del$Age <- abs(az.del$Age-126502)
tx.del$Age <- abs(tx.del$Age-126502)
df.del$Age <- abs(df.del$Age-126502)

private.del <- del2[-which(del2$Site=="Shared"),]
private.age <- private.del[which(private.del$variable=="Age"),]
private.maf1 <- private.del[which(private.del$variable=="MAF1"),]
private.maf2 <- private.del[which(private.del$variable=="MAF2"),]

p1 <- ggplot(df.del,aes(fill=Site, y=MAF, x=Age))+
  ylim(0,1)+
  geom_point(color=adjustcolor(col_az,alpha.f = 0.2), size=0.5, shape=16)+ 
  geom_point(color=adjustcolor(col_wtx,alpha.f = 0.2),  size=0.5, shape=16)+ 
  scale_x_reverse(limits=c(2000, 0), breaks = c(0,200,750,1500,2000))+
  geom_vline(xintercept = 200, linetype=2)+ geom_vline(xintercept = 1500, linetype=2)+
  theme_classic(base_size = 15)+
  labs(y="Derived Allele Frequency", x="Years before present")+
  theme(legend.title=element_blank(),
        legend.position = "none")

p2 <- p1 +geom_smooth(method=loess, span=0.80, aes(x=Age, y=MAF), data = az.del, color=col_az, fill="#69b3a2", se=F, size=2.5)+
  geom_smooth(method=loess, span=0.80, aes(x=Age, y=MAF), data = tx.del, color=col_wtx, fill="#69b3a2", se=F, size=2.5)+
  geom_smooth(method=loess, span=0.80, aes(x=Age, y=MAF), data = az.del, color="black", fill="#69b3a2", se=F, size=0.5)+
  geom_smooth(method=loess, span=0.80, aes(x=Age, y=MAF), data = tx.del, color="black", fill="#69b3a2", se=F, size=0.5)

p2

# Het and FROH for final gen

finalgen <- dataframe[which(dataframe$Gen==1499),]

ggplot(finalgen, aes(x=Pop,y=Het,fill=Pop))+
  geom_point(position=position_jitterdodge(), size=2,aes(color=Gen),color="black",pch=21)+
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  theme_classic(base_size = 22)+
  scale_fill_discrete(name = "Population", labels = c("AZ-like", "TX-like"))+
  scale_color_discrete(name = "Population", labels = c("AZ-like", "TX-like"))+
  labs(x="Population", 
       y=expression(paste(Het[exome])))+
  theme(legend.title=element_blank(),
        legend.position = "none")

ggplot(finalgen, aes(x=Pop,y=Froh,fill=Pop))+
  geom_point(position=position_jitterdodge(), size=2,aes(color=Gen),color="black",pch=21)+
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  theme_classic(base_size = 22)+
  scale_fill_discrete(name = "Population", labels = c("AZ-like", "TX-like"))+
  scale_color_discrete(name = "Population", labels = c("AZ-like", "TX-like"))+
  labs(x="Population", 
       y=expression(paste(F[ROH])))+
  theme(legend.title=element_blank(),
        legend.position = "none")

c(mean(finalgen$Het[which(finalgen$Pop=="P1")]),sd(finalgen$Het[which(finalgen$Pop=="P1")])) # 7.661958e-04 9.491077e-06
c(mean(finalgen$Het[which(finalgen$Pop=="P2")]),sd(finalgen$Het[which(finalgen$Pop=="P2")])) # 6.824299e-04 1.906859e-05

c(mean(finalgen$Froh[which(finalgen$Pop=="P1")]),sd(finalgen$Froh[which(finalgen$Pop=="P1")])) # 0.083645753 0.004909685
c(mean(finalgen$Froh[which(finalgen$Pop=="P2")]),sd(finalgen$Froh[which(finalgen$Pop=="P2")])) # 0.14551385 0.01332167
