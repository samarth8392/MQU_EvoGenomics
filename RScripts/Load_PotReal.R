###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 09/09/20                  Last Modified: 01/18/22 ###
###########################################################################
###########################################################################
###                   Load_segDrift.R          		                      ###
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
library(viridis)
library(scales)
library(gridExtra)
# Load colors

# Load colors
col_az <- carto_pal(12,"Vivid")[11]
col_wtx <- carto_pal(7,"Purp")[5]
col_etx <-carto_pal(7,"Peach")[5]
col_mx <- "#3598FF"
cols <- c(col_az,col_wtx,col_etx, col_mx)
colss <-c(rep(col_az,28),rep(col_mx,4),rep(col_wtx,31),rep(col_etx,3))
cols_l <- adjustcolor(cols, alpha.f = 0.5)
col_azwtx <- "#90CCC2"

#### PREREQUISITES: END #####

setwd("~/Documents/Thesis_Research/Final Results/Ch3/revise/load/")

datafile <- read.table("LoadDatafile.txt",header = T)

datafile$Site <- factor(datafile$Site , levels=c("AZ", "WTX", "CTX", "MX"))
datafile$Type <- factor(datafile$Type , levels=c("Deleterious","Tolerated","Nonsynonymous","Synonymous"))

df <- datafile[-which(datafile$Type=="Synonymous"),]
df2 <- datafile[which(datafile$Type=="Synonymous"),]

df <- cbind(df,df$Het/df$Snps,df$Althom/df$Snps)
df2 <- cbind(df2,df2$Het/df2$Snps,df2$Althom/df2$Snps)

colnames(df)[8:9] <- c("PropHet","PropHomo")
colnames(df2)[8:9] <- c("PropHet","PropHomo")


my_comparisons <- list( c("AZ", "WTX"), c("AZ", "CTX"), c("AZ", "MX"),
                        c("WTX", "CTX"),c("WTX", "MX"), c("CTX", "MX"))

my_comparisons <- list(c("AZ", "WTX"))



df3 <- df[-which(df$Site == "CTX" | df$Site == "MX"),]


pal <- c(cols[1],cols_l[1],cols[2],cols_l[2])

# SNPs
ggplot(df,aes(x=Type,y=Snps,fill=Site))+
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  geom_point(position=position_jitterdodge(), size=2,aes(color=Site),color="black",pch=21)+
  labs(x="Site", y = "No. of SNPs")+
  #scale_x_discrete(labels=c("Deleterious" = "Deleterious", "Tolerated" = "Tolerated",
  #                          "Nonsynonymous" = "Functional"))+
  scale_y_continuous(labels = scales::comma)+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:4])+
  scale_color_manual(values=cols[1:4])+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position = "none")

# Statistical test for snps
df.all <- df[which(df$Type=="Nonsynonymous"),]
ddply(df.all,"Site",summarise,Mean=mean(Althom),SD=sd(Althom))
kruskal.test(Snps ~ Site, data = df.all) # Kruskal-Wallis chi-squared = 53.295, df = 3, p-value = 1.586e-11
pairwise.wilcox.test(df.all$Snps, df.all$Site,
                     p.adjust.method = "fdr")

wilcox.test(df3$Snps[which(df3$Type=="Deleterious" & df3$Site=="AZ")],
            df3$Snps[which(df3$Type=="Deleterious" & df3$Site=="WTX")]) # W = 868, p-value = 4.687e-11

wilcox.test(df3$Snps[which(df3$Type=="Tolerated" & df3$Site=="AZ")],
            df3$Snps[which(df3$Type=="Tolerated" & df3$Site=="WTX")]) # W = 868,p-value < 2.2e-16

wilcox.test(df3$Snps[which(df3$Type=="Nonsynonymous" & df3$Site=="AZ")],
            df3$Snps[which(df3$Type=="Nonsynonymous" & df3$Site=="WTX")]) # W = 868,p-value < 2.2e-16

# Masked Load
ggplot(df,aes(x=Type,y=Het,fill=Site))+
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  geom_point(position=position_jitterdodge(), size=2,aes(color=Site),color="black",pch=21)+
  labs(x="Site", y = "No. of heterozygous SNPs")+
  theme_classic(base_size = 15)+
  #scale_x_discrete(labels=c("Deleterious" = "Deleterious", "Tolerated" = "Tolerated",
  #                          "Nonsynonymous" = "Functional"))+
  scale_fill_manual(values=cols[1:4])+
  scale_color_manual(values=cols[1:4])+
  scale_y_continuous(labels = scales::comma)+
  theme(axis.title.y=element_text(size=15),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none")

df.all <- df[which(df$Type=="Nonsynonymous"),]
kruskal.test(Het ~ Site, data = df.all) # Kruskal-Wallis chi-squared = 53.295, df = 3, p-value = 1.586e-11
pairwise.wilcox.test(df.all$Het, df.all$Site,
                     p.adjust.method = "fdr")

wilcox.test(df3$Het[which(df3$Type=="Deleterious" & df3$Site=="AZ")],
            df3$Het[which(df3$Type=="Deleterious" & df3$Site=="WTX")],conf.int = T) # W = 180, p-value = 0.000119

wilcox.test(df3$Het[which(df3$Type=="Tolerated" & df3$Site=="AZ")],
            df3$Het[which(df3$Type=="Tolerated" & df3$Site=="WTX")],conf.int = T) # W = 354.5, p-value = 0.2304

wilcox.test(df3$Het[which(df3$Type=="Nonsynonymous" & df3$Site=="AZ")],
            df3$Het[which(df3$Type=="Nonsynonymous" & df3$Site=="WTX")],conf.int = T) # W = 601, p-value = 0.01075

# PropHomo
ggplot(df,aes(x=Type,y=Althom,fill=Site))+
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  geom_point(position=position_jitterdodge(), size=2,aes(color=Site),color="black",pch=21)+
  labs(x="Site", y = "No. of derived homozyogous SNPs")+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:4])+
  scale_color_manual(values=cols[1:4])+
  scale_y_continuous(labels = scales::comma)+
  #scale_x_discrete(labels=c("Deleterious" = "Deleterious", "Tolerated" = "Tolerated",
  #                          "Nonsynonymous" = "Functional"))+
  theme(axis.title.y=element_text(size=15),
        axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position = "none")

kruskal.test(Althom ~ Site, data = df.all) # Kruskal-Wallis chi-squared = 53.295, df = 3, p-value = 1.586e-11
pairwise.wilcox.test(df.all$Althom, df.all$Site,
                     p.adjust.method = "fdr")

wilcox.test(df3$Althom[which(df3$Type=="Deleterious" & df3$Site=="AZ")],correct = T,
            df3$Althom[which(df3$Type=="Deleterious" & df3$Site=="WTX")],conf.int = T) # W = 12, p-value = 1.555e-10

wilcox.test(df3$Althom[which(df3$Type=="Tolerated" & df3$Site=="AZ")],correct = T,
            df3$Althom[which(df3$Type=="Tolerated" & df3$Site=="WTX")],conf.int = T) # W = 28, p-value = 7.466e-10

wilcox.test(df3$Althom[which(df3$Type=="Nonsynonymous" & df3$Site=="AZ")],correct = T,
            df3$Althom[which(df3$Type=="Nonsynonymous" & df3$Site=="WTX")],conf.int = T) # W = 28, p-value = 7.466e-10

# Calculate potential and realized load

# Potential load = # no. of deleterious SNPs / Total no. of functional SNPs
# Realized load = # no. of alt homo sites / deleterious SNPs

df1 <- cbind(df,df$Snps/df$Totalsnps)
df1 <- cbind(df1,df$Althom/(df$Snps))
colnames(df1)[10] <- "LoadP"
colnames(df1)[11] <- "LoadR"


df2 <- df1[-which(df1$Site == "CTX" | df1$Site == "MX"),]
df2 <- df2[-which(df2$Type == "Nonsynonymous"),]

ggplot(df2,aes(x=Type,y=LoadP,fill=Site))+
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  geom_point(position=position_jitterdodge(), size=2,aes(color=Site),color="black",pch=21)+
  labs(x="Site", y = expression(paste("Potential Load (",Load[P],")")))+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:4])+
  scale_color_manual(values=cols[1:4])+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position = "none")

df.tol <- df1[which(df1$Type=="Tolerated"),]
ddply(df.tol,"Site",summarise,Mean=mean(LoadR),SD=sd(LoadR))
kruskal.test(LoadR ~ Site, data = df.tol) # Kruskal-Wallis chi-squared = 53.295, df = 3, p-value = 1.586e-11
pairwise.wilcox.test(df.tol$LoadR, df.tol$Site,
                     p.adjust.method = "fdr")

wilcox.test(df2$LoadP[which(df2$Type=="Deleterious" & df2$Site=="AZ")],
            df2$LoadP[which(df2$Type=="Deleterious" & df2$Site=="WTX")],conf.int = T) #W = 356, p-value = 0.2414

wilcox.test(df2$LoadP[which(df2$Type=="Tolerated" & df2$Site=="AZ")],
            df2$LoadP[which(df2$Type=="Tolerated" & df2$Site=="WTX")],conf.int = T) # W = 868, p-value < 2.2e-16

# Mean and SD

mean(df2$LoadP[which(df2$Site=="AZ" & df2$Type=="Deleterious")]) # 0.07129103
mean(df2$LoadP[which(df2$Site=="WTX" & df2$Type=="Deleterious")]) # 0.08901395
sd(df2$LoadP[which(df2$Site=="AZ" & df2$Type=="Deleterious")]) # 0.0005090463
sd(df2$LoadP[which(df2$Site=="WTX" & df2$Type=="Deleterious")]) # 0.03690097

mean(df2$LoadP[which(df2$Site=="AZ" & df2$Type=="Tolerated")]) # 0.1539355
mean(df2$LoadP[which(df2$Site=="WTX" & df2$Type=="Tolerated" & df2$LoadP < 0.2)]) # 0.1301437
sd(df2$LoadP[which(df2$Site=="AZ" & df2$Type=="Tolerated")]) # 0.004947425
sd(df2$LoadP[which(df2$Site=="WTX" & df2$Type=="Tolerated")]) # 0.007359706

ggplot(df2,aes(x=Type,y=LoadR,fill=Site))+
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  geom_point(position=position_jitterdodge(), size=2,aes(color=Site),color="black",pch=21)+
  labs(x="Site", y = expression(paste("Realized Load (",Load[R],")")))+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:4])+
  scale_color_manual(values=cols[1:4])+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position = "none")
  
wilcox.test(df2$LoadR[which(df2$Type=="Deleterious" & df2$Site=="AZ")],
            df2$LoadR[which(df2$Type=="Deleterious" & df2$Site=="WTX")],conf.int = T) # W = 2, p-value < 2.2e-16

wilcox.test(df2$LoadR[which(df2$Type=="Tolerated" & df2$Site=="AZ")],
            df2$LoadR[which(df2$Type=="Tolerated" & df2$Site=="WTX")],conf.int = T) #  W = 4, p-value < 2.2e-16

mean(df2$LoadR[which(df2$Site=="AZ" & df2$Type=="Deleterious")]) # 0.0109189
mean(df2$LoadR[which(df2$Site=="WTX" & df2$Type=="Deleterious")]) # 0.01565053
sd(df2$LoadR[which(df2$Site=="AZ" & df2$Type=="Deleterious")]) # 0.0009535385
sd(df2$LoadR[which(df2$Site=="WTX" & df2$Type=="Deleterious")]) # 0.001753375

mean(df2$LoadR[which(df2$Site=="AZ" & df2$Type=="Tolerated")]) # 0.0186585
mean(df2$LoadR[which(df2$Site=="WTX" & df2$Type=="Tolerated")]) # 0.02569688
sd(df2$LoadR[which(df2$Site=="AZ" & df2$Type=="Tolerated")]) # 0.001424254
sd(df2$LoadR[which(df2$Site=="WTX" & df2$Type=="Tolerated")]) # 0.001297107

# Compare load with ROH
roh <- read.table("~/Documents/Thesis_Research/Final Results/Ch3/revise/rohs/ROH_individual.txt", header = T)
roh <- roh[-c(29:32,64:66),]
del <- df2[which(df2$Type=="Deleterious"),]
tol <- df2[which(df2$Type=="Tolerated"),]

df.del <- cbind(del,roh)
df.del <- df.del[,c(1,2,10,11,14,15,16,17)]

df.tol <- cbind(tol,roh)
df.tol <- df.tol[,c(1,2,10,11,14,15,16,17)]

df.roh <- rbind(df.del,df.tol)

# Compare it with Het
### Heterozygosity ###

het_all <- read.table("~/Documents/Thesis_Research/Final Results/Ch3/revise/diversity/best66.old.auto.noSing.nomiss.het", header = T)
total <- 960796788 # Total genomic length analyzed

het_all <- cbind(het_all,(het_all$N_SITES-het_all$O.HOM.)/total)
colnames(het_all)[6] <- "het"

az_het <- data.frame(Site="AZ", value = het_all$het[1:28])
wtx_het <- data.frame(Site="WTX", value = het_all$het[33:63])

het <- rbind(az_het,wtx_het)
het$Site <- factor(het$Site , levels=c("AZ", "WTX"))

df.del <- cbind(df2,rbind(het,het))
colnames(df.del)[13] <- "Het"

write.table(df.del,"DiversityLoadTable.txt",quote = F,row.names = F)


# Plot relationships

diversity <- read.table("DiversityLoadTable.txt",header = T)

div_del <- diversity[which(diversity$Type=="Deleterious"),]
div_tol <- diversity[which(diversity$Type=="Tolerated"),]

az.divDel <- div_del[which(div_del$Site=="AZ"),]
wtx.divDel <- div_del[which(div_del$Site=="WTX"),]

az.divTol <- div_tol[which(div_tol$Site=="AZ"),]
wtx.divTol <- div_tol[which(div_tol$Site=="WTX"),]

p1 <- ggplot(div_del,aes(x=Het.1,y=LoadP,fill=Site, color=Site))+
  geom_point(size=5, shape=16)+
  labs(x="Genome-wide Heterozygosity", y = expression(paste(Load[P])))+
  geom_smooth(method=lm, color=col_az, fill=adjustcolor(col_az,alpha.f = 0.1), se=T, aes(x=Het.1,y=LoadP),data=az.divDel)+
  geom_smooth(method=lm, color=col_wtx, fill=adjustcolor(col_wtx,alpha.f = 0.1), se=T, aes(x=Het.1,y=LoadP),data=wtx.divDel)+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:2])+
  scale_color_manual(values=cols[1:2])+
  theme(legend.title=element_blank(),
        legend.position = "none")

p2 <- ggplot(div_tol,aes(x=Het.1,y=LoadP,fill=Site, color=Site,shape=Type))+
  geom_point(size=5,shape=17)+
  labs(x="Genome-wide Heterozygosity", y = expression(paste(Load[P])))+
  geom_smooth(method=lm, color=col_az, fill=adjustcolor(col_az,alpha.f = 0.1), se=T, aes(x=Het.1,y=LoadP),data=az.divTol)+
  geom_smooth(method=lm, color=col_wtx, fill=adjustcolor(col_wtx,alpha.f = 0.1), se=T, aes(x=Het.1,y=LoadP),data=wtx.divTol)+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:2])+
  scale_color_manual(values=cols[1:2])+
  theme(legend.title=element_blank(),
        legend.position = "none")

p3 <- ggplot(div_del,aes(x=Het.1,y=LoadR,fill=Site, color=Site,shape=Type))+
  geom_point(size=5,shape=16)+
  labs(x="Genome-wide Heterozygosity", y = expression(paste(Load[R])))+
  geom_smooth(method=lm, color=col_az, fill=adjustcolor(col_az,alpha.f = 0.1), se=T, aes(x=Het.1,y=LoadR),data=az.divDel)+
  geom_smooth(method=lm, color=col_wtx, fill=adjustcolor(col_wtx,alpha.f = 0.1), se=T, aes(x=Het.1,y=LoadR),data=wtx.divDel)+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:2])+
  scale_color_manual(values=cols[1:2])+
  theme(legend.title=element_blank(),
        legend.position = "none")

p4 <- ggplot(div_tol,aes(x=Het.1,y=LoadR,fill=Site, color=Site,shape=Type))+
  geom_point(size=5,shape=17)+
  labs(x="Genome-wide Heterozygosity", y = expression(paste(Load[R])))+
  geom_smooth(method=lm, color=col_az, fill=adjustcolor(col_az,alpha.f = 0.1), se=T, aes(x=Het.1,y=LoadR),data=az.divTol)+
  geom_smooth(method=lm, color=col_wtx, fill=adjustcolor(col_wtx,alpha.f = 0.1), se=T, aes(x=Het.1,y=LoadR),data=wtx.divTol)+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:2])+
  scale_color_manual(values=cols[1:2])+
  theme(legend.title=element_blank(),
        legend.position = "none")

grid.arrange(p1,p2,p3,p4, nrow = 2)

summary(lm(az.divDel$LoadP ~ az.divDel$Het.1)) # Multiple R-squared:  0.01102;  p-value: 0.5949
summary(lm(wtx.divDel$LoadP ~ wtx.divDel$Het.1)) # Multiple R-squared:  0.1345;  p-value: 0.04239

summary(lm(az.divTol$LoadP ~ az.divTol$Het.1)) # Multiple R-squared:  0.2696;  p-value: 0.002729
summary(lm(wtx.divTol$LoadP ~ wtx.divTol$Het.1)) # Multiple R-squared:  0.01269;  p-value: 0.5463

summary(lm(az.divDel$LoadR ~ az.divDel$Het.1)) # Multiple R-squared:  0.1981;  p-value: 0.01022
summary(lm(wtx.divDel$LoadR ~ wtx.divDel$Het.1)) # Multiple R-squared:  0.07143;  p-value: 0.1461

summary(lm(az.divTol$LoadR ~ az.divTol$Het.1)) # Multiple R-squared:  0.01552 ;  p-value: 0.5276
summary(lm(wtx.divTol$LoadR ~ wtx.divTol$Het.1)) # Multiple R-squared:  0.04534;  p-value: 0.1303

summary(lm(wtx.div$LoadP[which(wtx.div$Type=="Deleterious")] ~
             wtx.div$Het[which(wtx.div$Type=="Deleterious")])) # Multiple R-squared:  0.01182;  p-value: 0.5604
summary(lm(wtx.div$LoadP[which(wtx.div$Type=="Tolerated")] ~
             wtx.div$Het[which(wtx.div$Type=="Tolerated")])) # Multiple R-squared:  0.3359;  p-value: 0.0003768

summary(lm(az.div$LoadR[which(az.div$Type=="Deleterious")] ~
             az.div$Het[which(az.div$Type=="Deleterious")])) # Multiple R-squared:  0.2899;  p-value: 0.003121
summary(lm(az.div$LoadR[which(az.div$Type=="Tolerated")] ~
             az.div$Het[which(az.div$Type=="Tolerated")])) # Multiple R-squared:  0.4106;  p-value: 0.0002391
summary(lm(wtx.div$LoadR[which(wtx.div$Type=="Deleterious")] ~
             wtx.div$Het[which(wtx.div$Type=="Deleterious")])) # Multiple R-squared:  0.4388;  p-value: 2.956e-05
summary(lm(wtx.div$LoadR[which(wtx.div$Type=="Tolerated")] ~
             wtx.div$Het[which(wtx.div$Type=="Tolerated")])) # Multiple R-squared:  0.62;  p-value: 1.475e-07




##### Load of different ROHs ####

roh.df <- read.table("ROhLoad.txt",header = T)
roh.df <- roh.df[-which(roh.df$Site == "CTX" | roh.df$Site =="MX"),]

del <- roh.df[which(roh.df$Type=="Deleterious"),]
tol <- roh.df[which(roh.df$Type=="Tolerated"),]
fun <- roh.df[which(roh.df$Type=="Functional"),]


# Get total alleles

tot.all <- datafile$All[which(datafile$Type=="Deleterious")]
tot.all <- tot.all[c(1:28,33:63)]

tot.alt <- datafile$Althom[which(datafile$Type=="Deleterious")]
tot.alt <- tot.alt[c(1:28,33:63)]

allAll <- NULL
allAlt <- NULL
for (i in c(1:59))
{
  allAll <- c(allAll,rep(tot.all[i],3))
  allAlt <- c(allAlt,rep(tot.alt[i],3))
}
allAll <- allAll[-58]
allAlt <- allAlt[-58]

del <- cbind(del,allAlt,allAll,del$All/allAlt,del$All/allAll)
colnames(del)[11:12] <- c("AltPer","AllPer")


ggplot(del,aes(x=SampleID,y=AltPer,fill=Length))+
  geom_bar(position="stack", stat="identity", color="black")+
  #geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  labs(x="Site", y = expression(paste("% of alternate homozygous SNPs")))+
  theme_classic(base_size = 22)+
  scale_fill_manual(values=c("#F8766D","#00BA38","#619CFF"))+
  #scale_fill_manual(values=cols[1:3])+
  #scale_color_manual(values=cols[1:3])+
  theme(axis.title.y=element_blank(),
        legend.title=element_blank())
# legend.position = "none")

df1 <- inner_join(del,tol,by=c("SampleID","Site","Length"))
df1 <- inner_join(df1,fun,by=c("SampleID","Site","Length"))

df2 <- cbind(df1[,c(1,2,3)],
             df1$All.x/df1$Snps.x,df1$All.y/df1$Snps.y)

colnames(df2)[4:5] <- c("LoadR_Del","LoadR_Tol")

df3 <- rbind(data.frame(SampleID=df2$SampleID,Site=df2$Site,Length=df2$Length,
                        Type="Deleterious",LoadR=df2$LoadR_Del),
             data.frame(SampleID=df2$SampleID,Site=df2$Site,Length=df2$Length,
                        Type="Tolerated",LoadR=df2$LoadR_Tol))

df3$Site <- factor(df3$Site , levels=c("AZ", "WTX"))
df3$Length <- factor(df3$Length , levels=c("Long", "Medium","Short"))


#df3 <- df3[-which(df3$LoadP == 0 |df3$LoadR ==0 ),]

del.df3 <- df3[which(df3$Type=="Deleterious"),]
tol.df3 <- df3[which(df3$Type=="Tolerated"),]

id <- NULL
for (i in c(1:27,29:59,28))
{
  id <- c(id,rep(i,3))
}
id <- id[-length(id)]

del.df3$SampleID <- id
tol.df3$SampleID <- id
#del.df3$SampleID <-  factor(del.df3$SampleID , levels=as.character(id))

my_comparisons <- list( c("Long", "Medium"), c("Long", "Short"), c("Medium", "Short"))
my_comparisons <- list( c("AZ", "WTX"))

del.df3 <- del.df3[-which(del.df3$LoadR==1),]
del.df3 <- del.df3[-which(is.na(del.df3$LoadR)),]
tol.df3 <- tol.df3[-which(is.na(tol.df3$LoadR)),]

p1 <- ggplot(del.df3,aes(x=Length,y=LoadR,fill=Site))+ # pdf dim = 3x12.7
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  theme_classic(base_size = 15)+
  #scale_fill_manual(values=c("#F8766D","#00BA38","#619CFF"))+
  scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.1))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",method = "wilcox",
                     bracket.size=1, size=5)+
  geom_point(position=position_jitterdodge(), size=2,aes(color=Length),color="black",pch=21)+
  labs(x="Site", y = expression(paste(Load[R]," within ROHs (Deleterious)")))+
  scale_fill_manual(values=cols[1:3])+
  scale_color_manual(values=cols[1:3])+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position = "none")


wilcox.test(del.df3$LoadR[which(del.df3$Length=="Short" & del.df3$Site == "AZ")],
            del.df3$LoadR[which(del.df3$Length=="Short" & del.df3$Site == "WTX")]) # W = 307, p-value = 0.05439
wilcox.test(del.df3$LoadR[which(del.df3$Length=="Medium" & del.df3$Site == "AZ")],
            del.df3$LoadR[which(del.df3$Length=="Medium" & del.df3$Site == "WTX")]) # W = 156, p-value = 0.0001398
wilcox.test(del.df3$LoadR[which(del.df3$Length=="Long" & del.df3$Site == "AZ")],
            del.df3$LoadR[which(del.df3$Length=="Long" & del.df3$Site == "WTX")]) # W = 334, p-value = 0.5234


p2 <- ggplot(tol.df3,aes(x=Length,y=LoadR,fill=Site))+ # pdf dim = 3x12.7
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  theme_classic(base_size = 15)+
  #scale_fill_manual(values=c("#F8766D","#00BA38","#619CFF"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",method = "wilcox",
                     bracket.size=1, size=5)+
  scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.1))+
  geom_point(position=position_jitterdodge(), size=2,aes(color=Length),color="black",pch=21)+
  labs(x="Site", y = expression(paste(Load[R]," within ROHs (Tolerated)")))+
  scale_fill_manual(values=cols[1:3])+
  scale_color_manual(values=cols[1:3])+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position = "none")
grid.arrange(p1,p2, nrow = 1)

wilcox.test(tol.df3$LoadR[which(tol.df3$Length=="Short" & tol.df3$Site == "AZ")],
            tol.df3$LoadR[which(tol.df3$Length=="Short" & tol.df3$Site == "WTX")]) # W = 310.5, p-value = 0.06189
wilcox.test(tol.df3$LoadR[which(tol.df3$Length=="Medium" & tol.df3$Site == "AZ")],
            tol.df3$LoadR[which(tol.df3$Length=="Medium" & tol.df3$Site == "WTX")]) #W = 275, p-value = 0.06488
wilcox.test(tol.df3$LoadR[which(tol.df3$Length=="Long" & tol.df3$Site == "AZ")],
            tol.df3$LoadR[which(tol.df3$Length=="Long" & tol.df3$Site == "WTX")]) # W = 369.5, p-value = 0.972


p3 <- ggplot(del.df3,aes(x=Length,y=LoadR,fill=Length))+ # pdf dim = 3x12.7
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=c("#F8766D","#00BA38","#619CFF"))+
  scale_color_manual(values=c("#F8766D","#00BA38","#619CFF"))+
  scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.1))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",method = "wilcox",
                     bracket.size=1, size=5)+
  geom_point(position=position_jitterdodge(), size=2,aes(color=Length),color="black",pch=21)+
  labs(x="Site", y = expression(paste(Load[R]," within ROHs (Deleterious)")))+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position = "none")

p4 <- ggplot(tol.df3,aes(x=Length,y=LoadR,fill=Length))+ # pdf dim = 3x12.7
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=c("#F8766D","#00BA38","#619CFF"))+
  scale_color_manual(values=c("#F8766D","#00BA38","#619CFF"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",method = "wilcox",
                     bracket.size=1, size=5)+
  scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.1))+
  geom_point(position=position_jitterdodge(), size=2,aes(color=Length),color="black",pch=21)+
  labs(x="Site", y = expression(paste(Load[R]," within ROHs (Tolerated)")))+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position = "none")

grid.arrange(p3,p4, nrow = 1)

kruskal.test(LoadR ~ Length, data = del.df3) # Kruskal-Wallis chi-squared = 20.149, df = 3, p-value = 0.0001581
pairwise.wilcox.test(del.df3$LoadR, del.df3$Length,
                     p.adjust.method = "fdr")

kruskal.test(LoadR ~ Length, data = tol.df3) # Kruskal-Wallis chi-squared = 20.149, df = 3, p-value = 0.0001581
pairwise.wilcox.test(tol.df3$LoadR, tol.df3$Length,
                     p.adjust.method = "fdr")

Type <- c(rep("Deleterious",nrow(del.df3)),rep("Tolerated",nrow(del.df3)))
df <- rbind(del.df3,tol.df3)
#df <- cbind(df,Type)

ggplot(df,aes(x=Length,y=LoadR,fill=Type))+ # pdf dim = 3x12.7
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  theme_classic(base_size = 15)+
  #scale_fill_manual(values=c("#F8766D","#00BA38","#619CFF"))+
  #scale_color_manual(values=c("#F8766D","#00BA38","#619CFF"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",method = "wilcox",
                     bracket.size=1, size=5)+
  scale_y_continuous(limits=c(0, 1), breaks = seq(0,1,0.1))+
  geom_point(position=position_jitterdodge(), size=2,aes(color=Length),color="black",pch=21)+
  labs(y = expression(paste(Load[R]," within ROHs")))+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position = "none")

kruskal.test(LoadR ~ Length, data = df) # Kruskal-Wallis chi-squared = 20.149, df = 3, p-value = 0.0001581
pairwise.wilcox.test(df$LoadR, df$Length,
                     p.adjust.method = "fdr")

wilcox.test(df$LoadR[which(df$Length=="Short" & df$Type == "Deleterious")],
            df$LoadR[which(df$Length=="Short" & df$Type == "Tolerated")]) # W = 1259.5, p-value = 0.009704
wilcox.test(df$LoadR[which(df$Length=="Medium" & df$Type == "Deleterious")],
            df$LoadR[which(df$Length=="Medium" & df$Type == "Tolerated")]) # W = 1015.5, p-value = 0.001316
wilcox.test(df$LoadR[which(df$Length=="Long" & df$Type == "Deleterious")],
            df$LoadR[which(df$Length=="Long" & df$Type == "Tolerated")]) # W = 2188, p-value = 4.473e-05
