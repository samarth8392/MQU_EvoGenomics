###########################################################################
###                          Samarth Mathur, PhD                     	  ###
###                        The Ohio State University                 	  ###
###                                                                     ###
###########################################################################
###########################################################################
###                   ROHs.R         		                                ###
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

# Load colors
#display_carto_all(colorblind_friendly = TRUE)

col_az <- carto_pal(12,"Vivid")[11]
col_wtx <- carto_pal(7,"Purp")[5]
col_etx <-carto_pal(7,"Peach")[5]
col_mx <- "#3598FF"
cols <- c(col_az,col_wtx,col_etx, col_mx)
colss <-c(rep(col_az,28),rep(col_mx,4),rep(col_wtx,31),rep(col_etx,3))
cols_l <- adjustcolor(cols, alpha.f = 0.5)
col_azwtx <- "#90CCC2"
#### PREREQUISITES: END #####

setwd("~/Documents/Thesis_Research/Final Results/Ch3/revise/rohs/")

roh <- read.table("best66.old.auto.poly.nomiss.hom.indiv", header = T)
a <- read.table("onlyAZ.recode.nomiss.hom.indiv", header = T)
b <- read.table("onlyWTX.recode.nomiss.hom.indiv", header = T)
total <- 960796.788 # Total genomic length analyzed

# Plot FROH
az <- data.frame(Site="AZ", froh = a$KB/total,nroh=a$NSEG,lroh=a$KB,avgkb=a$KBAVG)
wtx <- data.frame(Site="WTX", froh = b$KB/total,nroh=b$NSEG,lroh=b$KB,avgkb=b$KBAVG)
etx <- data.frame(Site="CTX", froh = roh$KB[64:66]/total,nroh=roh$NSEG[64:66],lroh=roh$KB[64:66],avgkb=roh$KBAVG[64:66])
mx <- data.frame(Site="MX", froh = roh$KB[29:32]/total,nroh=roh$NSEG[29:32],lroh=roh$KB[29:32],avgkb=roh$KBAVG[29:32])

SampleID <- c(roh$FID[1:28],roh$FID[29:32],roh$FID[33:63],roh$FID[64:66])

theta_all <- rbind(az,mx,wtx,etx)
theta_all <- cbind(SampleID,theta_all)
theta_all$Site <- factor(theta_all$Site , levels=c("AZ", "WTX","CTX","MX"))

write.table(theta_all,"ROH_individual.txt",quote=F,row.names = F)

kruskal.test(froh ~ Site, data = theta_all) # Kruskal-Wallis chi-squared = 18.399, df = 3, p-value = 0.0003639
pairwise.wilcox.test(theta_all$froh, theta_all$Site,
                     p.adjust.method = "fdr")

mu_all <- ddply(theta_all, "Site", summarise, grp.mean=mean(froh))
sd_all <- ddply(theta_all, "Site", summarise, grp.mean=sd(froh))

my_comparisons <- list( c("AZ", "WTX"), c("AZ", "CTX"), c("AZ", "MX"),
                        c("WTX", "CTX"),c("WTX", "MX"), c("CTX", "MX"))

theta_all$Site <- factor(theta_all$Site , levels=c("AZ", "WTX", "CTX", "MX"))

p1<- ggplot(theta_all, aes(x=Site, y=froh, fill=Site))+
  geom_boxplot(lwd=1)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",method = "wilcox",
                     bracket.size=1, size=5)+ 
  theme_classic(base_size = 15)+ scale_fill_manual(values=cols)+
  labs(x="Site", y =expression(paste(F[ROH])))+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position = "none")
p1

# ROH vs coverage

cov <- read.table("~/Documents/Thesis_Research/Final Results/Ch3/revise/Final_Figs/SI/DataS1.txt",header = T)
het <- read.table("~/Documents/Thesis_Research/Final Results/Ch3/revise/diversity/best66.auto.noSing.nomiss.het", header = T)

total <- 960796788 # Total genomic length analyzed

het <- cbind(het,(het$N_SITES-het$O.HOM.)/total)
colnames(het)[c(1,6)] <- c("SampleID","het")
df <- inner_join(theta_all,cov,by="SampleID")
df <- inner_join(df,het,by="SampleID")

df$froh[1:28] <- sample(theta_all$froh[which(theta_all$Site=="AZ")],28)
df$froh[33:63] <- sample(theta_all$froh[which(theta_all$Site=="WTX")],31)
#df <- cbind(theta_all,cov)
#df <- cbind(df,het$value)
#colnames(df)[9] <- "hets"

p1 <- ggplot(df,aes(x=het,y=froh, fill=Site))+
  geom_point(size=5,  shape=21)+
  geom_smooth(method=lm , aes(x=het, y=froh), fill=adjustcolor("#69b3a2",alpha.f = 0.1), se=T) +
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:4])+
  labs(y=expression(paste(F[ROH])), x="Genome-wide heterozygosity")

p2 <- ggplot(df,aes(x=MeanDepth,y=het, fill=Site))+
  geom_point(size=5,  shape=21)+
  geom_smooth(method=lm , aes(x=MeanDepth, y=het), fill=adjustcolor("#69b3a2",alpha.f = 0.1), se=T) +
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:4])+
  labs(x="Mean depth of coverage (X)", y="Genome-wide heterozygosity")

p3 <- ggplot(df,aes(x=MeanDepth,y=froh, fill=Site))+
  geom_point(size=5,  shape=21)+
  geom_smooth(method=lm , aes(x=MeanDepth, y=froh), fill=adjustcolor("#69b3a2",alpha.f = 0.1), se=T) +
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:4])+
  labs(y=expression(paste(F[ROH])), x="Mean depth of coverage (X)")

grid.arrange(p1,p2, p3, nrow = 1)

summary(lm(df$froh ~ df$het)) # Adjusted R-squared:   0.143 ; p-value: 0.001021
summary(lm(df$het ~ df$MeanDepth)) # Adjusted R-squared:  0.01453 ; p-value: 0.1665
summary(lm(df$froh ~ df$MeanDepth)) # Multiple R-squared:  0.1027; p-value: 0.008724

mean(df$froh[which(df$Site=="AZ")])

# Average length of ROH

ggplot(theta_all,aes(fill=Site,x=Site,y=avgkb))+
  geom_boxplot(alpha = 0.5,outlier.shape = NA)+
  geom_point(position=position_jitterdodge(), size=2,aes(color=group),color="black",pch=21)+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:4])+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",method = "wilcox",
                     bracket.size=1, size=5)+ 
  labs(y=expression("Average Length of ROHs (Kb)"))+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position = "none")

### Plot ROHs ###

hom <- readExternalRuns("best66.old.auto.poly.nomiss.hom", program = "plink")
#group
a <- noquote(unique(hom$group))
r <- data.frame(from = c(a[1:28],a[29:32],a[33:63],a[64:66]),
                to=c(rep("AZ",28),rep("MX",4),rep("WTX",31),rep("CTX",3)))
hom_new <- FindReplace(data = hom, Var = "group", replaceData = r,
                       from = "from", to = "to", exact = T)

hom_new <- hom_new[-which(hom_new$group == "AZ" | hom_new$group == "WTX"),]

az <- readExternalRuns("onlyAZ.recode.nomiss.hom", program = "plink")
wtx <- readExternalRuns("onlyWTX.recode.nomiss.hom", program = "plink")
hom1 <- rbind(az,wtx)
hom_new2 <- FindReplace(data = hom1, Var = "group", replaceData = r,
                       from = "from", to = "to", exact = T)

hom_new3 <- rbind(hom_new2,hom_new)


hom_new3$lengthBps <- hom_new3$lengthBps/1000000

# Divide by roh length
# Short - 200-300kb
# Medium - 300-500kb
# Long - > 500kb
short <- hom_new3[which(hom_new3$lengthBps >= 0.2 & hom_new3$lengthBps < 0.3),]
medium <- hom_new3[which(hom_new3$lengthBps >= 0.3 & hom_new3$lengthBps < 0.6),]
long <- hom_new3[which(hom_new3$lengthBps >= 0.6),]

s <- ddply(short,"id",summarise,Sum=sum(lengthBps))
m <- ddply(medium,"id",summarise,Sum=sum(lengthBps))
l <- ddply(long,"id",summarise,Sum=sum(lengthBps))

df <- rbind(data.frame(SampleID=s$id,Type="Short",Length=s$Sum),
            data.frame(SampleID=m$id,Type="Medium",Length=m$Sum),
            data.frame(SampleID=l$id,Type="Long",Length=l$Sum))

Site <- c(rep(c(rep("AZ",28),rep("MX",4),rep("WTX",31),rep("CTX",3)),2),
          c(rep("AZ",27),rep("MX",4),rep("WTX",31),rep("CTX",3)))
          

df <- cbind(Site,df)
df$Site <-  factor(df$Site , levels=c("AZ", "WTX", "CTX", "MX"))
df$SampleID <- c(1:66,1:66,1:27,29:66)
df$SampleID <-  factor(df$SampleID , levels=c(1:66))
ggplot(df, aes(fill=Type, y=Length, x=SampleID)) + 
  geom_bar(position="stack", stat="identity", color="black")+
  theme_classic(base_size = 22)+ 
  labs(x="Site", y ="Length (Mb)")+
  #scale_x_continuous(limits=c(1, 66), breaks = seq(1,66,1))+
  #scale_fill_manual(values=cols_l)+
  #scale_color_manual(values=cols)+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        #axis.text.x =element_blank(),
        axis.text.x =element_text(color=colss,size=13),
        legend.position = "none")

# Extract ROHs for each individual

for (len in c("short","medium","long"))
{
  for (ind in s$id)
  {
    a <- get(paste(len))
    df <- a[which(a$id == ind),]
    pos <- df[,c(3,5,6)]
    colnames(pos) <- c("CHROM","START","END")
    write.table(pos,paste(ind,".",len,"ROHs.bed",sep=""),quote=F,row.names=F)
  }
}


