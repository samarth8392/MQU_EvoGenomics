###########################################################################
###                          Samarth Mathur, PhD                     	  ###
###                        The Ohio State University                 	  ###
###                                                                     ###
###########################################################################
###########################################################################
###                   Rxy.R          		                                ###
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

# Load colors
col_az <- carto_pal(12,"Vivid")[11]
col_wtx <- carto_pal(7,"Purp")[5]
col_etx <-carto_pal(7,"Peach")[5]
col_mx <- "#3598FF"
colss <-c(rep(col_az,28),rep(col_mx,4),rep(col_wtx,31),rep(col_etx,3))
cols_l <- adjustcolor(cols, alpha.f = 0.5)
col_azwtx <- "#90CCC2"
cols <- c(col_az,col_wtx,col_etx, col_mx,col_azwtx)
#### PREREQUISITES: END #####

setwd("~/Documents/Thesis_Research/Final Results/Ch3/revise/load/")

#To calculate Rxy statistics

# From Xue et al. (2015) 10.1126/science.aaa4484

# We began by calculating a statistic which
#compares two populations, given a particular category of sites, in terms of the number of
#derived alleles found at sites within that category in one population rather than the other.
#Specifically, at each site i we write the observed derived allele frequency in population A and Population B,
#Then if C is a particular category of protein-coding sites and I a set of intergenic sites, we define
# L(A,B) = sum(fiA(1-fiB))/sum(fjA(1-fjB))
#where i = deleterious site and j = neutral sites (synonymous sites)

#R(A,B) = L(A,B)/L(B/A)

AZ_del <- read.table("freq/az.best66.deleterious.freq",header = T)
WTX_del <- read.table("freq/wtx.best66.deleterious.freq",header = T)
ETX_del <- read.table("freq/etx.best66.deleterious.freq",header = T)
MX_del <- read.table("freq/mx.best66.deleterious.freq",header = T)

AZ_wdel <- read.table("freq/az.best66.tolerated.freq",header = T)
WTX_wdel <- read.table("freq/wtx.best66.tolerated.freq",header = T)
ETX_wdel <- read.table("freq/etx.best66.tolerated.freq",header = T)
MX_wdel <- read.table("freq/mx.best66.tolerated.freq",header = T)

AZ_non <- read.table("freq/az.best66.nonsynonymous.freq",header = T)
WTX_non <- read.table("freq/wtx.best66.nonsynonymous.freq",header = T)
ETX_non <- read.table("freq/etx.best66.nonsynonymous.freq",header = T)
MX_non <- read.table("freq/mx.best66.nonsynonymous.freq",header = T)


AZ_syn <- read.table("freq/az.best66.synonymous.freq",header = T)
WTX_syn <- read.table("freq/wtx.best66.synonymous.freq",header = T)
ETX_syn <- read.table("freq/etx.best66.synonymous.freq",header = T)
MX_syn <- read.table("freq/mx.best66.synonymous.freq",header = T)

# Calculate Rxy and jacknife it for each population pair

pops <- c("AZ","WTX","ETX","MX")
chr <- read.table("chick_chr.list",header = F)

count1 <- 1
count2 <- 1

Rxy <- NULL
r_df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(r_df) <- c("PopPairs","Rdel","Rtol","Rfun")
for (i in c(1:4))
{
  count1 <- count1+1
  count2<- 1
  for (j in c(1:4))
  {
    count2 <- count2+1
    if (count2 > count1)
    {
      if (i != j)
      {
        pop1 <- pops[i]
        pop2 <- pops[j]
        del_name <- paste("del_",pop1,"_",pop2,sep="")
        
        pop1_df <- get(paste(pop1,"_del",sep = ""))
        pop2_df <- get(paste(pop2,"_del",sep = ""))
        
        df <- inner_join(pop1_df,pop2_df,by=c("CHROM","POS"))
        df <- df[,c(1,2,4,6)]
        colnames(df) <- c("Chrom","Pos","Pop1_altfreq","Pop2_altfreq")
        df <- df[-which(df$Pop1_altfreq==0 & df$Pop2_altfreq==0),] 
        assign(del_name,df)
        
        tol_name <- paste("tol_",pop1,"_",pop2,sep="")
        
        pop1_df <- get(paste(pop1,"_wdel",sep = ""))
        pop2_df <- get(paste(pop2,"_wdel",sep = ""))
        
        df <- inner_join(pop1_df,pop2_df,by=c("CHROM","POS"))
        df <- df[,c(1,2,4,6)]
        colnames(df) <- c("Chrom","Pos","Pop1_altfreq","Pop2_altfreq")
        df <- df[-which(df$Pop1_altfreq==0 & df$Pop2_altfreq==0),] # N = 36329
        assign(tol_name,df)
        
        fun_name <- paste("fun_",pop1,"_",pop2,sep="")
        
        pop1_df <- get(paste(pop1,"_non",sep = ""))
        pop2_df <- get(paste(pop2,"_non",sep = ""))
        
        df <- inner_join(pop1_df,pop2_df,by=c("CHROM","POS"))
        df <- df[,c(1,2,4,6)]
        colnames(df) <- c("Chrom","Pos","Pop1_altfreq","Pop2_altfreq")
        df <- df[-which(df$Pop1_altfreq==0 & df$Pop2_altfreq==0),] # N = 36329
        assign(fun_name,df)
        
        
        syn_name <- paste("syn_",pop1,"_",pop2,sep="")
        
        pop1_df <- get(paste(pop1,"_syn",sep = ""))
        pop2_df <- get(paste(pop2,"_syn",sep = ""))
        
        df <- inner_join(pop1_df,pop2_df,by=c("CHROM","POS"))
        df <- df[,c(1,2,4,6)]
        colnames(df) <- c("Chrom","Pos","Pop1_altfreq","Pop2_altfreq")
        df <- df[-which(df$Pop1_altfreq==0 & df$Pop2_altfreq==0),] # N = 206862
        assign(syn_name,df)
        
        #Calculate Rxy
        
        del <- get(del_name)
        tol <- get(tol_name)
        fun <- get(fun_name)
        syn <- get(syn_name)
        
        l_del_xy <- sum(del[,3]*(1-del[,4]))/sum(syn[,3]*(1-syn[,4]))
        l_del_yx <- sum(del[,4]*(1-del[,3]))/sum(syn[,4]*(1-syn[,3]))
        r_del_xy <- l_del_xy/l_del_yx # 1.10245
        
        l_tol_xy <- sum(tol[,3]*(1-tol[,4]))/sum(syn[,3]*(1-syn[,4]))
        l_tol_yx <- sum(tol[,4]*(1-tol[,3]))/sum(syn[,4]*(1-syn[,3]))
        r_tol_xy <- l_tol_xy/l_tol_yx # 1.075284
        
        l_fun_xy <- sum(fun[,3]*(1-fun[,4]))/sum(syn[,3]*(1-syn[,4]))
        l_fun_yx <- sum(fun[,4]*(1-fun[,3]))/sum(syn[,4]*(1-syn[,3]))
        r_fun_xy <- l_fun_xy/l_fun_yx # 1.075284
        
        # Divide by chromosome
        
        r_del <- vector()
        r_wdel <- vector()
        r_fun <- vector()
        for (k in 1:33)
        {
          a <- del[which(del$Chrom != chr$V1[k]),]
          b <- tol[which(tol$Chrom != chr$V1[k]),]
          d <- fun[which(fun$Chrom != chr$V1[k]),]
          c <- syn[which(syn$Chrom != chr$V1[k]),]
          l_del_1 <- sum(a[,3]*(1-a[,4]))/sum(c[,3]*(1-c[,4]))
          l_del_2 <- sum(a[,4]*(1-a[,3]))/sum(c[,4]*(1-c[,3]))
          
          l_wdel_1 <- sum(b[,3]*(1-b[,4]))/sum(c[,3]*(1-c[,4]))
          l_wdel_2 <- sum(b[,4]*(1-b[,3]))/sum(c[,4]*(1-c[,3]))
          
          l_fun_1 <- sum(d[,3]*(1-d[,4]))/sum(c[,3]*(1-c[,4]))
          l_fun_2 <- sum(d[,4]*(1-d[,3]))/sum(c[,4]*(1-c[,3]))
          
          r_del[k] <- l_del_1/l_del_2
          r_wdel[k] <- l_wdel_1/l_wdel_2
          r_fun[k] <- l_fun_1/l_fun_2
        }
        r_df <- data.frame(PopPairs=paste(pop1,"-",pop2,sep=""),Rdel=r_del,Rtol=r_wdel,Rfun=r_fun)
      }
      Rxy <- rbind(Rxy,r_df)
    }
  }
}

# Plot Ryx = 1/Rxy
Rxy$PopPairs <- factor(Rxy$PopPairs, levels = c("AZ-WTX","AZ-ETX","AZ-MX","WTX-ETX","WTX-MX","ETX-MX"))

Ryx <- as.data.frame(cbind(1/Rxy$Rdel,1/Rxy$Rtol,1/Rxy$Rfun))
rev.pops <- c(rep("WTX-AZ",33),rep("ETX-AZ",33),rep("MX-AZ",33),rep("ETX-WTX",33),rep("MX-WTX",33),rep("MX-ETX",33))

Ryx <- cbind(rev.pops,Ryx)

df.ryx <- rbind(data.frame(PopPairs=Ryx$rev.pops,Type="Deleterious",Rxy=Ryx$V1),
                data.frame(PopPairs=Ryx$rev.pops,Type="Tolerated",Rxy=Ryx$V2),
                data.frame(PopPairs=Ryx$rev.pops,Type="Functional",Rxy=Ryx$V3))

df.ryx$PopPairs <- factor(df.ryx$PopPairs, levels = c("WTX-AZ","ETX-AZ","MX-AZ","ETX-WTX","MX-WTX","MX-ETX"))

df.ryx.del <- df.ryx[which(df.ryx$Type=="Deleterious"),]
df.ryx.tol <- df.ryx[which(df.ryx$Type=="Tolerated"),]
df.ryx.fun <- df.ryx[which(df.ryx$Type=="Functional"),]

mu_del <- ddply(df.ryx.del, "PopPairs", summarise, Rxy.mean=mean(Rxy),Rxy.sd=sd(Rxy))
mu_tol <- ddply(df.ryx.tol, "PopPairs", summarise, Rxy.mean=mean(Rxy),Rxy.sd=sd(Rxy))
mu_fun <- ddply(df.ryx.fun, "PopPairs", summarise, Rxy.mean=mean(Rxy),Rxy.sd=sd(Rxy))

mean.rxy <- rbind(data.frame(PopPairs=mu_del$PopPairs,Type="Deleterious",MeanRxy=mu_del$Rxy.mean,SDRxy=mu_del$Rxy.sd),
                  data.frame(PopPairs=mu_tol$PopPairs,Type="Tolerated",MeanRxy=mu_tol$Rxy.mean,SDRxy=mu_tol$Rxy.sd),
                  data.frame(PopPairs=mu_fun$PopPairs,Type="Functional",MeanRxy=mu_fun$Rxy.mean,SDRxy=mu_fun$Rxy.sd))


df <- df.ryx[which(df.ryx$PopPairs =="WTX-AZ"),]
df$Type <- factor(df$Type , levels=c("Deleterious","Tolerated","Functional"))
p1 <- ggplot(df)+
  geom_boxplot(aes(x=Type,y=Rxy),position=position_dodge())+
  theme_classic(base_size = 15)+ ylim(0.95,1.15)+
  geom_hline(aes(yintercept=1),color="Black",linetype="dashed",size=1.1)+
  labs(x="Population Pairs (X-Y)", y = expression(paste(R[XY])))+
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=15),
        legend.text = element_blank(),
        legend.title=element_blank(),
        legend.position = "none")


p1

ggplot(Ryx)+
  geom_boxplot(aes(x=PopPairs,y=Rtol, fill=PopPairs))+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=carto_pal(12,"Vivid"))+
  geom_hline(aes(yintercept=1),color="Black",linetype="dashed",size=1.1)+
  labs(x="Population Pairs (X-Y)", y = expression(paste(R[XY])))+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_blank(),
        legend.title=element_blank(),
        legend.position = "none")


#### Draw SFS ####

az_del <- AZ_del[-which(AZ_del$MAF==0),]
wtx_del <- WTX_del[-which(WTX_del$MAF==0),]

#df.del <- inner_join(az_del,wtx_del,by=c("CHROM","POS"))

del_freq <- rbind(data.frame(Site="AZ",freq=az_del$MAF,All=round(az_del$MAF*az_del$N_CHR)),
                  data.frame(Site="WTX",freq=wtx_del$MAF,All=round(wtx_del$MAF*wtx_del$N_CHR)))

ggplot(del_freq,aes(x=freq, fill=Site)) +
  geom_histogram(alpha=0.7, position = 'identity', bins = 20, color="black") +
  scale_fill_manual(values=c(col_az,col_wtx)) +
  theme_classic(base_size = 20)+
  labs(x="Derived allele frequency", y = "No. of deleterious SNPs")+
  theme(axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.text = element_blank(),
        legend.title=element_blank(),
        legend.position = "none")


az_tol <- AZ_wdel[-which(AZ_wdel$MAF==0),]
wtx_tol <- WTX_wdel[-which(WTX_wdel$MAF==0),]

df.tol <- inner_join(az_tol,wtx_del,by=c("CHROM","POS"))

tol_freq <- rbind(data.frame(Site="AZ",freq=df.tol$Alt_freq.x),
                  data.frame(Site="WTX",freq=df.tol$Alt_freq.y))

ggplot(del_freq,aes(x=freq, fill=Site)) +
  geom_density(alpha=0.7, position = 'identity', color="black") +
  geom_histogram(alpha=0.7, position = 'identity', bins = 20, color="black") +
  scale_fill_manual(values=c(col_az,col_wtx)) +
  theme_classic(base_size = 20)+
  labs(x="Derived allele frequency", y = "No. of deleterious SNPs")+
  theme(axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.text = element_blank(),
        legend.title=element_blank(),
        legend.position = "none")

ggplot(tol_freq,aes(x=freq, fill=Site)) +
  geom_histogram(alpha=0.7, position = 'identity', bins = 20, color="black") +
  scale_fill_manual(values=c(col_az,col_wtx)) +
  theme_classic(base_size = 15)+
  labs(x="Derived allele frequency", y = "No. of tolerated SNPs")+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_blank(),
        legend.title=element_blank(),
        legend.position = "none")

#### Mean +/SD Freq ####

par(mar=c(4.5,7.1,2.5,5.5))
plot(NULL, xlim=c(0,2), ylim=c(0,1), 
     ylab=expression("Derived allele freqzuency"), xlab="Site", xaxt='n')
axis(side = 1, at = c(0.5,1.5), labels = c("AZ", "WTX"))

points(0.3,mean(az_del$Alt_freq),pch= 21,bg=col_az, col="black", cex=1.5)
arrows(0.3, mean(az_del$Alt_freq)-sd(az_del$Alt_freq), 0.3, mean(az_del$Alt_freq)+sd(az_del$Alt_freq), length=0.05, angle=90, code=3, lwd=1.5)

points(0.7,mean(az_tol$Alt_freq),pch= 21,bg=cols_l[1], col="black", cex=1.5)
arrows(0.7, mean(az_tol$Alt_freq)-sd(az_tol$Alt_freq), 0.7, mean(az_tol$Alt_freq)+sd(az_tol$Alt_freq), length=0.05, angle=90, code=3, lwd=1.5)

points(1.3,mean(wtx_del$Alt_freq),pch= 21,bg=col_az, col="black", cex=1.5)
arrows(1.3, mean(wtx_del$Alt_freq)-sd(wtx_del$Alt_freq), 1.3, mean(wtx_del$Alt_freq)+sd(wtx_del$Alt_freq), length=0.05, angle=90, code=3, lwd=1.5)

points(0.7,mean(az_tol$Alt_freq),pch= 21,bg=cols_l[1], col="black", cex=1.5)
arrows(0.7, mean(az_tol$Alt_freq)-sd(az_del$Alt_freq), 0.7, mean(az_del$Alt_freq)+sd(az_del$Alt_freq), length=0.05, angle=90, code=3, lwd=1.5)

