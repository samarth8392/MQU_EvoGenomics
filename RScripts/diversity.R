###########################################################################
###                          Samarth Mathur, PhD                     	  ###
###                        The Ohio State University                 	  ###
###                                                                     ###
###########################################################################
###########################################################################
###                   diversity.R         		                          ###
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

setwd("~/Documents/Thesis_Research/Final Results/Ch3/revise/diversity/")

### Heterozygosity ###

het_all <- read.table("best66.old.auto.poly.nomiss.het", header = T)
total <- 960796788 # Total genomic length analyzed

het_all <- cbind(het_all,(het_all$N_SITES-het_all$O.HOM.)/total)
colnames(het_all)[6] <- "het"

az_het <- data.frame(Site="AZ", value = het_all$het[1:28])
mx_het <- data.frame(Site="MX", value = het_all$het[29:32])
wtx_het <- data.frame(Site="WTX", value = het_all$het[33:63])
etx_het <- data.frame(Site="CTX", value = het_all$het[64:66])

het <- rbind(az_het,mx_het,wtx_het,etx_het)

het$Site <- factor(het$Site , levels=c("AZ", "WTX", "CTX", "MX"))

kruskal.test(value ~ Site, data = het) # Kruskal-Wallis chi-squared = 20.149, df = 3, p-value = 0.0001581
pairwise.wilcox.test(het$value, het$Site,
                     p.adjust.method = "fdr")


mu_all <- ddply(het, "Site", summarise, grp.mean=mean(value))
sd_all <- ddply(het, "Site", summarise, grp.mean=sd(value))


my_comparisons <- list( c("AZ", "WTX"), c("AZ", "CTX"), c("AZ", "MX"),
                        c("WTX", "CTX"),c("WTX", "MX"), c("CTX", "MX"))
#my_comparisons <- list( c("AZ", "WTX"),  c("AZ", "CTX"),c("WTX", "CTX") )
par(mar=c(5,10,2,10))
p1<- ggplot(het, aes(x=Site, y=value, fill=Site))+
  geom_boxplot(lwd=1)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",method = "wilcox",
                     bracket.size=1, size=5)+ 
  theme_classic(base_size = 15)+ scale_fill_manual(values=cols[1:4])+
  labs(x="Site", y ="Heterozygosity")+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position = "none")
p1
  
#### Heterozygosit per kb ####

# Create bins
setwd("hetperkb/bins100kb/")
ideo <- getIdeo("galGal6")
chroms <- NULL
for (i in c(1:28,30:33))
{
  chr <- paste("chr",i,sep="")
  chrom_bins <- getBins(chr,ideo,stepSize=100*1000)
  write.table(chrom_bins,paste(chr,"_bin100kb.txt",sep=""),quote = F,row.names = F, sep="\t")
  a <- read.table(paste(chr,"_bin100kb.txt",sep=""),header = T)
  a <- a[,c(1,2,3)]
  colnames(a) <- c("CHROM","START","END")
  write.table(a,paste(chr,"_bin100kb.txt",sep=""),quote = F,row.names = F, sep="\t")
}


bychr <- read.table("../Meanhetperkb.txt",header = T)
chrom <- unique(bychr$CHROM)
df2 <- NULL
for (i in chrom[17:32])
{
  a <- bychr[which(bychr$CHROM==i),]
  df <- rbind(data.frame(Site="AZ",mean=a$Mean_AZ,sd=a$SD_AZ,chrom=a$CHROM),
              data.frame(Site="WTX",mean=a$Mean_WTX,sd=a$SD_WTX,chrom=a$CHROM),
              data.frame(Site="ETX",mean=a$Mean_ETX,sd=a$SD_ETX,chrom=a$CHROM),
              data.frame(Site="MX",mean=a$Mean_MX,sd=a$SD_MX,chrom=a$CHROM))
  
  df2 <- rbind(df2,df)
}

df2$chrom <- factor(df2$chrom,levels=rev(unique(df2$chrom)))
df2$Site <- factor(df2$Site,levels=c("AZ","WTX","ETX","MX"))
  
p <- ggplot(df2, aes(x=chrom, y=mean, fill=Site)) + 
    #geom_point(position=position_dodge(.9))+
  geom_bar(stat="identity", color="black",position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9))+
    theme_classic(base_size = 15)+coord_flip()+
    labs(x="Chromosome", y ="Mean het/kb")+scale_fill_manual(values=cols)+
  theme(legend.position = "none")
  p

