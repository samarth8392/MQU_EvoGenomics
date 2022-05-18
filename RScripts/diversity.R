###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 09/09/20                  Last Modified: 01/11/22 ###
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




#### PI: START #####


az <-read.table("az.best66.old.auto.noSing.nomiss_100kb.windowed.pi",header = T)
wtx <-read.table("wtx.best66.old.auto.noSing.nomiss_100kb.windowed.pi", header = T)


pi <-inner_join(az,wtx,by=c("CHROM","BIN_START","BIN_END"))
colnames(pi) <- c("CHROM","BIN_START","BIN_END","SNP_AZ","PI_AZ","SNP_WTX","PI_WTX")

#plot by chr
chroms <- unique(pi$CHROM)
for (i in chroms)
{
pi_chr <- pi[which(pi$CHROM == i),]
pdf(file = paste("../Tables&Figures/Pi/Pi_",i,".pdf",sep = ""), width = 12, height = 3)
plot(pi_chr$BIN_END/1e6,pi_chr$PI_AZ/pi_chr$SNP_AZ*1e5, col=adjustcolor(col_wtx,alpha.f = 1),type="l", lwd=1,
     xlab="Position (Mb)", ylab=expression(paste("Per site nucleotide diversity (",pi,x10^5,")")))
lines(pi_chr$BIN_END/1e6,pi_chr$PI_WTX/pi_chr$SNP_WTX*1e5, col=adjustcolor(col_az,alpha.f = 1),type="l", lwd=1)
dev.off()
}

az <- data.frame(Site="AZ", value = pi$PI_AZ)
wtx <- data.frame(Site="WTX", value = pi$PI_WTX)

pi_all <- rbind(az,wtx)
mu_all <- ddply(pi_all, "Site", summarise, grp.mean=mean(value))
# Compare means
my_comparisons <- list( c("AZ", "WTX"))
par(mar=c(5,10,2,10))
cp1 <- ggplot(pi_all, aes(x=value, fill=Site)) +
  geom_density(alpha=0.4) + xlim(0,0.003) + xlab(expression(paste("Per window", (pi))))+
  geom_vline(data=mu_all, aes(xintercept=grp.mean, color=Site),linetype="dashed") + theme_classic(base_size = 22)

cp1
#### PI: ENDS #####


#### FST: START #####
setwd("~/Documents/Thesis_Research/Final Results/Ch3/diversity/fst/")

p12 <- read.table("az.best68.txt_wtx.best68.txt.best68.auto.longchr.angsd.windowed.weir.fst", header = T)
p13 <- read.table("az.best68.txt_etx.best68.txt.best68.auto.longchr.angsd.windowed.weir.fst", header = T)
p14 <- read.table("az.best68.txt_nm.best68.txt.best68.auto.longchr.angsd.windowed.weir.fst", header = T)
p15 <- read.table("az.best68.txt_mx.best68.txt.best68.auto.longchr.angsd.windowed.weir.fst", header = T)
p23 <- read.table("wtx.best68.txt_etx.best68.txt.best68.auto.longchr.angsd.windowed.weir.fst", header = T)
p24 <- read.table("wtx.best68.txt_nm.best68.txt.best68.auto.longchr.angsd.windowed.weir.fst", header = T)
p25 <- read.table("wtx.best68.txt_mx.best68.txt.best68.auto.longchr.angsd.windowed.weir.fst", header = T)
p34 <- read.table("etx.best68.txt_nm.best68.txt.best68.auto.longchr.angsd.windowed.weir.fst", header = T)
p35 <- read.table("etx.best68.txt_mx.best68.txt.best68.auto.longchr.angsd.windowed.weir.fst", header = T)
p45 <- read.table("nm.best68.txt_mx.best68.txt.best68.auto.longchr.angsd.windowed.weir.fst", header = T)
pall <- read.table("best68.auto.longchr.angsd.windowed.weir.fst", header = T)


#Remove outliers
Q <- quantile(p12$WEIGHTED_FST, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(p12$WEIGHTED_FST)
p12_noout <- subset(p12, p12$WEIGHTED_FST > (Q[1] - 2*iqr) & p12$WEIGHTED_FST < (Q[2]+2*iqr))

Q <- quantile(p13$WEIGHTED_FST, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(p13$WEIGHTED_FST)
p13_noout <- subset(p13, p13$WEIGHTED_FST > (Q[1] - 2*iqr) & p13$WEIGHTED_FST < (Q[2]+2*iqr))

Q <- quantile(p14$WEIGHTED_FST, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(p14$WEIGHTED_FST)
p14_noout <- subset(p14, p14$WEIGHTED_FST > (Q[1] - 2*iqr) & p14$WEIGHTED_FST < (Q[2]+2*iqr))

Q <- quantile(p15$WEIGHTED_FST, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(p15$WEIGHTED_FST)
p15_noout <- subset(p15, p15$WEIGHTED_FST > (Q[1] - 2*iqr) & p15$WEIGHTED_FST < (Q[2]+2*iqr))

Q <- quantile(p23$WEIGHTED_FST, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(p23$WEIGHTED_FST)
p23_noout <- subset(p23, p23$WEIGHTED_FST > (Q[1] - 2*iqr) & p23$WEIGHTED_FST < (Q[2]+2*iqr))

Q <- quantile(p24$WEIGHTED_FST, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(p24$WEIGHTED_FST)
p24_noout <- subset(p24, p24$WEIGHTED_FST > (Q[1] - 2*iqr) & p24$WEIGHTED_FST < (Q[2]+2*iqr))

Q <- quantile(p25$WEIGHTED_FST, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(p25$WEIGHTED_FST)
p25_noout <- subset(p25, p25$WEIGHTED_FST > (Q[1] - 2*iqr) & p25$WEIGHTED_FST < (Q[2]+2*iqr))

Q <- quantile(p34$WEIGHTED_FST, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(p34$WEIGHTED_FST)
p34_noout <- subset(p34, p34$WEIGHTED_FST > (Q[1] - 2*iqr) & p34$WEIGHTED_FST < (Q[2]+2*iqr))

Q <- quantile(p35$WEIGHTED_FST, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(p25$WEIGHTED_FST)
p35_noout <- subset(p35, p35$WEIGHTED_FST > (Q[1] - 2*iqr) & p35$WEIGHTED_FST < (Q[2]+2*iqr))

Q <- quantile(p45$WEIGHTED_FST, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(p45$WEIGHTED_FST)
p45_noout <- subset(p45, p45$WEIGHTED_FST > (Q[1] - 2*iqr) & p45$WEIGHTED_FST < (Q[2]+2*iqr))

#Compare means
c_p12 <- data.frame(Site="AZ-WTX", value = p12_noout$WEIGHTED_FST)
c_p13 <- data.frame(Site="AZ-ETX", value = p13_noout$WEIGHTED_FST)
c_p14 <- data.frame(Site="AZ-NM", value = p14_noout$WEIGHTED_FST)
c_p15 <- data.frame(Site="AZ-MX", value = p15_noout$WEIGHTED_FST)
c_p23 <- data.frame(Site="WTX-ETX", value = p23_noout$WEIGHTED_FST)
c_p24 <- data.frame(Site="WTX-NM", value = p24_noout$WEIGHTED_FST)
c_p25 <- data.frame(Site="WTX-MX", value = p25_noout$WEIGHTED_FST)
c_p34 <- data.frame(Site="ETX-NM", value = p34_noout$WEIGHTED_FST)
c_p35 <- data.frame(Site="ETX-MX", value = p35_noout$WEIGHTED_FST)
c_p45 <- data.frame(Site="NM-MX", value = p45_noout$WEIGHTED_FST)

fst_all <- rbind(c_p12,c_p13,c_p14,c_p15,c_p23,c_p24,c_p25,c_p34,c_p35,c_p45)
mu_all <- ddply(fst_all, "Site", summarise, grp.mean=mean(value))
mu_all

res.aov <- aov(value ~ Site, data = fst_all)
plot(res.aov, 2)
summary(res.aov)
TukeyHSD(res.aov)

#Plot hist
par(mfrow=c(2,5))
hist(p12$WEIGHTED_FST, breaks = 100, xlab = "Weighted Fst", main = "AZ-WTX")
abline(v=mean(p12$WEIGHTED_FST), lty=2,col="red", lwd=2)
mtext(noquote(paste("Mean Fst=", round(mean(p12$WEIGHTED_FST),digits = 3))))

hist(p13$WEIGHTED_FST, breaks = 100, xlab = "Weighted Fst", main = "AZ-ETX")
abline(v=mean(p13$WEIGHTED_FST), lty=2,col="red", lwd=2)
mtext(noquote(paste("Mean Fst=", round(mean(p13$WEIGHTED_FST),digits = 3))))

hist(p14$WEIGHTED_FST, breaks = 100, xlab = "Weighted Fst", main = "AZ-NM")
abline(v=mean(p14$WEIGHTED_FST), lty=2,col="red", lwd=2)
mtext(noquote(paste("Mean Fst=", round(mean(p14$WEIGHTED_FST),digits = 3))))

hist(p15$WEIGHTED_FST, breaks = 100, xlab = "Weighted Fst", main = "AZ-MX")
abline(v=mean(p15$WEIGHTED_FST), lty=2,col="red", lwd=2)
mtext(noquote(paste("Mean Fst=", round(mean(p15$WEIGHTED_FST),digits = 3))))

hist(p23$WEIGHTED_FST, breaks = 100, xlab = "Weighted Fst", main = "WTX-ETX")
abline(v=mean(p23$WEIGHTED_FST), lty=2,col="red", lwd=2)
mtext(noquote(paste("Mean Fst=", round(mean(p23$WEIGHTED_FST),digits = 3))))

hist(p24$WEIGHTED_FST, breaks = 100, xlab = "Weighted Fst", main = "WTX-NM")
abline(v=mean(p24$WEIGHTED_FST), lty=2,col="red", lwd=2)
mtext(noquote(paste("Mean Fst=", round(mean(p24$WEIGHTED_FST),digits = 3))))

hist(p25$WEIGHTED_FST, breaks = 100, xlab = "Weighted Fst", main = "WTX-MX")
abline(v=mean(p25$WEIGHTED_FST), lty=2,col="red", lwd=2)
mtext(noquote(paste("Mean Fst=", round(mean(p25$WEIGHTED_FST),digits = 3))))

hist(p34$WEIGHTED_FST, breaks = 100, xlab = "Weighted Fst", main = "ETX-NM")
abline(v=mean(p34$WEIGHTED_FST), lty=2,col="red", lwd=2)
mtext(noquote(paste("Mean Fst=", round(mean(p34$WEIGHTED_FST),digits = 3))))

hist(p35$WEIGHTED_FST, breaks = 100, xlab = "Weighted Fst", main = "ETX-MX")
abline(v=mean(p35$WEIGHTED_FST), lty=2,col="red", lwd=2)
mtext(noquote(paste("Mean Fst=", round(mean(p35$WEIGHTED_FST),digits = 3))))

hist(p45$WEIGHTED_FST, breaks = 100, xlab = "Weighted Fst", main = "NM-MX")
abline(v=mean(p45$WEIGHTED_FST), lty=2,col="red", lwd=2)
mtext(noquote(paste("Mean Fst=", round(mean(p45$WEIGHTED_FST),digits = 3))))

#Get no outlier regions
p1213 <-inner_join(p12_noout,p13_noout,by=c("CHROM","BIN_START", "BIN_END"))
p121314 <-inner_join(p1213,p14_noout,by=c("CHROM","BIN_START", "BIN_END"))
p12131415 <-inner_join(p121314,p15_noout,by=c("CHROM","BIN_START", "BIN_END"))
p1213141523 <-inner_join(p12131415,p23_noout,by=c("CHROM","BIN_START", "BIN_END"))
p121314152324 <-inner_join(p1213141523,p24_noout,by=c("CHROM","BIN_START", "BIN_END"))
p12131415232425 <-inner_join(p121314152324,p25_noout,by=c("CHROM","BIN_START", "BIN_END"))
p1213141523242534 <-inner_join(p12131415232425,p34_noout,by=c("CHROM","BIN_START", "BIN_END"))
p121314152324253435 <-inner_join(p1213141523242534,p35_noout,by=c("CHROM","BIN_START", "BIN_END"))
p12131415232425343545 <-inner_join(p121314152324253435,p45_noout,by=c("CHROM","BIN_START", "BIN_END"))

noout <- p12131415232425343545[,c(1,2,3)]
write.csv(file = "No_Fst_outliers.txt", quote = F, noout)
