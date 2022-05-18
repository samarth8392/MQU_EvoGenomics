###########################################################################
###                          Samarth Mathur, PhD                     	  ###
###                        The Ohio State University                 	  ###
###                                                                     ###
###########################################################################
###########################################################################
###                   heterozygosity_norel.R         		                ###
###########################################################################

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
library(tidyr)
library(DataCombine)
library(DescTools)
library(rcartocolor)
library(car)
library(detectRUNS)
library(gplots)
library(IdeoViz)
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

# All Samples heterozygosity

het_all <- read.table("~/Documents/Thesis_Research/Final Results/Ch3/revise/diversity/best66.auto.noSing.nomiss.het", header = T)
total <- 960796788 # Total genomic length analyzed

het_all <- cbind(het_all,(het_all$N_SITES-het_all$O.HOM.)/total)
colnames(het_all)[6] <- "het"

az_het <- data.frame(Site="AZ", value = het_all$het[1:28])
wtx_het <- data.frame(Site="WTX", value = het_all$het[33:63])

het <- rbind(az_het,wtx_het)
het$Site <- factor(het$Site , levels=c("AZ", "WTX"))

mu_all <- ddply(het, "Site", summarise, MeanHet=mean(value), sdHet=sd(value))


# # independent 2-group Mann-Whitney U Test
wilcox.test(az_het$value,wtx_het$value) # W = 676, p-value = 0.0001571

# Randomization (all samples)
var.test(az_het$value, wtx_rand) # F = 5.9168, num df = 27, denom df = 30, p-value = 6.715e-06
# This means Variance is non-homogeneous

### Testing significance using independent 2-group Mann-Whitney U Test ###
wtx_28 <- wtx_het[1:28,]
my.w <- wilcox.test(az_het$value,wtx_28$value)$statistic

# Randomization of length data to create null W-distribution #
my.null.w <-vector(length=9999)
my.w <-vector(length=9999)
s.size.paz <- length(az_het$value)
s.size.pwtx <- length(wtx_het$value)


for(w in 1:9999)
{
  wtx_rand <- sample(wtx_het$value, 28, replace=TRUE)
  my.w[w] <- wilcox.test(az_het$value,wtx_28$value)$statistic
  het_aztwx <- c(az_het$value,wtx_rand)
  my.rand.het <- sample(het_aztwx, length(het_aztwx), replace=FALSE)
  my.rand.data <- matrix(data=my.rand.het, nrow=28, ncol=2, byrow=FALSE)
  colnames(my.rand.data) <- c("AZ", "WTX")
  rownames (my.rand.data) <- seq(1:28)
  my.null.w[w] <- wilcox.test(my.rand.data[,1],my.rand.data[,2])$statistic
}

hist(my.null.w, col="blue", density = 8, 
     main = "Null W-Distribution",
     xlab = "W Statistic")
arrows(my.w,-280,my.w,0, xpd = TRUE, code= 2, length=0.15, col="red")
arrows(-my.t,-280,-my.t,0, xpd = TRUE, code= 2, length=0.15, col="red")
text(my.t,-350, col="red", labels = sprintf("%.2f", my.t), xpd = TRUE)
text(-my.t,-350, col="red", labels = sprintf("%.2f", -my.t), xpd = TRUE)
pval.t <- (length(which(my.null.t >= my.t))+
             + length(which(my.null.t <= -my.t)))/10000
pval.t  



# Related samples
#Sample 1	Sample 2
#E9031	E9036
#E9032	E9034
#E9038	E9039
#E9042	E9043
#E9046	E9048
#E9046	E9050
#E9048	E9051
#E9050	E9051
# All related = E9031, E9032, 
# E9046,E9048,E9050,E9051 = 1 Family group
# Selection choices = 2C1.2C1.2C1.2C1.4C1 = 2x2x2x2x4 = 64
wtxID <- het_all$INDV[33:63]
wtx_norel <- wtx_het[c(1,2,5,7,9,12,13,16,17,19,21,24:31),]
wtx_norelID <- wtxID[c(1,2,5,7,9,12,13,16,17,19,21,24:31)]

a <- c(3,8)
b <- c(4,6)
c <- c(10,11)
d <- c(14,15)
e <- c(18,20,22,23)
comb <- as.matrix(crossing(a, b, c, d, e))

wtxcomb <- NULL
for (i in 1:64)
{
  wtx_select <- c(wtxID[comb[i,]],wtx_norelID)
  wtxcomb <- as.data.frame(cbind(wtxcomb,wtx_select))
  colnames(wtxcomb)[i] <- paste("Comb",i,sep="")
}

write.table(wtxcomb,"../revise/Norel_WTX_combinations.txt",quote = F,row.names = F)

pval.t <- vector(length=64)
my.null.t <-vector(length=9999)
my.t <-vector(length=64)
h_wtx <-vector(length=64)
h_az <-vector(length=64)
iterations <- 10000
for (i in 1:64)
{
  wtx_select <- wtx_het[c(comb[i,]),]
  wtx_allnonrel <- rbind(wtx_norel,wtx_select)
  h_wtx[i] <- mean(wtx_allnonrel$value)
  resample_az <- sample(az_het$value, 24, replace=TRUE)
  h_az[i] <- mean(resample_az)
  my.t[i] <- t.test(resample_az,wtx_allnonrel$value, var.equal = F, alternative = "greater")$statistic
  #my.t[i] <- h_az[i]-h_wtx[i]
  #my.t[i] <- wilcox.test(resample_az,wtx_allnonrel$value, alternative = "greater", correct = T)$statistic
  #my.t[i] <- kruskal.test(resample_az,wtx_allnonrel$value, alternative = "greater", correct = T)$statistic
  for (t in 1 : iterations) {
    het_aztwx <- c(resample_az,wtx_allnonrel$value)
    my.rand.het <- sample(het_aztwx, length(het_aztwx), replace=FALSE)
    my.rand.data <- matrix(data=my.rand.het, nrow=24, ncol=2, byrow=FALSE)
    colnames(my.rand.data) <- c("AZ", "WTX")
    rownames (my.rand.data) <- seq(1:24)
    my.null.t[t] <- t.test(my.rand.data[,1],my.rand.data[,2], var.equal = F, alternative = "greater")$statistic
    #my.null.t[t] <- mean(my.rand.data[,1])-mean(my.rand.data[,2])
    #my.null.t[t] <- wilcox.test(my.rand.data[,1],my.rand.data[,2], alternative = "greater", correct = T)$statistic
    #my.null.t[t] <- kruskal.test(my.rand.data[,1],my.rand.data[,2], alternative = "greater", correct = T)$statistic
  }
  pval.t[i] <- (length(which(my.null.t >= my.t[i]))+ length(which(my.null.t <= -my.t[i])))/10000
}

hist(pval.t, density = 40, col="orange", border = 2, ylim=c(0,60), breaks = 20,
     main = "Difference in mean genome-wide heterozygosity (AZ-WTX)", 
     xlab = "p-value", xaxt='n')
axis(1,at=seq(0,0.015,by=0.001), labels =seq(0,0.015,by=0.001) )
box()

az_random <- data.frame(Site="AZ", value = h_az)
wtx_unrelated <- data.frame(Site="WTX", value = h_wtx)
het_sub <- rbind(az_random,wtx_unrelated)
my_comparisons <- list( c("AZ", "WTX"))
p3<- ggplot(het_sub, aes(x=Site, y=value, fill=Site))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif",method = "wilcox",
                     bracket.size=1.2, size=10)+ ylim(0.004,0.006)+
  theme_classic(base_size = 15)+ scale_fill_manual(values=c(col_az,col_wtx))+
  labs(x="Site", y ="Heterozygosity")+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position = "none")
p3
