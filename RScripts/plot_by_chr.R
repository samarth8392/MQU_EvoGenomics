#install.packages("IdeoViz", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")

library(IdeoViz, lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs" )
library(dplyr)
cols <- c("#f8b58b","#FF0000","#3a7c89","#A16928")

setwd("/scratch/snyder/m/mathur20/MQU/ch3/angsd/hwe")

ideo <- read.table("Chicken_ideo")

az <- read.table("az.best66.Hobs.rename", header=T)
wtx <- read.table("wtx.best66.Hobs.rename", header=T)
etx <- read.table("etx.best66.Hobs.rename", header=T)
#mx <- read.table("mx.best66.Hobs.rename", header=T)

het1 <-inner_join(az,wtx,by=c("CHROM","POS"))
het2 <-inner_join(het1,etx,by=c("CHROM","POS"))

chroms <- unique(het2$CHROM)
for (i in chroms)
{
chrom_bins <- getBins(i,ideo,stepSize=1*1*1000)
avg_het <- avgByBin(data.frame(value=het2[,c(3:4)]),
                    pi[,1:2], chrom_bins)
pdf(file = paste(i,".pdf",sep = ""), width = 12, height = 3)
plotOnIdeo(chrom=seqlevels(chrom_bins),
           ideoTable=ideo,
           values_GR=avg_het, value_cols=colnames(mcols(avg_het))[2:4],
           val_range=c(0,0.01), ylab="Het/kb"),
           col=c(cols))
dev.off()
}