###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 09/09/20                  Last Modified: 09/09/20 ###
###########################################################################
###########################################################################
###                   admixture.R         		                          ###
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

setwd("~/Documents/Thesis_Research/Final Results/Ch3/revise/admix/")
#make plots
#choose colors
pal <- c(cols[2],cols[1],cols[3],cols[4],carto_pal(7,"PurpOr")[1],carto_pal(7,"Peach")[1])
pal <- c(cols[2],cols[1],cols[3],cols[4],carto_pal(7,"PurpOr")[1],carto_pal(7,"Peach")[1])

pal <-c(pal[2],pal[1],"#4472C4")


#ADMIXTURE
# Current labeling: 28AZ, 4MX, 32WTX, 3ETX
for (k in c(2:10))
{
  name <- paste("best66.old.auto.noSing.nomiss.rename.",k, sep = "")
  f_name <- paste("best66.old.auto.noSing.nomiss.rename.",k,".Q",sep = "")
  f <- as.data.frame(read.table(f_name, header = F))
  assign(name,f)
}
par(mar=c(1,10,1,10))
par(mfrow=c(5,1))
for (k in c(2:6))
{
  name <- get(noquote(paste("best66.old.auto.noSing.nomiss.rename.",k, sep = "")))
  #f_new <- rbind(name[33:63,],name[64:66,],name[1:28,],name[29:32,])
  barplot(t(name),width = 1, col=pal[1:k],cex.names=0.75, space = 0, border = "gray25",
          ylab = "Ancestry", cex.axis = 1.2,cex.lab=1.4,xaxt='n')
  abline(v=28,lwd=3)
  abline(v=32,lwd=3)
  abline(v=63,lwd=3)
}

