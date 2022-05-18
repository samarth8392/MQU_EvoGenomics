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

setwd("~/Documents/Thesis_Research/Final Results/Ch3/revise/load/freq/")

az_del_freq <- read.table("az.best66.deleterious.frq", header = T)
az_tol_freq <- read.table("az.best66.tolerated.frq", header = T)

wtx_del_freq <- read.table("wtx.best66.deleterious.frq", header = T)
wtx_tol_freq <- read.table("wtx.best66.tolerated.frq", header = T)


# Segregating load = Proportion of deleterious mutations present in a population that are segregating (i.e. MAF!=0 and MAF!=1)
az_del_tot <- length(which(az_del_freq$Alt_freq !=0)) # 15028
az_tol_tot <- length(which(az_tol_freq$Alt_freq !=0)) # 29272
wtx_del_tot <- length(which(wtx_del_freq$Alt_freq !=0)) # 14000
wtx_tol_tot <- length(which(wtx_tol_freq$Alt_freq !=0)) # 26452

segL_del_az <- length(az_del_freq$Alt_freq[which(az_del_freq$Alt_freq != 0 & az_del_freq$Alt_freq != 1)])/az_del_tot # 1
segL_del_wtx <- c(length(wtx_del_freq$Alt_freq[which(wtx_del_freq$Alt_freq !=0 & wtx_del_freq$Alt_freq != 1)]))/wtx_del_tot # 1

drift_del_az <- c(length(az_del_freq$Alt_freq[which(az_del_freq$Alt_freq == 1)]))/az_del_tot # 0
drift_del_wtx <- c(length(wtx_del_freq$Alt_freq[which(wtx_del_freq$Alt_freq == 1)]))/az_del_tot # 0

segL_tol_az <- length(az_tol_freq$Alt_freq[which(az_tol_freq$Alt_freq != 0 & az_tol_freq$Alt_freq != 1)])/az_tol_tot # 1
segL_tol_wtx <- length(wtx_tol_freq$Alt_freq[which(wtx_tol_freq$Alt_freq != 0 & wtx_tol_freq$Alt_freq != 1)])/wtx_tol_tot # 1

drift_tol_az <- c(length(az_tol_freq$Alt_freq[which(az_tol_freq$Alt_freq == 1)]))/az_tol_tot # 0
drift_tol_wtx <- c(length(wtx_tol_freq$Alt_freq[which(wtx_tol_freq$Alt_freq == 1)]))/wtx_tol_tot # 0

# To Draw SFS

daf <- seq(0.1,0.9,by=0.1)

Pazdel <- NULL
Pwtxdel <- NULL
Paztol <- NULL
Pwtxtol <- NULL

for (j in 1:length(daf))
{
  Pazdel[j] <- length(az_del_freq$Alt_freq[which(az_del_freq$Alt_freq== daf[j])])/az_del_tot
  Pwtxdel[j] <- length(wtx_del_freq$Alt_freq[which(wtx_del_freq$Alt_freq== daf[j])])/az_tol_tot
  Paztol[j] <- length(az_tol_freq$Alt_freq[which(az_tol_freq$Alt_freq== daf[j])])/wtx_del_tot
  Pwtxtol[j] <- length(wtx_tol_freq$Alt_freq[which(wtx_tol_freq$Alt_freq== daf[j])])/wtx_tol_tot
}

saf <- data.frame(DAF=daf,freq_AZ_del=Pazdel,freq_WTX_del=Pwtxdel,freq_AZ_tol=Paztol,freq_WTX_tol=Pwtxtol)
barplot(t(saf[,2:3]),beside=T, col=c(col_az,col_wtx))

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

