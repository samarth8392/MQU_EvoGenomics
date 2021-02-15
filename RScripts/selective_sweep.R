###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 09/09/20                  Last Modified: 09/09/20 ###
###########################################################################
###########################################################################
###                   selective_sweep.R         		                ###
###########################################################################

#### PREREQUISITES #####
#Install packages
#install.packages("ggsignif", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")
#install.packages("tidyr", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")
#install.packages("reshape2", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")
#install.packages("ggpubr", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")
#install.packages("cowplot", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")
#install.packages("plotrix", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")
#install.packages("wesanderson", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")
#install.packages("reshape2", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")
#install.packages("DataCombine", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")
#install.packages("DescTools", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")
#install.packages("rcartocolor", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")
#install.packages("car", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")


# load packages
library(ggplot2)
library(ggsignif, lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs")
library(reshape2, lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs")
#library(cowplot, lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs")
#library(plotrix, lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs")
library(wesanderson, lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs")
#library(plyr)
library(dplyr)
library(DataCombine, lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs")
#library(DescTools, lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs")
library(rcartocolor, lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs")
#library(car, lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs")

# Load colors

display_carto_all(colorblind_friendly = TRUE)
col_az <- carto_pal(7,"Teal")[7]
col_wtx <- carto_pal(7,"ag_Sunset")[3]
col_etx <-carto_pal(7,"Safe")[5]
col_nm <- carto_pal(7,"Peach")[7]
col_mx <- carto_pal(7,"BrwnYl")[7]
cols <- c(col_az,col_wtx,col_etx, col_nm, col_mx)
cols_l <- adjustcolor(cols, alpha.f = 0.5)
#### PREREQUISITES: END #####


##PLOT PI per KB
#Get_Pi
setwd("/scratch/snyder/m/mathur20/MQU/ch3/diversity/pi/per_kb/")

az_all <-read.table("az.best.angsd.final.1kb.windowed.pi", header = T)
wtx_all <-read.table("wtx.best.angsd.final.1kb.windowed.pi", header = T)
etx_all <-read.table("etx.best.angsd.final.1kb.windowed.pi", header = T)
nm_all <-read.table("nm.best.angsd.final.1kb.windowed.pi", header = T)
mx_all <-read.table("mx.best.angsd.best.1kb.windowed.pi", header = T)

az <-az_all[which(az_all$N_VARIANTS >=50),]
wtx <-wtx_all[which(az_all$N_VARIANTS >=50),]
etx <-etx_all[which(az_all$N_VARIANTS >=50),]
nm <-nm_all[which(az_all$N_VARIANTS >=50),]
mx <-mx_all[which(az_all$N_VARIANTS >=50),]

scaff_labes <- c(seq(1,28),32,33,30,31)

scaffs_az <- sort(as.vector(unique(az$CHROM)))
scaffs_wtx <- sort(as.vector(unique(wtx$CHROM)))
scaffs_etx <- sort(as.vector(unique(etx$CHROM)))
scaffs_nm <- sort(as.vector(unique(nm$CHROM)))
scaffs_mx <- sort(as.vector(unique(mx$CHROM)))

# Replace RefSeq names to chr names
#AZ
r <- data.frame(from = c(scaffs_az),
                to = c(scaff_labes))

az_new <- FindReplace(data = az, Var = "CHROM", replaceData = r,
					 from = "from", to = "to", exact = T)

wtx_new <- FindReplace(data = wtx, Var = "CHROM", replaceData = r,
                      from = "from", to = "to", exact = T)

etx_new <- FindReplace(data = etx, Var = "CHROM", replaceData = r,
                       from = "from", to = "to", exact = T)

nm_new <- FindReplace(data = nm, Var = "CHROM", replaceData = r,
                       from = "from", to = "to", exact = T)

mx_new <- FindReplace(data = mx, Var = "CHROM", replaceData = r,
                      from = "from", to = "to", exact = T)


# Order to 1-33

ord <- c(scaff_labes[1:28],scaff_labes[31],scaff_labes[32],scaff_labes[29],scaff_labes[30])

az_ord  = az_new[FALSE,]
for (i in 1:length(ord)){
  a <- az_new[which(az_new$CHROM == ord[i]),]
  az_ord <- rbind(az_ord,a)
}

wtx_ord  = wtx_new[FALSE,]
for (i in 1:length(ord)){
  a <- wtx_new[which(wtx_new$CHROM == ord[i]),]
  wtx_ord <- rbind(wtx_ord,a)
}

etx_ord  = etx_new[FALSE,]
for (i in 1:length(ord)){
  a <- etx_new[which(etx_new$CHROM == ord[i]),]
  etx_ord <- rbind(etx_ord,a)
}

nm_ord  = nm_new[FALSE,]
for (i in 1:length(ord)){
  a <- nm_new[which(nm_new$CHROM == ord[i]),]
  nm_ord <- rbind(nm_ord,a)
}

mx_ord  = mx_new[FALSE,]
for (i in 1:length(ord)){
  a <- mx_new[which(mx_new$CHROM == ord[i]),]
  mx_ord <- rbind(mx_ord,a)
}

# Pairwise ratio

#AZ-WTX

az_wtx <- inner_join(az_ord, wtx_ord,by=c("CHROM","BIN_START","BIN_END"), suffix = c("_az", "_wtx"))

pi_azwtx <- cbind.data.frame(az_wtx$CHROM,az_wtx$BIN_START,az_wtx$BIN_END,az_wtx$PI_az/az_wtx$PI_wtx)
colnames(pi_azwtx) <- c("CHROM","BIN_START","BIN_END","PI_azbywtx")
head(pi_azwtx)