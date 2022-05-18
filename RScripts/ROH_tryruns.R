###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 09/09/20                  Last Modified: 01/11/22 ###
###########################################################################
###########################################################################
###                   ROH_tryruns.R  		                                ###
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

setwd("~/Documents/Thesis_Research/Final Results/Ch3/revise/rohs/")


##### Sensitivity functions #####
# From: https://complementarytraining.net/simple-sensitivity-analysis-with-r/

source("../../../../Codes and Pipelines/simple-sensitivity-analysis-with-r/sensitivity-analysis.R")


##### TRY RUNS #####
library(sensitivity)

setwd("tryruns/")
### ROUND 1 ###
# Get parameters
columns <- c("snp","den","gap","het","wsnp","whet","wmiss","wthresh")
params <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(params) = columns

for (snp in c(50,100,500))
{
  for (den in c(10,50,100))
  {
    for (gap in c(500,1000,5000))
    {
      for (het in c(0,2,5))
      {
        for (wsnp in c(10,50,100))
        {
          for (whet in c(0,1,5))
          {
            for (wmiss in c(0,5,10))
            {
              for (wthresh in c(0.01,0.05,0.1))
              {
                p <- c(snp,den,gap,het,wsnp,whet,wmiss,wthresh)
                params <- rbind(params,p)
              }
            }
          }
        }
      }
    }
  }
}
colnames(params) <- columns
# Get files and calculate mean Froh for each population for each set of parameters

filenames <- list.files("tryruns", pattern="*indiv", full.names=TRUE)

total <- 960796788 # Total genomic length analyzed
rohs <- NULL

for (i in 1:length(filenames))
{
  a <- read.table(filenames[i], header = T)
  roh <- a$KB/total
  m_az <- mean(roh[1:28])
  m_mx <- mean(roh[29:32])
  m_wtx <- mean(roh[33:63])
  m_etx <- mean(roh[64:66])
  means <- cbind(m_az,m_mx,m_wtx,m_etx)
  rohs <- rbind(rohs,means)
}
popmeans <- cbind(rohs,params[1:nrow(rohs),])
write.table(popmeans,"popmeans.tryruns.txt",quote=F,row.names=F)

popmeans <- read.table("popmeans.tryruns.txt", header = T)
popmeans <- read.table("popmeans.tryruns2.txt", header = T)
popmeans <- read.table("popmeans.tryruns3.txt", header = T)
popmeans <- read.table("popmeans.tryruns4.txt", header = T)

#lm_az <- lm(m_az ~ snp + den + gap + het + wsnp + whet + wmiss + wthresh, data=popmeans) 
#lm_mx <- lm(m_mx ~ snp + den + gap + het + wsnp + whet + wmiss + wthresh, data=popmeans) 
#lm_wtx <- lm(m_wtx ~ snp + den + gap + het + wsnp + whet + wmiss + wthresh, data=popmeans) 
#lm_etx <- lm(m_etx ~ snp + den + gap + het + wsnp + whet + wmiss + wthresh, data=popmeans) 

# Standardized Rank Regression Coefficients
targets <- popmeans[,1:4]
predictors <- popmeans[,5:12]
#predictors <- popmeans[,5:9]
predictors <- popmeans[,5:7]
predictors <- as.data.frame(popmeans[,5])
colnames(predictors) <- "snp"

src_az <- src(predictors,targets[,1],rank = TRUE, logistic = TRUE, nboot = 1000, conf = 0.95)
src_mx <- src(predictors,targets[,2],rank = TRUE, logistic = TRUE, nboot = 1000, conf = 0.95)
src_wtx <- src(predictors,targets[,3],rank = TRUE, logistic = TRUE, nboot = 1000, conf = 0.95)
src_etx <- src(predictors,targets[,4],rank = TRUE, logistic = TRUE, nboot = 1000, conf = 0.95)

# Plot SRC indices
par(mar=c(4.5,7.1,2.5,5.5))
plot(NULL, xlim=c(0,3), ylim=c(-0.85,0.45), 
     ylab=expression("Standardized Regression Coefficients"), xlab="Parameter", xaxt='n')
#axis(side = 1, at = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5), labels = c("SNP", "den", "gap","het","wsnp","whet","wmiss","wthresh"))
#axis(side = 1, at = c(0.5,1.5,2.5,3.5,4.5), labels = c("SNP", "het", "wsnp","whet","wmiss"))
axis(side = 1, at = c(0.5,1.5,2.5), labels = c("SNP", "het","wmiss"))

for (i in c(1:3))
{
  points(i-0.6,src_az$SRC[,1][i],pch= 21,bg=col_az, col="black", cex=1.5)
  arrows(i-0.6, src_az$SRC[,4][i], i-0.6, src_az$SRC[,5][i], length=0.05, angle=90, code=3, lwd=1.5)
  
  points(i-0.5,src_wtx$SRC[,1][i],pch= 21,bg=col_wtx, col="black", cex=1.5)
  arrows(i-0.5, src_wtx$SRC[,4][i], i-0.5, src_wtx$SRC[,5][i], length=0.05, angle=90, code=3, lwd=1.5)
  
  points(i-0.4,src_etx$SRC[,1][i],pch= 21,bg=col_etx, col="black", cex=1.5)
  arrows(i-0.4, src_etx$SRC[,4][i], i-0.4, src_etx$SRC[,5][i], length=0.05, angle=90, code=3, lwd=1.5)
  
  points(i-0.3,src_mx$SRC[,1][i],pch= 21,bg=col_mx, col="black", cex=1.5)
  arrows(i-0.3, src_mx$SRC[,4][i], i-0.3, src_mx$SRC[,5][i], length=0.05, angle=90, code=3, lwd=1.5)
}
abline(h=0, lty=2, col="black")


# Plot Sensitivity plots (for mean of means)
target <- rowMeans(targets)
ctrl <- trainControl(method="cv", allowParallel = TRUE, number = 100, verboseIter = TRUE)
model <- train(y = targets[,1], x = predictors, 
               method = "gam",
               preProcess = c("center", "scale"),
               trControl = ctrl)
results <- sensitivityAnalysisCaret(model, targets[,1], predictors, targetPrediction = "ratio", from = 0.5,to =1.5)
plot(results$ggplot + ylab(expression(paste("Normalized change in target (",F[ROH],")"))))
