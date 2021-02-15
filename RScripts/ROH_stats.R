###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 09/09/20                  Last Modified: 09/09/20 ###
###########################################################################
###########################################################################
###                   ROH_plot.R         		                        ###
###########################################################################

#### PREREQUISITES #####
# load packages

#install.packages("detectRUNS", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")
#install.packages("DataCombine", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")

#install.packages("crayon", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")
#install.packages("pillar", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")

#library("detectRUNS", lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs" )
#library("DataCombine", lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs" )
#library(crayon, lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs" )
library("optparse")
# Load colors

#Option Parser

option_list = list(
  make_option(c("-f","--file"), type="character",default=NULL,help="Plink hom file",metavar="character"),
  make_option(c("-o","--out"), type="character",default=NULL,help="Output file",metavar="character"),
  make_option(c("-s","--homsnp"), type="numeric",default=NULL,help="only runs of homozygosity containing at least -s SNPs",metavar="numeric"),
  #make_option(c("-z","--gap"), type="numeric",default=NULL,help="if two consecutive SNPs are more than -g kb apart, they cannot be in the same ROH",metavar="numeric"),
  make_option(c("-d","--density"), type="numeric",default=NULL,help="a ROH must have at least one SNP per -d kb on average",metavar="numeric"),
  make_option(c("-k","--kb"), type="numeric",default=NULL,help="only runs of homozygosity of total length ≥ k kilobases are noted",metavar="numeric"),
  #make_option(c("-m","--miss"), type="numeric",default=NULL,help=" scanning window hit can contain at most -m missing calls",metavar="numeric"),
  make_option(c("-w","--win"), type="numeric",default=NULL,help=" the scanning window contains -w SNPs;",metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


#setwd("/scratch/snyder/m/mathur20/MQU/ch3/roh/plink1.9/stats/tryruns")

roh_ind <- read.csv(opt$file,header=F)
params <- cbind(rep(opt$homsnp, dim(roh_ind)[1]),
                #rep(opt$gap, dim(roh_ind)[1]),
                rep(opt$density, dim(roh_ind)[1]),
                rep(opt$kb, dim(roh_ind)[1]),
                #rep(opt$missing, dim(roh_ind)[1]),
                rep(opt$win, dim(roh_ind)[1]))
head(params)
ind_sum <-as.data.frame(cbind(roh_ind,params))

write.table(ind_sum,file=opt$out, quote=F, append=T)
#setwd("/scratch/snyder/m/mathur20/MQU/ch3/roh/plink1.9/final/")
#runs <- readExternalRuns(inputFile = "homsnp50.gap100.dens50.kb100.miss0.win20.hom", program = "plink")
#a <- noquote(unique(runs$group))
#r <- data.frame(from = c(a[1:28],a[29:30],a[31:34],a[35:65],a[66:68]),
#                to=c(rep("AZ",28),rep("NM",2),rep("MX",4),rep("WTX",31),rep("ETX",3)))
#runs_new <- FindReplace(data = runs, Var = "group", replaceData = r,
#                      from = "from", to = "to", exact = T)
#plot_Runs(runs = runs_new, suppressInds = F, separatePlots = F, savePlots = T)
