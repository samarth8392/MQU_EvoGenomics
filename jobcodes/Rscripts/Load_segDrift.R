##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/06/20                  Last Modified: 01/18/22 ###
###########################################################################
###########################################################################
###                     Load_segDrift.R                   		        ###
###########################################################################

# get segregating and drift load

library(optparse)

### Creating an argument parser
option_list = list(
  make_option(c("-d","--del"), type="character",default=NULL,help="path to pop1 deleterious count file (no header)",metavar="character"),
  make_option(c("-t","--tol"), type="character",default=NULL,help="path to pop1 tolerated count file (no header)",metavar="character"),
  make_option(c("-n","--chr"), type="numeric",default=NULL,help="sample size of pop1",metavar="numeric"),
  make_option(c("-o","--out"), type="character",default=NULL,help="Output name",metavar="character")
  )
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

del_count <- scan(opt$d)
tol_count <- scan(opt$t)

daf <- seq(0,2*opt$n, by=1)

# Make SFS for each mutation type
Pd <- NULL
Pt <- NULL
for (i in 1:length(daf))
{
	Pd[i] <- length(del_count[which(del_count== daf[i])])
	Pt[i] <- length(tol_count[which(tol_count== daf[i])])
}

sfs <- data.frame(MAC=daf,Nd=Pd,Nt=Pt)
write.table(sfs,filename=paste(opt$out,".SFS.txt",sep=""),quote=F,row.names=F)

# calculate segregating and drift load

segLd <- length(del_count[which(del_count != 0 & del_count != 2*opt$n)])
driftLd <- length(del_count[which(del_count == 2*opt$n)])
segLt <- length(tol_count[which(tol_count != 0 & tol_count != 2*opt$n)])
driftLt <- length(tol_count[which(tol_count == 2*opt$n)])

load <- data.frame(rbind(c(segLd,segLt),c(driftLd,driftLt)))
colnames(load) <- c("Deleterious","Tolerated")
rownames(load) <- c("Segregating Load","Drift Load")
write.table(load,filename=paste(opt$out,".segDrift.txt",sep=""),quote=F,row.names=T)