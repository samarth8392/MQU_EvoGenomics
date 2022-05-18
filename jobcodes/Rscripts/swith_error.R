##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/06/20                  Last Modified: 01/28/22 ###
###########################################################################
###########################################################################
###                     switch_error.R                   		        ###
###########################################################################

# get switch error between phasing runs

setwd("/scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/switch")

count1<-0
count2<-0
switch <- NULL
count3 <-1
chrmean <- NULL
for (i in c(1:10))
{
	count1 <- count1+1
	count2 <-0
	for (j in c(1:10))
	{
		count2 <- count2+1
		if (count2 > count1)
		{
			if (i != j)
			{
				for (chr in c(1:15,17:28))
				{
					file <- read.table(paste("chr",chr,".run",i,"-run",j,".diff.indv.switch",sep=""),header=T)
					a <- file$SWITCH
					switch <- cbind(switch,a)
				}
				switch_mean <- rowMeans(switch)
				chrmean <- rbind(chrmean,switch_mean)
			}
		}
	}
}

error <- data.frame(SampleID=file$INDV, error=colMeans(chrmean)*100)