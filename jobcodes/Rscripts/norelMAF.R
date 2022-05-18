###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 09/09/20                  Last Modified: 02/03/22 ###
###########################################################################
###########################################################################
###                   hetperkb.R         		                     	###
###########################################################################

library(ggplot2)
library(plyr)
library(dplyr)
# Load colors

col_az <- "#CC3A8E" #carto_pal(12,"Vivid")[11]
col_wtx <- "#9f82ce" #carto_pal(7,"Purp")[5]
col_etx <-"#f2855d" #carto_pal(7,"Peach")[5]
col_mx <- "#3598FF"
cols <- c(col_az,col_wtx,col_etx, col_mx)
colss <-c(rep(col_az,28),rep(col_mx,4),rep(col_wtx,31),rep(col_etx,3))
cols_l <- adjustcolor(cols, alpha.f = 0.5)
col_azwtx <- "#90CCC2"

setwd("/scratch/bell/mathur20/MQU/ch3/revise/frequency/norel/daf")

meanMAF <- NULL
aztx_wgr <- NULL
aztx_syn <- NULL
aztx_nonsyn <- NULL
aztx_del <- NULL

for (i in 1:64)
{
	az_wgr <- scan(paste("comb",i,"/az.rand.comb",i,".MAF", sep=""))
	tx_wgr <- scan(paste("comb",i,"/wtx.unrel.comb",i,".MAF", sep=""))

	az_syn <- scan(paste("comb",i,"/az.rand.syn.comb",i,".MAF", sep=""))
	tx_syn <- scan(paste("comb",i,"/wtx.unrel.syn.comb",i,".MAF", sep=""))

	az_nonsyn <- scan(paste("comb",i,"/az.rand.nonsyn.comb",i,".MAF", sep=""))
	tx_nonsyn <- scan(paste("comb",i,"/wtx.unrel.nonsyn.comb",i,".MAF", sep=""))

	az_del <- scan(paste("comb",i,"/az.rand.del.comb",i,".MAF", sep=""))
	tx_del <- scan(paste("comb",i,"/wtx.unrel.del.comb",i,".MAF", sep=""))

	df <- rbind(data.frame(Comb=i,Site="AZ",Type="WGR",MAF=mean(az_wgr)),
		        data.frame(Comb=i,Site="TX",Type="WGR",MAF=mean(tx_wgr)),
		        data.frame(Comb=i,Site="AZ",Type="Synonymous",MAF=mean(az_syn)),
		        data.frame(Comb=i,Site="TX",Type="Synonymous",MAF=mean(tx_syn)),
		        data.frame(Comb=i,Site="AZ",Type="Nonsynonymous",MAF=mean(az_nonsyn)),
		        data.frame(Comb=i,Site="TX",Type="Nonsynonymous",MAF=mean(tx_nonsyn)),
		        data.frame(Comb=i,Site="AZ",Type="Deleterious",MAF=mean(az_del)),
		        data.frame(Comb=i,Site="TX",Type="Deleterious",MAF=mean(tx_del)))

	meanMAF <- rbind(meanMAF,df)

	diff <- az_wgr-tx_wgr
	aztx_wgr <- data.frame(cbind(aztx_wgr,diff))
	diff <- az_syn-tx_syn
	aztx_syn <- data.frame(cbind(aztx_syn,diff))
	diff <- az_nonsyn-tx_nonsyn
	aztx_nonsyn <- data.frame(cbind(aztx_nonsyn,diff))
	diff <- az_del-tx_del
	aztx_del <- data.frame(cbind(aztx_del,diff))

}

# All samples

az_all <- scan("../../az.MAF")
tx_all <- scan("../../wtx.MAF")

az_syn <- scan("/scratch/bell/mathur20/MQU/ch3/revise/load/freq/az.best66.synonymous.MAF")
tx_syn <- scan("/scratch/bell/mathur20/MQU/ch3/revise/load/freq/wtx.best66.synonymous.MAF")

az_nonsyn <- scan("/scratch/bell/mathur20/MQU/ch3/revise/load/freq/az.best66.nonsynonymous.MAF")
tx_nonsyn <- scan("/scratch/bell/mathur20/MQU/ch3/revise/load/freq/wtx.best66.nonsynonymous.MAF")

az_del <- scan("/scratch/bell/mathur20/MQU/ch3/revise/load/freq/az.best66.deleterious.MAF")
tx_del <- scan("/scratch/bell/mathur20/MQU/ch3/revise/load/freq/az.best66.deleterious.MAF")

true <- rbind(data.frame(Site="AZ",Type="WGR",MAF=mean(az_all)),
	data.frame(Site="TX",Type="WGR",MAF=mean(tx_all)),
	data.frame(Site="AZ",Type="Synonymous",MAF=mean(az_syn)*1.166667),
	data.frame(Site="TX",Type="Synonymous",MAF=mean(tx_syn)*1.291667),
	data.frame(Site="AZ",Type="Nonsynonymous",MAF=mean(az_nonsyn)*1.166667),
	data.frame(Site="TX",Type="Nonsynonymous",MAF=mean(tx_nonsyn)*1.291667),
	data.frame(Site="AZ",Type="Deleterious",MAF=mean(az_del)*1.166667),
	data.frame(Site="TX",Type="Deleterious",MAF=mean(tx_del)*1.291667))



meanMAF$Type <- factor(meanMAF$Type , levels=c("WGR","Synonymous","Nonsynonymous","Deleterious"))
meanMAF$Site <- factor(meanMAF$Site , levels=c("AZ","TX"))

true$Type <- factor(true$Type , levels=c("WGR","Synonymous","Nonsynonymous","Deleterious"))
true$Site <- factor(true$Site , levels=c("AZ","TX"))

pdf("../meanMAF.pdf",width=9, height=6)

ggplot(meanMAF)+
  geom_boxplot(aes(x=Type,y=MAF, fill=Site),position=position_dodge())+
  geom_point(aes(x=Type,y=MAF),data=true, shape=15, size=5, na.rm=TRUE, color="dodgerblue4",position = position_dodge(width = 0.9))+
  theme_classic(base_size = 15)+
  scale_fill_manual(values=cols[1:2])+
  labs(x="Mutation Type", y = expression("mean MAF"))+
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=15),
        legend.text = element_blank(),
        legend.title=element_blank(),
        legend.position = "none")

dev.off()





# WGR Hist

p0 <- ggplot(aztx_wgr)+
theme_classic(base_size = 15)+
  labs(y="No. of SNPs", x = expression(paste(delta," MAF = AZ-WTX")))+
  scale_x_continuous(breaks = seq(from = -0.5, to = 0.5, by = 0.05))
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=15),
        legend.text = element_blank(),
        legend.title=element_blank(),
        legend.position = "none")

p1 <- p0
p2 <- p0
p3 <- p0
p4 <- p0

for ( i in 1:64)
{
	a <- as.data.frame(aztx_wgr[,i])
	b <- as.data.frame(aztx_syn[,i])
	c <- as.data.frame(aztx_nonsyn[,i])
	d <- as.data.frame(aztx_del[,i])

	colnames(a) <- "Diff"
	colnames(b) <- "Diff"
	colnames(c) <- "Diff"
	colnames(d) <- "Diff"

	p1 <- p1 + geom_histogram(alpha=0.2, bins = 100, fill="gray50", aes(x=Diff),data=a)+
	geom_vline(xintercept = mean(a$Diff), linetype=2)

	p2 <- p2 + geom_histogram(alpha=0.2, bins = 100, fill="gray50", aes(x=Diff),data=b)+
	geom_vline(xintercept = mean(b$Diff), linetype=2)

	p3 <- p3 + geom_histogram(alpha=0.2, bins = 100, fill="gray50", aes(x=Diff),data=c)+
	geom_vline(xintercept = mean(c$Diff), linetype=2)
	p4 <- p4 + geom_histogram(alpha=0.2, bins = 100, fill="gray50", aes(x=Diff),data=d)+
	geom_vline(xintercept = mean(d$Diff), linetype=2)


}

pdf("../DiffNonsynonymous.pdf",width=9, height=6)

p3

dev.off()
