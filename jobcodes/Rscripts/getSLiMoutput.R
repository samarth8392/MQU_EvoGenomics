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


# Calculate Rxy
# We began by calculating a statistic which
#compares two populations, given a particular category of sites, in terms of the number of
#derived alleles found at sites within that category in one population rather than the other.
#Specifically, at each site i we write the observed derived allele frequency in population A and Population B,
#Then if C is a particular category of protein-coding sites and I a set of intergenic sites, we define
# L(A,B) = sum(fiA(1-fiB))/sum(fjA(1-fjB))
#where i = deleterious site and j = neutral sites (synonymous sites)

setwd("/scratch/bell/mathur20/MQU/ch3/revise/slim/1k/final_output/geno")

gens <- c(125002,125252,125502,125752,126002,126302,126502)
inds <- NULL
for (i in 0:199)
{
	inds[i+1] <- paste("i",i,sep="")
}




df.rxy <- NULL

for (gen in gens[6:7])
{
	for (i in 1:100)
	{
		vcf <- read.table(paste("gen",gen,".run",i,".txt",sep=""),header=F)
		s <- scan(paste("gen",gen,".run",i,".selcoef.txt",sep=""))
		p1 <- scan(paste("gen",gen,".run",i,".p1.freq",sep=""))
		p2 <- scan(paste("gen",gen,".run",i,".p2.freq",sep=""))
		df <- rbind(data.frame(Chr=vcf$V1,Pos=vcf$V2,S=s,Pop="P1",MAF=p1),
			data.frame(Chr=vcf$V1,Pos=vcf$V2,S=s,Pop="P2",MAF=p2))
		
		df2 <- data.frame(Chr=vcf$V1,Pos=vcf$V2,S=s,MAF1=p1,MAF2=p2)

		neut <- df2[which(df2$S == 0),]
		del <- df2[which(df2$S < 0),]
		ben <- df2[which(df2$S > 0),]

		neut.sites <- neut[,c(1,2)]
		del.sites <- del[,c(1,2)]
		ben.sites <- ben[,c(1,2)]

		write.table(neut.sites,paste("../sites/gen",gen,".run",i,".neutral.sites",sep=""),quote=F,row.names=F)
		write.table(del.sites,paste("../sites/gen",gen,".run",i,".del.sites",sep=""),quote=F,row.names=F)
		write.table(ben.sites,paste("../sites/gen",gen,".run",i,".ben.sites",sep=""),quote=F,row.names=F)

		l_del_12 <- sum(del$MAF1*(1-del$MAF2))/sum(neut$MAF1*(1-neut$MAF2))
		l_del_21 <- sum(del$MAF2*(1-del$MAF1))/sum(neut$MAF2*(1-neut$MAF1))
		r_del_12 <- l_del_12/l_del_21

		l_ben_12 <- sum(ben$MAF1*(1-ben$MAF2))/sum(neut$MAF1*(1-neut$MAF2))
		l_ben_21 <- sum(ben$MAF2*(1-ben$MAF1))/sum(neut$MAF2*(1-neut$MAF1))
		r_ben_12 <- l_ben_12/l_ben_21

		df4<- data.frame(Gen=gen,Run=i,Rxy_Del=r_del_12,Rxy_Ben=r_ben_12)

		df.rxy <- rbind(df.rxy,df4)
	}
}

df.ryx <- data.frame(Gen=df.rxy$Gen,Run=df.rxy$Run,Ryx_Del=1/df.rxy$Rxy_Del,Ryx_Ben=1/df.rxy$Rxy_Ben)

write.table(df.rxy,"../Rxy.sims.txt",quote=F,row.names=F)
write.table(df.ryx,"../Ryx.sims.txt",quote=F,row.names=F)



# Age distribution of last VCF

df.age <- NULL
for (i in 1:100)
{
	vcf <- read.table(paste("gen126502.run",i,".txt",sep=""),header=F)
	s <- scan(paste("gen126502.run",i,".selcoef.txt",sep=""))
	p1 <- scan(paste("gen126502.run",i,".p1.freq",sep=""))
	p2 <- scan(paste("gen126502.run",i,".p2.freq",sep=""))
	age <- scan(paste("gen126502.run",i,".age",sep=""))

	#df <- rbind(data.frame(Run=i, Chr=vcf$V1,Pos=vcf$V2,S=s,Age=age,Pop="P1",MAF=p1),
	#		    data.frame(Run=i, Chr=vcf$V1,Pos=vcf$V2,S=s,Age=age,Pop="P2",MAF=p2))

	df2 <- data.frame(Run=i, Chr=vcf$V1, Pos=vcf$V2,S=s,Age=age, MAF1=p1, MAF2=p2)

	del <- df2[which(df2$S < 0),]
	ben <- df2[which(df2$S > 0),]

	shared.del <- del[which(del$MAF1 != 0 & del$MAF2 != 0),]
	priv.az.del <- del[which(del$MAF1 != 0 & del$MAF2 == 0),]
	priv.tx.del <- del[which(del$MAF1 == 0 & del$MAF2 != 0),]

	shared.ben <- ben[which(ben$MAF1 != 0 & ben$MAF2 != 0),]
	priv.az.ben <- ben[which(ben$MAF1 != 0 & ben$MAF2 == 0),]
	priv.tx.ben <- ben[which(ben$MAF1 == 0 & ben$MAF2 != 0),]

	age_df <- rbind(data.frame(Run=i, Type="Deleterious", Site="Shared", Age = shared.del$Age, MAF1=shared.del$MAF1,MAF2=shared.del$MAF2),
                data.frame(Run=i, Type="Deleterious", Site="Private_AZ", Age = priv.az.del$Age,MAF1=priv.az.del$MAF1,MAF2=priv.az.del$MAF2),
                data.frame(Run=i, Type="Deleterious", Site="Private_TX", Age = priv.tx.del$Age,MAF1=priv.tx.del$MAF1,MAF2=priv.tx.del$MAF2),

                data.frame(Run=i, Type="Beneficial", Site="Shared", Age = shared.ben$Age, MAF1=shared.ben$MAF1,MAF2=shared.ben$MAF2),
                data.frame(Run=i, Type="Beneficial", Site="Private_AZ", Age = priv.az.ben$Age,MAF1=priv.az.ben$MAF1,MAF2=priv.az.ben$MAF2),
                data.frame(Run=i, Type="Beneficial", Site="Private_TX", Age = priv.tx.ben$Age,MAF1=priv.tx.ben$MAF1,MAF2=priv.tx.ben$MAF2))

	df.age <- rbind(df.age,age_df)

}


write.table(df.age,"../Age.sims.txt",quote=F,row.names=F)





####### EXTRA ##########

# Calculate admixture
# Admixture was run for each run and for each generation at K=2
# To get mean admixture proportions from each run, we averaged the admixture proportions

setwd("/scratch/bell/mathur20/MQU/ch3/revise/slim/1k/final_output/admix")

df.admix <- NULL
for (gen in gens)
{
	for (i in 1:100)
	{
		
		
		file <- read.table(paste("1k.gen",gen,".run",i,".2.Q",sep=""))
		df <- data.frame(Gen=gen,Run=i,Ind=inds,P1=file$V1,P2=file$V2)

		df.admix <- rbind(df.admix,df)
	}
}


for (gen in gens)
{
	for (run in 1:100)
	{
		df <- df.admix[which(df.admix$Gen == gen),]
		admix.mean <- ddply(df, "Run", summarise, meanP1=mean(P1),meanP2=mean(P2))
		pdf(file=paste("../Admix.gen",gen,".pdf",sep=""),width=6,height=4)
		barplot(t(admix.mean),width = 1, col=cols[1:2],cex.names=0.75, space = 0, border = "gray25",
          ylab = "Ancestry", cex.axis = 1.2,cex.lab=1.4,xaxt='n')
		dev.off()
	}
}










