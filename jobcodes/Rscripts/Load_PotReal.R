##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/06/20                  Last Modified: 01/18/22 ###
###########################################################################
###########################################################################
###                     Load_PotRead.R                   		        ###
###########################################################################

# Load colors

col_az <- "#CC3A8E" #carto_pal(12,"Vivid")[11]
col_wtx <- "#9f82ce" #carto_pal(7,"Purp")[5]
col_etx <-"#f2855d" #carto_pal(7,"Peach")[5]
col_mx <- "#3598FF"
cols <- c(col_az,col_wtx,col_etx, col_mx)
colss <-c(rep(col_az,28),rep(col_mx,4),rep(col_wtx,31),rep(col_etx,3))
cols_l <- adjustcolor(cols, alpha.f = 0.5)
col_azwtx <- "#90CCC2"


# Convert genotypes into allele numbers

#devtools::install_github('christophergandrud/DataCombine')
#devtools::install_github("kassambara/ggpubr")

library(DataCombine)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("/scratch/bell/mathur20/MQU/ch3/revise/load/geno")

r <- data.frame(from = c("0/0","0/1","1/1","./."), 
                to = c(0,1,2,NA))


# Get files
for (pop in c("az","wtx","etx","mx"))
{
	for (mut in c("deleterious","tolerated","synonymous","nonsynonymous"))
	{
		file <- read.table(paste(pop,".best66.",mut,".GT.FORMAT",sep=""),header=T)
		for (j in colnames(file)[3:dim(file)[2]])
		{
			file <- FindReplace(data = file, Var = j, replaceData = r,
                                from = "from", to = "to", exact = T)
		}
		f <- paste(pop,"_",mut,sep="")
		assign(f,file)
	}
}

#No. of sites

for (pop in c("az","wtx","etx","mx"))
{
	for (mut in c("deleterious","tolerated","synonymous","nonsynonymous"))
	{
		freq <- scan(paste("../freq/",pop,".best66.",mut,".MAF",sep=""))
		file <- get(paste(pop,"_",mut,sep=""))
		file <- cbind(file,freq)

		file <- file[-which(file$freq ==0),]
		name <- paste(pop,"_",mut,"_sites",sep="")
		sites <- vector()
		count <- 1
		for (j in 3:(dim(file)[2]-1))
		{
			a <- file[,j]
			a <- a[-which(is.na(a))]
			sites[count] <- length(a)
			count=count+1
		}
		assign(name,sites)
	}
}

#  No. of hetero

for (pop in c("az","wtx","etx","mx"))
{
	for (mut in c("deleterious","tolerated","synonymous","nonsynonymous"))
	{
		file <- get(paste(pop,"_",mut,sep=""))
		name <- paste(pop,"_",mut,"_hetero",sep="")
		sites <- vector()
		count <- 1
		for (j in 3:dim(file)[2])
		{
			a <- file[,j]
			a <- a[-which(is.na(a))]
			a <- a[which(a == "1")]
			sites[count] <- length(a)
			count=count+1
		}
		assign(name,sites)
	}
}

#  No. of alt homo

for (pop in c("az","wtx","etx","mx"))
{
	for (mut in c("deleterious","tolerated","synonymous","nonsynonymous"))
	{
		file <- get(paste(pop,"_",mut,sep=""))
		name <- paste(pop,"_",mut,"_homo",sep="")
		sites <- vector()
		count <- 1
		for (j in 3:dim(file)[2])
		{
			a <- file[,j]
			a <- a[-which(is.na(a))]
			a <- a[which(a == "2")]
			sites[count] <- length(a)
			count=count+1
		}
		assign(name,sites)
	}
}

# No.of alleles

for (pop in c("az","wtx","etx","mx"))
{
	for (mut in c("deleterious","tolerated","synonymous","nonsynonymous"))
	{
		file <- get(paste(pop,"_",mut,sep=""))
		name <- paste(pop,"_",mut,"_all",sep="")
		sites <- vector()
		count <- 1
		for (j in 3:dim(file)[2])
		{
			a <- file[,j]
			a <- a[-which(is.na(a))]
			b <- a[which(a == "1")]
			c <- a[which(a == "2")]
			sites[count] <- length(b) + 2*length(c)
			count=count+1
		}
		assign(name,sites)
	}
}


# Combine dataset of 
# type of mutation (type)
# No. of sites (snps)
# No. of hetero (het)
# No. of alt homo (althom)
# No.of alleles (all)

df2 <- data.frame(matrix(ncol = 6, nrow = 0))
df3 <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df2) <- c("Site","type","snps","het","althom","all")
colnames(df3) <- c("Site","type","snps","het","althom","all")

for (pop in c("az","wtx","etx","mx"))
{
	for (mut in c("deleterious","tolerated","synonymous","nonsynonymous"))
	{
		df <- rbind(data.frame(Site=pop,type=mut,snps=get(paste(pop,"_",mut,"_sites",sep="")),
			het=get(paste(pop,"_",mut,"_hetero",sep="")),althom=get(paste(pop,"_",mut,"_homo",sep="")),
			all=get(paste(pop,"_",mut,"_all",sep=""))))

		df2 <- rbind(df2,df)
	}
		df3 <- rbind(df3,df2)

}

df3 <- df3[c(1:112,225:348,585:596,845:860),]

# total number of sites
az_total <- az_deleterious_sites+az_tolerated_sites+az_nonsynonymous_sites
wtx_total <- wtx_deleterious_sites+wtx_tolerated_sites+wtx_nonsynonymous_sites
etx_total <- etx_deleterious_sites+etx_tolerated_sites+etx_nonsynonymous_sites
mx_total <- mx_deleterious_sites+mx_tolerated_sites+mx_nonsynonymous_sites

df4 <- cbind(df3,c(az_total,wtx_total,etx_total,mx_total))
colnames(df4)[7] <- "TotalSnps"

write.table(df4,"../LoadDatafile.txt",quote=F,row.names=F)


# Get frequency for Rxy statistic

setwd("/scratch/bell/mathur20/MQU/ch3/revise/load/freq")
for (pop in c("az","wtx","etx","mx"))
{
	for (mut in c("deleterious","tolerated","synonymous","nonsynonymous"))
	{
		file <- read.table(paste(pop,".best66.",mut,".frq",sep=""), header=F)
		freq <- scan(paste(pop,".best66.",mut,".MAF",sep=""))
		file <- cbind(file,freq)
		file <-  file[,c(1,2,4,7)]
		colnames(file) <- c("CHROM","POS","N_CHR","MAF")
		write.table(file,paste(pop,".best66.",mut,".freq",sep=""),quote=F,row.names=F)
	}
}


############# SHARED PRIVATE ###################

# get genotypes for shared and private mutations

setwd("/scratch/bell/mathur20/MQU/ch3/revise/age/geno")

r <- data.frame(from = c("0/0","0/1","1/1","./."), 
                to = c(0,1,2,NA))


# Get files
for (f in c("az.del.","wtx.del.","priv.az.del.","priv.wtx.del.","shared.del."))
{
	file <- read.table(paste(f,"GT.FORMAT",sep=""),header=T)	
	name <- paste(f,"geno",sep="")
	for (j in colnames(file)[3:dim(file)[2]])
	{
		file <- FindReplace(data = file, Var = j, replaceData = r,
                            from = "from", to = "to", exact = T)
	}
	assign(name,file)
}


# No. of heterozygotes
for (f in c("az.del.","wtx.del.","priv.az.del.","priv.wtx.del.","shared.del."))
{
	file <- get(paste(f,"geno",sep=""))
	name <- paste(f,"hetero",sep="")
	sites <- vector()
	count <- 1
	for (j in 3:dim(file)[2])
	{
		a <- file[,j]
		a <- a[-which(is.na(a))]
		a <- a[which(a == "1")]
		sites[count] <- length(a)
		count=count+1
	}
	assign(name,sites)
}

# No. of alt homo sites per individual
for (f in c("az.del.","wtx.del.","priv.az.del.","priv.wtx.del.","shared.del."))
{
	file <- get(paste(f,"geno",sep=""))
	name <- paste(f,"homo",sep="")
	sites <- vector()
	count <- 1
	for (j in 3:dim(file)[2])
	{
		a <- file[,j]
		a <- a[-which(is.na(a))]
		b <- a[which(a == "1")]
		c <- a[which(a == "2")]
		sites[count] <- length(b) + 2*length(c)
		count=count+1
	}
	assign(name,sites)
}


# Get df

az.shared <- shared.del.hetero[1:28]
wtx.shared <- shared.del.hetero[29:59]

df.heter <- rbind(data.frame(Site="AZ",TotalDel=az.del.hetero,PrivDel=priv.az.del.hetero,SharedDel=az.shared),
	data.frame(Site="WTX",TotalDel=wtx.del.hetero,PrivDel=priv.wtx.del.hetero,SharedDel=wtx.shared))

az.shared <- shared.del.homo[1:28]
wtx.shared <- shared.del.homo[29:59]

df.homo <- rbind(data.frame(Site="AZ",TotalDel=az.del.homo,PrivDel=priv.az.del.homo,SharedDel=az.shared),
	data.frame(Site="WTX",TotalDel=wtx.del.homo,PrivDel=priv.wtx.del.homo,SharedDel=wtx.shared))


write.table(df.heter,"Heterozygotes.txt",quote=F,row.names=F)
write.table(df.homo,"AltHomos.txt",quote=F,row.names=F)


# Combine dataset of 
# type of mutation (type)
# No. of sites (snps)
# No. of hetero (het)
# No. of alt homo (althom)
# No.of alleles (all)

df2 <- data.frame(matrix(ncol = 6, nrow = 0))
df3 <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df2) <- c("Site","type","snps","het","althom","all")
colnames(df3) <- c("Site","type","snps","het","althom","all")

df <- data.frame()
		df <- rbind(data.frame(Site=pop,type=mut,snps=get(paste(pop,"_",mut,"_sites",sep="")),
			het=get(paste(pop,"_",mut,"_hetero",sep="")),althom=get(paste(pop,"_",mut,"_homo",sep="")),
			all=get(paste(pop,"_",mut,"_all",sep=""))))

		df2 <- rbind(df2,df)
	}
		df3 <- rbind(df3,df2)

}

df3 <- df3[c(1:84,169:261,439:447,634:645),]



############# Load in ROHs of different length ###################


setwd("/scratch/bell/mathur20/MQU/ch3/revise/load/rohs/geno")

r <- data.frame(from = c("0/0","0/1","1/1","./."), 
                to = c(0,1,2,NA))

SampleID <- read.table("/scratch/bell/mathur20/MQU/ch3/revise/lists/best66.samples.list",header=F)
colnames(SampleID) <- "SampleID"

# Get files

a <- SampleID$SampleID[-20]
a <- "E7969"
for (ind in a)
{
	
	for (len in c("short","medium","long"))
	{
			file <- read.table(paste(ind,".",len,".GT.FORMAT",sep=""),header=T)
			for (j in colnames(file)[3:dim(file)[2]])
			{
				file <- FindReplace(data = file, Var = j, replaceData = r,
                                	from = "from", to = "to", exact = T)
			}
		f <- paste(ind,".",len,sep="")
		assign(f,file)
	}
}

# Get only deleterious/tolerated/nonsynonymous sites

del <- read.table("/scratch/bell/mathur20/MQU/ch3/revise/annotation/snps/all66/best66.deleterious.SNPs",header=F)
tol <- read.table("/scratch/bell/mathur20/MQU/ch3/revise/annotation/snps/all66/best66.tolerated.SNPs",header=F)
fun <- read.table("/scratch/bell/mathur20/MQU/ch3/revise/annotation/snps/all66/best66.nonsynonymous.SNPs",header=F)

colnames(del) <- c("CHROM","POS")
colnames(tol) <- c("CHROM","POS")
colnames(fun) <- c("CHROM","POS")

roh.df <- NULL

for (ind in a)
{
	
	for (len in c("short","medium"))
	{

		file <- get(paste(ind,".",len,sep=""))
		del.df <- inner_join(del,file,by=c("CHROM","POS"))
		tol.df <- inner_join(tol,file,by=c("CHROM","POS"))
		fun.df <- inner_join(fun,file,by=c("CHROM","POS"))

		del.snps <- dim(del.df)[1]
		tol.snps <- dim(tol.df)[1]
		fun.snps <- dim(fun.df)[1]

		del.het <- length(which(del.df[,3] == 1))
		tol.het <- length(which(tol.df[,3] == 1))
		fun.het <- length(which(fun.df[,3] == 1))

		del.alt <- length(which(del.df[,3] == 2))
		tol.alt <- length(which(tol.df[,3] == 2))
		fun.alt <- length(which(fun.df[,3] == 2))

		del.all <- del.het+ 2*del.alt
		tol.all <- tol.het+ 2*tol.alt
		fun.all <- fun.het+ 2*fun.alt

		df <- rbind(data.frame(SampleID=ind, Length=len, Type="Deleterious", 
			Snps=del.snps,Het=del.het,AltHomo=del.alt,All=del.all),
		data.frame(SampleID=ind, Length=len, Type="Tolerated", 
			Snps=tol.snps,Het=tol.het,AltHomo=tol.alt,All=tol.all),
		data.frame(SampleID=ind, Length=len, Type="Functional", 
			Snps=fun.snps,Het=fun.het,AltHomo=fun.alt,All=fun.all))
		roh.df <- rbind(roh.df,df)
	}
}

write.table(roh.df,"../ROhLoad.txt",quote=F,row.names=F)







