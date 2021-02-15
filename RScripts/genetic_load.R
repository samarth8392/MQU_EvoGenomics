###########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 09/09/20                  Last Modified: 09/09/20 ###
###########################################################################
###########################################################################
###                   genetic_load.R         		                    ###
###########################################################################

#install.packages("ggpubr", lib="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs", repos="http://cran.us.r-project.org")


library("optparse")

library("plyr")
library("dplyr")
library("ggpubr", lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs" )
library(ggplot2)
library("DataCombine", lib.loc ="/scratch/snyder/m/mathur20/MQU/ch3/Rlibs" )

option_list = list(
  make_option(c("-f","--file"), type="character",default=NULL,help="Geno file",metavar="character"),
  make_option(c("-o","--out"), type="character",default=NULL,help="Output file",metavar="character"),
  make_option(c("-i","--ind"), type="numeric",default=NULL,help="No. of individuals",metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

geno <- read.table(opt$file,header=F)
ind <- opt$ind

for (i in 1:ind)
{
	geno[,i] <- gsub("0|0","0",geno[,i], fixed = T)
	geno[,i] <- gsub("1|0","1",geno[,i], fixed = T)
	geno[,i] <- gsub("0|1","1",geno[,i], fixed = T)
	geno[,i] <- gsub("1|1","0",geno[,i], fixed = T)
}

for (i in 1:ind)
{
	geno[,i] <- as.numeric(geno[,i])
}

n_all <- colSums(geno)
pro_all <- n_all/(2*dim(geno)[1])

#cat("No.of derived alleles per individual",file=opt$out,sep="\n")
#cat(n_all,file=opt$out,append=TRUE,sep="\n")
#cat("Proportion of derived alleles per individual",file=opt$out,append=TRUE,sep="\n")
#cat(pro_all,file=opt$out,append=TRUE,sep="\n")

#Hetero
hom<- vector()
het <- vector()
total <- vector()

for (i in 1:ind)
{
	het[i] <- as.numeric(table(geno[,i])[2])
	hom[i] <- as.numeric(table(geno[,i])[1])
	total[i] <- hom[i] + het[i]
}


cat("Proportion of homozygotes",file=opt$out,sep="\n")
cat(hom/total,file=opt$out,append=TRUE,sep="\n")
cat("Proportion of heterozygotes",file=opt$out,append=TRUE,sep="\n")
cat(het/total,file=opt$out,append=TRUE,sep="\n")

################################################################

######FINAL########

setwd("/scratch/snyder/m/mathur20/MQU/ch3/load/geno")

az_del_s <- read.table("az.shared.del.GT.FORMAT",header=T)
wtx_del_s <- read.table("wtx.shared.del.GT.FORMAT",header=T)

az_wdel_s <- read.table("az.shared.weakdel.GT.FORMAT",header=T)
wtx_wdel_s <- read.table("wtx.shared.weakdel.GT.FORMAT",header=T)

az_syn_s <- read.table("az.shared.syn.GT.FORMAT",header=T)
wtx_syn_s <- read.table("wtx.shared.syn.GT.FORMAT",header=T)


az_del_p <- read.table("az.private.del.GT.FORMAT",header=T)
wtx_del_p <- read.table("wtx.private.del.GT.FORMAT",header=T)

az_wdel_p <- read.table("az.private.weakel.GT.FORMAT",header=T)
wtx_wdel_p <- read.table("wtx.private.weakdel.GT.FORMAT",header=T)

az_syn_p <- read.table("az.private.syn.GT.FORMAT",header=T)
wtx_syn_p <- read.table("wtx.private.syn.GT.FORMAT",header=T)



#Convert genotypes

#for (i in c("az_del","wtx_del","az_weakdel","wtx_weakdel", "az_syn", "wtx_syn"))
for (i in c("az_del_s","wtx_del_s","az_wdel_s","wtx_wdel_s","az_syn_s" ,"wtx_syn_s",
	"az_del_p","wtx_del_p","az_wdel_p","wtx_wdel_p","az_syn_p","wtx_syn_p"))
{
	file <- get(i)
	r <- data.frame(from = c("0/0","0/1","1/1","./."), 
                to = c(0,1,2,NA))

	for (j in colnames(file)[3:dim(file)[2]])
	{
		file <- FindReplace(data = file, Var = j, replaceData = r,
                                from = "from", to = "to", exact = T)
	}
	assign(i,file)
}

#No. of sites

#for (i in c("az_del","wtx_del","az_weakdel","wtx_weakdel", "az_syn", "wtx_syn"))
for (i in c("az_del_s","wtx_del_s","az_wdel_s","wtx_wdel_s","az_syn_s" ,"wtx_syn_s",
	"az_del_p","wtx_del_p","az_wdel_p","wtx_wdel_p","az_syn_p","wtx_syn_p"))
{
	file <- get(i)
	name <- paste(i,"_sites",sep="")
	sites <- vector()

	count=1
	for (j in 3:dim(file)[2])
	{
		a <- file[,j]
		a <- a[-which(is.na(a))]
		#a <- a[which(a == "1" | a == "2")]
		#a <- a[which(a == "1" | a == "2")]
		sites[count] <- length(a)
		count=count+1
	}
	assign(name,sites)
}

# No. of hetero


#for (i in c("az_del","wtx_del","az_weakdel","wtx_weakdel", "az_syn", "wtx_syn"))
for (i in c("az_del_s","wtx_del_s","az_wdel_s","wtx_wdel_s","az_syn_s" ,"wtx_syn_s",
	"az_del_p","wtx_del_p","az_wdel_p","wtx_wdel_p","az_syn_p","wtx_syn_p"))
{
	file <- get(i)
	name <- paste(i,"_alt",sep="")
	sites <- vector()

	count=1
	for (j in 3:dim(file)[2])
	{
		a <- file[,j]
		a <- a[-which(is.na(a))]
		a <- a[which(a == "1")]
		#a <- a[which(a == "1" | a == "2")]
		sites[count] <- length(a)
		count=count+1
	}
	assign(name,sites)
}

# ALt homo

#for (i in c("az_del","wtx_del","az_weakdel","wtx_weakdel", "az_syn", "wtx_syn"))
for (i in c("az_del_s","wtx_del_s","az_wdel_s","wtx_wdel_s","az_syn_s" ,"wtx_syn_s",
	"az_del_p","wtx_del_p","az_wdel_p","wtx_wdel_p","az_syn_p","wtx_syn_p"))
{
	file <- get(i)
	name <- paste(i,"_hom",sep="")
	sites <- vector()

	count=1
	for (j in 3:dim(file)[2])
	{
		a <- file[,j]
		a <- a[-which(is.na(a))]
		a <- a[which(a == "2")]
		#a <- a[which(a == "1" | a == "2")]
		sites[count] <- length(a)
		count=count+1
	}
	assign(name,sites)
}

# No.of alleles

#for (i in c("az_del","wtx_del","az_weakdel","wtx_weakdel", "az_syn", "wtx_syn"))
for (i in c("az_del_s","wtx_del_s","az_wdel_s","wtx_wdel_s","az_syn_s" ,"wtx_syn_s",
	"az_del_p","wtx_del_p","az_wdel_p","wtx_wdel_p","az_syn_p","wtx_syn_p"))
{
	file <- get(i)
	name <- paste(i,"_all",sep="")
	sites <- vector()

	count=1
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

# Combine
del_s_sites <- rbind(data.frame(Site="AZ", value = az_del_s_sites),
data.frame(Site="WTX", value = wtx_del_s_sites))

weakdel_s_sites <- rbind(data.frame(Site="AZ", value = az_wdel_s_sites),
data.frame(Site="WTX", value = wtx_wdel_s_sites))

syn_s_sites <- rbind(data.frame(Site="AZ", value = az_syn_s_sites),
data.frame(Site="WTX", value = wtx_syn_s_sites))

del_p_sites <- rbind(data.frame(Site="AZ", value = az_del_p_sites),
data.frame(Site="WTX", value = wtx_del_p_sites))

weakdel_p_sites <- rbind(data.frame(Site="AZ", value = az_wdel_p_sites),
data.frame(Site="WTX", value = wtx_wdel_p_sites))

syn_p_sites <- rbind(data.frame(Site="AZ", value = az_syn_p_sites),
data.frame(Site="WTX", value = wtx_syn_p_sites))

del_s_alt <- rbind(data.frame(Site="AZ", value = az_del_s_alt),
data.frame(Site="WTX", value = wtx_del_s_alt))

weakdel_s_alt <- rbind(data.frame(Site="AZ", value = az_wdel_s_alt),
data.frame(Site="WTX", value = wtx_wdel_s_alt))

syn_s_alt <- rbind(data.frame(Site="AZ", value = az_syn_s_alt),
data.frame(Site="WTX", value = wtx_syn_s_alt))

del_p_alt <- rbind(data.frame(Site="AZ", value = az_del_p_alt),
data.frame(Site="WTX", value = wtx_del_p_alt))

weakdel_p_alt <- rbind(data.frame(Site="AZ", value = az_wdel_p_alt),
data.frame(Site="WTX", value = wtx_wdel_p_alt))

syn_p_alt <- rbind(data.frame(Site="AZ", value = az_syn_p_alt),
data.frame(Site="WTX", value = wtx_syn_p_alt))

del_s_hom <- rbind(data.frame(Site="AZ", value = az_del_s_hom),
data.frame(Site="WTX", value = wtx_del_s_hom))

weakdel_s_hom <- rbind(data.frame(Site="AZ", value = az_wdel_s_hom),
data.frame(Site="WTX", value = wtx_wdel_s_hom))

syn_s_hom <- rbind(data.frame(Site="AZ", value = az_syn_s_hom),
data.frame(Site="WTX", value = wtx_syn_s_hom))

del_p_hom <- rbind(data.frame(Site="AZ", value = az_del_p_hom),
data.frame(Site="WTX", value = wtx_del_p_hom))

weakdel_p_hom <- rbind(data.frame(Site="AZ", value = az_wdel_p_hom),
data.frame(Site="WTX", value = wtx_wdel_p_hom))

syn_p_hom <- rbind(data.frame(Site="AZ", value = az_syn_p_hom),
data.frame(Site="WTX", value = wtx_syn_p_hom))

del_s_all <- rbind(data.frame(Site="AZ", value = az_del_s_all),
data.frame(Site="WTX", value = wtx_del_s_all))

weakdel_s_all <- rbind(data.frame(Site="AZ", value = az_wdel_s_all),
data.frame(Site="WTX", value = wtx_wdel_s_all))

syn_s_all <- rbind(data.frame(Site="AZ", value = az_syn_s_all),
data.frame(Site="WTX", value = wtx_syn_s_all))

del_p_all <- rbind(data.frame(Site="AZ", value = az_del_p_all),
data.frame(Site="WTX", value = wtx_del_p_all))

weakdel_p_all <- rbind(data.frame(Site="AZ", value = az_wdel_p_all),
data.frame(Site="WTX", value = wtx_wdel_p_all))

syn_p_all <- rbind(data.frame(Site="AZ", value = az_syn_p_all),
data.frame(Site="WTX", value = wtx_syn_p_all))

#Output files

#SNPS ((del_sites, wdel_sites, syn_sites))
# Het (del_alt, wdel_alt, syn_alt)
# Hom (del_hom, wdel_hom, syn_hom)
# Alleles (del_all, wdel_all, syn_all)

p1 <- cbind(del_s_sites,weakdel_s_sites,syn_s_sites,
del_p_sites,weakdel_p_sites,syn_p_sites,
del_s_alt,weakdel_s_alt,syn_s_alt,
del_p_alt,weakdel_p_alt,syn_p_alt,
del_s_hom,weakdel_s_hom,syn_s_hom,
del_p_hom,weakdel_p_hom,syn_p_hom,
del_s_all,weakdel_s_all,syn_s_all,
del_p_all,weakdel_p_all,syn_p_all)



# SNPs (shared_az, shared_wtx, priv_az, priv_wtx)
# Het
# Hom
# ALleles
# Potential Load
# Realized Load 


#p1 <- p1[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36)]
p1 <- p1[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46, 48)]

#colnames(p1) <- c("Site","Del_sites","Weakdel_sites","Syn_sites",
#	"Del_het","Weakdel_het","Syn_het",
#	"Del_homo","Weakdel_homo","Syn_homo",
#	"Del_alleles","Weakdel_alleles","Syn_alleles",
#	"Del_potential","Weakdel_potential","Syn_potential",
#	"Del_realized","Weakdel_realized","Syn_realized")


colnames(p1) <- c("Site","del_s_sites","weakdel_s_sites","syn_s_sites",
"del_p_sites","weakdel_p_sites","syn_p_sites",
"del_s_alt","weakdel_s_alt","syn_s_alt",
"del_p_alt","weakdel_p_alt","syn_p_alt",
"del_s_hom","weakdel_s_hom","syn_s_hom",
"del_p_hom","weakdel_p_hom","syn_p_hom",
"del_s_all","weakdel_s_all","syn_s_all",
"del_p_all","weakdel_p_all","syn_p_all")

write.csv(p1,"Individual_Shared_Private_load.csv",quote=F,row.names=F)



#az_del_s <- data.frame(Site="AZ", value = az_del_sites)
#wtx_del_s <- data.frame(Site="WTX", value = wtx_del_sites)

#az_wdel_s <- data.frame(Site="AZ", value = az_weakdel_sites)
#wtx_wdel_s <- data.frame(Site="WTX", value = wtx_weakdel_sites)

#az_syn_s <- data.frame(Site="AZ", value = az_syn_sites)
#wtx_syn_s <- data.frame(Site="WTX", value = wtx_syn_sites)


# Combine
wtx_del_shared_s <- data.frame(Site="WTX", value = wtx_del_shared_sites)

az_del_priv_s <- data.frame(Site="AZ", value = az_del_private_sites)
wtx_del_priv_s <- data.frame(Site="WTX", value = wtx_del_private_sites)


del_shared_sites <- rbind(az_del_shared_s,wtx_del_shared_s)
mu_del_shared_sites <- ddply(del_shared_sites, "Site", summarise, grp.mean=mean(value))
sd_del_shared_sites <- ddply(del_shared_sites, "Site", summarise, grp.mean=sd(value))

del_private_sites <- rbind(az_del_priv_s,wtx_del_priv_s)
mu_del_shared_sites <- ddply(del_private_sites, "Site", summarise, grp.mean=mean(value))
sd_del_shared_sites <- ddply(del_private_sites, "Site", summarise, grp.mean=sd(value))

#wdel_sites <- rbind(az_wdel_s,wtx_wdel_s)
#mu_wdel_sites <- ddply(wdel_sites, "Site", summarise, grp.mean=mean(value))
#sd_wdel_sites <- ddply(wdel_sites, "Site", summarise, grp.mean=sd(value))

#syn_sites <- rbind(az_syn_s,wtx_syn_s)
#mu_syn_sites <- ddply(syn_sites, "Site", summarise, grp.mean=mean(value))
#sd_syn_sites <- ddply(syn_sites, "Site", summarise, grp.mean=sd(value))





my_comparisons <- list( c("AZ", "WTX"))

az_del_s_a <- data.frame(Site="AZ", value = az_del_shared_alt)
wtx_del_s_a <- data.frame(Site="WTX", value = wtx_del_shared_alt)

az_del_p_a <- data.frame(Site="AZ", value = az_del_private_alt)
wtx_del_p_a <- data.frame(Site="WTX", value = wtx_del_private_alt)

#az_syn_a <- data.frame(Site="AZ", value = az_syn_alt)
#wtx_syn_a <- data.frame(Site="WTX", value = wtx_syn_alt)


del_s_alt <- rbind(az_del_s_a,wtx_del_s_a)
mu_del_alt <- ddply(del_s_alt, "Site", summarise, grp.mean=mean(value))
sd_del_alt <- ddply(del_s_alt, "Site", summarise, grp.mean=sd(value))

del_p_alt <- rbind(az_del_p_a,wtx_del_p_a)
mu_del_alt <- ddply(del_p_alt, "Site", summarise, grp.mean=mean(value))
sd_del_alt <- ddply(del_p_alt, "Site", summarise, grp.mean=sd(value))


#wdel_alt <- rbind(az_wdel_a,wtx_wdel_a)
#mu_wdel_alt <- ddply(wdel_alt, "Site", summarise, grp.mean=mean(value))
#sd_wdel_alt <- ddply(wdel_alt, "Site", summarise, grp.mean=sd(value))

#syn_alt <- rbind(az_syn_a,wtx_syn_a)
#mu_syn_alt <- ddply(syn_alt, "Site", summarise, grp.mean=mean(value))
#sd_syn_alt <- ddply(syn_alt, "Site", summarise, grp.mean=sd(value))


# Alt homo


my_comparisons <- list( c("AZ", "WTX"))

az_del_s_h <- data.frame(Site="AZ", value = az_del_shared_hom)
wtx_del_s_h <- data.frame(Site="WTX", value = wtx_del_shared_hom)

az_del_p_h <- data.frame(Site="AZ", value = az_del_private_hom)
wtx_del_p_h <- data.frame(Site="WTX", value = wtx_del_private_hom)

#az_wdel_h <- data.frame(Site="AZ", value = az_weakdel_hom)
#wtx_wdel_h <- data.frame(Site="WTX", value = wtx_weakdel_hom)

#az_syn_h <- data.frame(Site="AZ", value = az_syn_hom)
#wtx_syn_h <- data.frame(Site="WTX", value = wtx_syn_hom)


del_homo_s <- rbind(az_del_s_h,wtx_del_s_h)
mu_del_homo <- ddply(del_homo_s, "Site", summarise, grp.mean=mean(value))
sd_del_homo <- ddply(del_homo_s, "Site", summarise, grp.mean=sd(value))

del_homo_p <- rbind(az_del_p_h,wtx_del_p_h)
mu_del_homo <- ddply(del_homo_p, "Site", summarise, grp.mean=mean(value))
sd_del_homo <- ddply(del_homo_p, "Site", summarise, grp.mean=sd(value))

#wdel_homo <- rbind(az_wdel_h,wtx_wdel_h)
#mu_wdel_homo <- ddply(wdel_homo, "Site", summarise, grp.mean=mean(value))
#sd_wdel_homo <- ddply(wdel_homo, "Site", summarise, grp.mean=sd(value))

#syn_homo <- rbind(az_syn_h,wtx_syn_h)
#mu_syn_homo <- ddply(syn_homo, "Site", summarise, grp.mean=mean(value))
#sd_syn_homo <- ddply(syn_homo, "Site", summarise, grp.mean=sd(value))




my_comparisons <- list( c("AZ", "WTX"))

az_del_s_al <- data.frame(Site="AZ", value = az_del_shared_all)
wtx_del_s_al <- data.frame(Site="WTX", value = wtx_del_shared_all)

az_del_p_al <- data.frame(Site="AZ", value = az_del_private_all)
wtx_del_p_al <- data.frame(Site="WTX", value = wtx_del_private_all)

#az_wdel_al <- data.frame(Site="AZ", value = az_weakdel_all)
#wtx_wdel_al <- data.frame(Site="WTX", value = wtx_weakdel_all)

#az_syn_al <- data.frame(Site="AZ", value = az_syn_all)
#wtx_syn_al <- data.frame(Site="WTX", value = wtx_syn_all)


del_s_all <- rbind(az_del_s_al,wtx_del_s_al)
mu_del_all <- ddply(del_s_all, "Site", summarise, grp.mean=mean(value))
sd_del_all <- ddply(del_s_all, "Site", summarise, grp.mean=sd(value))

del_p_all <- rbind(az_del_p_al,wtx_del_p_al)
mu_del_all <- ddply(del_p_all, "Site", summarise, grp.mean=mean(value))
sd_del_all <- ddply(del_p_all, "Site", summarise, grp.mean=sd(value))

#wdel_all <- rbind(az_wdel_al,wtx_wdel_al)
#mu_wdel_all <- ddply(wdel_all, "Site", summarise, grp.mean=mean(value))
#sd_wdel_all <- ddply(wdel_all, "Site", summarise, grp.mean=sd(value))

#syn_all <- rbind(az_syn_al,wtx_syn_al)
#mu_syn_all <- ddply(syn_all, "Site", summarise, grp.mean=mean(value))
#sd_syn_all <- ddply(syn_all, "Site", summarise, grp.mean=sd(value))



# Potential load = # no. of deleterious allele / Total no. of SNPs

az_del_s_pot <- data.frame(Site="AZ", value = az_del_shared_all/az_del_shared_sites)
wtx_del_s_pot <- data.frame(Site="WTX", value = wtx_del_shared_all/wtx_del_shared_sites)

az_del_p_pot <- data.frame(Site="AZ", value = az_del_private_all/az_del_private_sites)
wtx_del_p_pot <- data.frame(Site="WTX", value = wtx_del_private_all/wtx_del_private_sites)

#az_wdel_pot <- data.frame(Site="AZ", value = az_weakdel_all/az_weakdel_sites)
#wtx_wdel_pot <- data.frame(Site="WTX", value = wtx_weakdel_all/wtx_weakdel_sites)

#az_syn_pot <- data.frame(Site="AZ", value = az_syn_all/az_syn_sites)
#wtx_syn_pot <- data.frame(Site="WTX", value = wtx_syn_all/wtx_syn_sites)

del_s_pot <- rbind(az_del_s_pot,wtx_del_s_pot)
mu_del_pot <- ddply(del_s_pot, "Site", summarise, grp.mean=mean(value))
sd_del_pot <- ddply(del_s_pot, "Site", summarise, grp.mean=sd(value))

del_p_pot <- rbind(az_del_p_pot,wtx_del_p_pot)
mu_del_pot <- ddply(del_p_pot, "Site", summarise, grp.mean=mean(value))
sd_del_pot <- ddply(del_p_pot, "Site", summarise, grp.mean=sd(value))



#wdel_pot <- rbind(az_wdel_pot,wtx_wdel_pot)
#mu_wdel_pot <- ddply(wdel_pot, "Site", summarise, grp.mean=mean(value))
#sd_wdel_pot <- ddply(wdel_pot, "Site", summarise, grp.mean=sd(value))

#syn_pot <- rbind(az_syn_pot,wtx_syn_pot)
#mu_syn_pot <- ddply(syn_pot, "Site", summarise, grp.mean=mean(value))
#sd_syn_pot <- ddply(syn_pot, "Site", summarise, grp.mean=sd(value))

# Realized load = # no. of alt hom / no. of deleterious allele

az_del_s_real <- data.frame(Site="AZ", value = az_del_shared_hom/az_del_shared_all)
wtx_del_s_real <- data.frame(Site="WTX", value = wtx_del_shared_hom/wtx_del_shared_all)

az_del_p_real <- data.frame(Site="AZ", value = az_del_private_hom/az_del_private_all)
wtx_del_p_real <- data.frame(Site="WTX", value = wtx_del_private_hom/wtx_del_private_all)

#az_wdel_real <- data.frame(Site="AZ", value = az_weakdel_hom/az_weakdel_all)
#wtx_wdel_real <- data.frame(Site="WTX", value = wtx_weakdel_hom/wtx_weakdel_all)

#az_syn_real <- data.frame(Site="AZ", value = az_syn_hom/az_syn_all)
#wtx_syn_real <- data.frame(Site="WTX", value = wtx_syn_hom/wtx_syn_all)

del_s_real <- rbind(az_del_s_real,wtx_del_s_real)
mu_del_real <- ddply(del_s_real, "Site", summarise, grp.mean=mean(value))
sd_del_real <- ddply(del_s_real, "Site", summarise, grp.mean=sd(value))


del_p_real <- rbind(az_del_p_real,wtx_del_p_real)
mu_del_real <- ddply(del_p_real, "Site", summarise, grp.mean=mean(value))
sd_del_real <- ddply(del_p_real, "Site", summarise, grp.mean=sd(value))
#wdel_real <- rbind(az_wdel_real,wtx_wdel_real)
#mu_wdel_real <- ddply(wdel_real, "Site", summarise, grp.mean=mean(value))
#sd_wdel_real <- ddply(wdel_real, "Site", summarise, grp.mean=sd(value))

#syn_real <- rbind(az_syn_real,wtx_syn_real)
#mu_syn_real <- ddply(syn_real, "Site", summarise, grp.mean=mean(value))
#sd_syn_real <- ddply(syn_real, "Site", summarise, grp.mean=sd(value))




# Plot Results
#SNP sites comparison (del_sites, wdel_sites, syn_sites)

p1<- ggplot(az_roh, aes(x=id, y=lengthBps, fill=id))+
  geom_violin()+
  geom_boxplot(width=0.1)+
  theme_classic(base_size = 15)+ scale_fill_manual(values=distinctColorPalette(28))+
  labs(x="Site", y = "Distribution of ROH length (kBp)")+
  theme(axis.title.x=element_blank(),
        legend.title=element_blank(),
        axis.text.x =element_text(angle = 45,vjust = 0.75),
        legend.position = "none")