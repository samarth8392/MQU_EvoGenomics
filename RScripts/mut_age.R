library(dplyr)


#make batch files

setwd("/scratch/snyder/m/mathur20/MQU/ch3/mut_age/aztwx/sites")

az_sites <- read.table("az.syn.sites.rename.2", header=T)
wtx_sites <- read.table("wtx.syn.sites.rename.2", header=T)


setwd("/scratch/snyder/m/mathur20/MQU/ch3/mut_age/aztwx/batchfiles")

for (i in unique(wtx_sites$Chromosome))
{
	f <- wtx_sites[which(wtx_sites == i),]
	name <- paste("wtx.syn.chr",i,".BATCH",sep="")
	write.csv(f,name, quote=F, row.names=F)
}


setwd("/scratch/snyder/m/mathur20/MQU/ch3/mut_age/run4/out")

for (i in c(1:15,17:28))
{
	name <- paste("best66.chr",i,".marker.txt", sep = "")
	f <- read.table(name,header=T)
	f <- f[,c(1,2,3)]
	assign(name,f)
}

setwd("/scratch/snyder/m/mathur20/MQU/ch3/mut_age/aztwx/weakdel/sites2")

for (i in c(1:10,13:15,17:21,23:27))
{
	#name1 <- paste("az.weakdel.chr",i,".sites2.txt", sep = "")
	name2 <- paste("wtx.weakdel.chr",i,".sites2.txt", sep = "")
	#f1 <- read.table(name1,header=T)
	f2 <- read.table(name2,header=T)
	assign(name1,f1)
	assign(name2,f2)
}

setwd("/scratch/snyder/m/mathur20/MQU/ch3/mut_age/aztwx/final/weakdel")

for (i in c(1:10,13:15,17:21,23:27))
{
	name1 <- paste("best66.chr",i,".marker.txt", sep = "")
	name2 <- paste("az.weakdel.chr",i,".sites2.txt", sep = "")
	name3 <- paste("wtx.weakdel.chr",i,".sites2.txt", sep = "")
	markers <- get(name1)
	az_age <- get(name2)
	wtx_age <- get(name3)

	az_comb <- inner_join(az_age,markers,by="MarkerID")
	wtx_comb <- inner_join(wtx_age,markers,by="MarkerID")
	name4 <- paste("az.weakdel.chr",i,".age",sep="")
	name5 <- paste("wtx.weakdel.chr",i,".age",sep="")
	write.csv(az_comb,name4,quote=F, sep=" ", row.names=F)
	write.csv(wtx_comb,name5,quote=F, sep=" ", row.names=F)
}

setwd("/scratch/snyder/m/mathur20/MQU/ch3/mut_age/aztwx/sites")

#del_sites_shared_azwtx <- read.table("del_sites_shared_azwtx.txt.rename",header=T)
#del_sites_private_az <- read.table("del_sites_private_az.txt.rename",header=T)
#del_sites_private_wtx <- read.table("del_sites_private_wtx.txt.rename",header=T)

del_az_sites <- read.table("az.weakdel.sites.rename.2", header=T)
del_wtx_sites <- read.table("wtx.weakdel.sites.rename.2", header=T)

#weakdel_sites_shared_azwtx <- read.table("weakdel_sites_shared_azwtx.txt",header=T)
#weakdel_sites_private_az <- read.table("weakdel_sites_private_az.txt",header=T)
#weakdel_sites_private_wtx <- read.table("weakdel_sites_private_wtx.txt",header=T)

#syn_sites_shared_azwtx<- read.table("syn_sites_shared_azwtx.txt",header=T)
#syn_sites_private_az <- read.table("syn_sites_private_az.txt",header=T)
#syn_sites_private_wtx <- read.table("syn_sites_private_wtx.txt",header=T)

#for (i in c("del_sites_shared_azwtx","del_sites_private_az","del_sites_private_wtx",
#	"weakdel_sites_shared_azwtx","weakdel_sites_private_az","weakdel_sites_private_wtx",
#	"syn_sites_shared_azwtx","syn_sites_private_az","syn_sites_private_wtx"))
#{
#	f <- get(i)
#	colnames(f) <- c("Chromosome", "Position")
#	assign(i,f)
#	write.table(f,i,quote=F)
#}
del_az_sites$Chromosome <- as.numeric(del_az_sites$Chromosome)
del_wtx_sites$Chromosome <- as.numeric(del_wtx_sites$Chromosome)
#del_sites_private_wtx$Chromosome <- as.numeric(del_sites_private_wtx$Chromosome)



setwd("/scratch/snyder/m/mathur20/MQU/ch3/mut_age/aztwx/final/weakdel")

age1 <- read.csv("az.weakdel.allchr.age", header=T)
age2 <- read.csv("wtx.weakdel.allchr.age", header=T)
#age3 <- read.csv("del.private.wtx.allchr.age", header=T)

az_del_age <- inner_join(del_az_sites,age1, by=c("Chromosome","Position")) 
wtx_del_age <- inner_join(del_wtx_sites,age2, by=c("Chromosome","Position")) 



#shared_del_age <- inner_join(del_sites_shared_azwtx,age1, by=c("Chromosome","Position")) 
#priv_az_del_age <- inner_join(del_sites_private_az,age2, by=c("Chromosome","Position")) 
#priv_wtx_del_age <- inner_join(del_sites_private_wtx,age3, by=c("Chromosome","Position")) 

#shared_del_age <- shared_del_age[which(shared_del_age$Clock == "J"),]
#priv_az_del_age <- priv_az_del_age[which(priv_az_del_age$Clock == "J"),]
#priv_wtx_del_age <- priv_wtx_del_age[which(priv_wtx_del_age$Clock == "J"),]


write.csv(az_del_age,"az_weakdel.age",quote=F, sep=" ", row.names=F)
write.csv(wtx_del_age,"wtx_weakdel.age",quote=F, sep=" ", row.names=F)
#write.csv(priv_wtx_del_age,"priv_wtx_del.age",quote=F, sep=" ", row.names=F)



setwd("/scratch/snyder/m/mathur20/MQU/ch3/mut_age/freq")

del_az_freq <- read.table("az.best66.weakdel.frq.rename.rename2", header=T)
del_wtx_freq <- read.table("wtx.best66.weakdel.frq.rename.rename2", header=T)

az_del_age_freq <- inner_join(az_del_age,del_az_freq, by=c("Chromosome","Position")) 
wtx_del_age_freq <- inner_join(wtx_del_age,del_wtx_freq, by=c("Chromosome","Position")) 

#priv_az_del_age_freq <- inner_join(priv_az_del_age,del_az_freq, by=c("Chromosome","Position")) 


#priv_wtx_del_age_freq <- inner_join(priv_wtx_del_age,del_wtx_freq, by=c("Chromosome","Position")) 


#del_shared_freq <- inner_join(shared_del_age,del_wtx_freq, by=c("Chromosome","Position")) 

#colnames(del_shared_freq)[10] <- "freq_A_wtx"
#colnames(del_shared_freq)[11] <- "freq_B_wtx" 

#del_shared_freq2 <- inner_join(del_shared_freq,del_az_freq, by=c("Chromosome","Position")) 

#colnames(del_shared_freq2)[14] <- "freq_A_az"
#colnames(del_shared_freq2)[15] <- "freq_B_az" 

setwd("/scratch/snyder/m/mathur20/MQU/ch3/mut_age/aztwx/final/weakdel")

write.csv(az_del_age_freq,"az.weakdel.age.withFreq",quote=F, sep=" ", row.names=F)
write.csv(wtx_del_age_freq,"wtx.weakdel.age.withFreq",quote=F, sep=" ", row.names=F)


#write.csv(priv_az_del_age_freq,"priv.az.del.age.withFreq",quote=F, sep=" ", row.names=F)
#write.csv(priv_wtx_del_age_freq,"priv.wtx.del.age.withFreq",quote=F, sep=" ", row.names=F)





setwd("/scratch/snyder/m/mathur20/MQU/ch3/mut_age/sites")

for (i in c("del_sites_shared_azwtx","del_sites_private_az","del_sites_private_wtx",
	"weakdel_sites_shared_azwtx","weakdel_sites_private_az","weakdel_sites_private_wtx",
	"syn_sites_shared_azwtx","syn_sites_private_az","syn_sites_private_wtx"))
{
	name <- paste(i,".txt.rename",sep="")
	f <- read.table(name,header=T)
	for (j in unique(f$Chromosome))
	{
		a <- f[which(f$Chromosome == j),2]
		name2 <- paste(i,"chr",j,".batch",sep="")
		write.csv(a,name2, row.names=F)
	}
}