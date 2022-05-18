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

# Load colors

col_az <- "#CC3A8E" #carto_pal(12,"Vivid")[11]
col_wtx <- "#9f82ce" #carto_pal(7,"Purp")[5]
col_etx <-"#f2855d" #carto_pal(7,"Peach")[5]
col_mx <- "#3598FF"
cols <- c(col_az,col_wtx,col_etx, col_mx)
colss <-c(rep(col_az,28),rep(col_mx,4),rep(col_wtx,31),rep(col_etx,3))
cols_l <- adjustcolor(cols, alpha.f = 0.5)
col_azwtx <- "#90CCC2"

############# FOR POPULATION #######################
# Estimate the number of heterozygous sites in each population per 1kb


# Get number of heterozygote individuals at each SNP in each population

setwd("/scratch/bell/mathur20/MQU/ch3/revise/hetperkb/bypop/")
#pops <- c("az","wtx","etx","mx")
#for (chr in c(1:33))
#{
#	a <- read.table(paste("chr",chr,"/","az.best66.hwe",sep=""), header=T)
#	a <- a[,c(1,2)]
#	for (i in c(1:4))
#	{
#		het <- read.table(paste("chr",chr,"/",pops[i],".best66.het",sep=""), header=T)
#		a <- cbind(a,het)
#		colnames(a)[i+2] <- pops[i]
#	}
#	write.table(a,paste("chr",chr,"/","allpop.het",sep=""),quote=F,row.names=F)
#}
#
# Get no. of heterozygous sites in each 1kb bin for each chromosome

hetsite = function(x){
	p <- length(which(x !=0))
    return(p)
}

means <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(means) <- c("CHROM","Mean_AZ","SD_AZ","Mean_WTX","SD_WTX","Mean_ETX","SD_ETX","Mean_MX","SD_MX")

pvals <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(pvals) <- c("AZ-WTX","AZ-ETX","AZ-MX","WTX-ETX","WTX-MX","ETX-MX")

count <- 0
for (chr in c(1:28,30:33))
{
	count <- count + 1
	bins <- read.table(paste("../bins/bins100kb/chr",chr,"_bin100kb.txt",sep=""),header=T)
	het <- read.table(paste("chr",chr,"/allpop.het",sep=""),header=T)
	colnames(het)[1] <- "CHROM"
	df <- data.frame(matrix(ncol = 7, nrow = 0))
	colnames(df) <- c("CHROM","START","END","AZ","WTX","ETX","MX")

	for (i in 1:nrow(bins))
	{
		a <- het[which(het$POS > bins$START[i] & het$POS < bins$END[i]),]
		hets <- sapply(a[3:6],hetsite)
		df[i,] <- c(bins$CHROM[i],bins$START[i],bins$END[i],hets[1],hets[2],hets[3],hets[4])
	}
	df <- na.omit(df)
	#write.table(df,paste("chr",chr,"/","allpop.hetper100kb",sep=""),quote=F,row.names=F)
	
	# Draw barplot
	pdf(file = paste("../plots/chr",chr,".het100kb.pdf",sep = ""), width = 12, height = 12)
  	par(mfrow=c(4,1))
  	barplot(df$AZ/1e3,col=col_az, border=col_az,xlab="Position (Mb)", 
  		ylab="Het/kb", ylim=c(0,2), main="AZ", yaxt = "n")
  	axis(2,at=seq(0,2,0.5))
  	abline(h=mean(df$AZ/1e3),lwd=2,col=col_az)
  	
  	barplot(df$WTX/1e3,col=col_wtx,border=col_wtx,xlab="Position (Mb)", 
  		ylab="Het/kb",ylim=c(0,2), main="WTX", yaxt = "n")
  	axis(2,at=seq(0,2,0.5))
  	abline(h=mean(df$WTX/1e3),lwd=2,col=col_wtx)

  	barplot(df$ETX/1e3,col=col_etx,border=col_etx, xlab="Position (Mb)", 
  		ylab="Het/kb",ylim=c(0,2), main="CTX", yaxt = "n")
  	axis(2,at=seq(0,2,0.5))
  	abline(h=mean(df$ETX/1e3),lwd=2,col=col_etx)

  	barplot(df$MX/1e3,col=col_mx,border=col_mx,xlab="Position (Mb)", 
  		ylab="Het/kb",ylim=c(0,2), main="MX", yaxt = "n")
  	axis(2,at=seq(0,2,0.5))
  	abline(h=mean(df$MX/1e3),lwd=2,col=col_mx)

  	dev.off()

  	#Draw histogram
  	#pdf(file = paste("../plots/chr",chr,".het100kb.hist.pdf",sep = ""), width = 12, height = 12)
  	#par(mfrow=c(4,1))
  	#h1 <- hist(df$AZ/1e3, plot=F)
  	#h2 <- hist(df$WTX/1e3, plot=F)
  	#h3 <- hist(df$ETX/1e3, plot=F)
  	#h4 <- hist(df$MX/1e3, plot=F)
  	#ylim=max(h1$counts,h2$counts,h3$counts,h4$counts)+50
  	#xlim=max(h1$breaks,h2$breaks,h3$breaks,h4$breaks)+0.5

  	#hist(df$AZ/1e3,col=col_az, ylab="No. of windows", xlab="Het/kb", main="AZ",ylim=c(0,ylim),xlim=c(0,xlim))
  	#abline(v=mean(df$AZ/1e3), lty=2, lwd=4,col=col_az)
  	#hist(df$WTX/1e3,col=col_wtx,ylab="No. of windows", xlab="Het/kb", main="WTX",ylim=c(0,ylim),xlim=c(0,xlim))
  	#abline(v=mean(df$WTX/1e3), lty=2, lwd=4,col=col_wtx)
  	#hist(df$ETX/1e3,col=col_etx,ylab="No. of windows", xlab="Het/kb", main="CTX",ylim=c(0,ylim),xlim=c(0,xlim))
  	#abline(v=mean(df$ETX/1e3), lty=2, lwd=4,col=col_etx)
  	#hist(df$MX/1e3,col=col_mx,ylab="No. of windows", xlab="Het/kb", main="MX",ylim=c(0,ylim),xlim=c(0,xlim))
  	#abline(v=mean(df$MX/1e3), lty=2, lwd=4,col=col_mx)
  	#dev.off()

  	#Get means and compare 
  	#means[count,] <- c(paste("chr",chr,sep=""),mean(df$AZ/1e3),sd(df$AZ/1e3),
  	#	mean(df$WTX/1e3),sd(df$WTX/1e3),mean(df$ETX/1e3),sd(df$ETX/1e3),
  	#	mean(df$MX/1e3),sd(df$MX/1e3))

  	#df2 <- rbind(data.frame(Site="AZ",value=df$AZ),
  	#	data.frame(Site="WTX",value=df$WTX),
  	#	data.frame(Site="ETX",value=df$ETX),
  	#	data.frame(Site="MX",value=df$MX))

  	#test <- pairwise.wilcox.test(df2$value, df2$Site,
    #                 p.adjust.method = "fdr")
  	#pvals[count,] <- c(paste("chr",chr,sep=""),as.numeric(na.omit(test$p.value[1:9])))



}

write.table(means,"../Meanhetperkb.txt",quote=F,row.names=F)
write.table(pvals,"../pValshetperkb.txt",quote=F,row.names=F)


############# FOR INDIVIDUAL #######################


# Estimate the number of heterozygous sites in each individual per 1kb


# Get number of heterozygote individuals at each SNP in each population

setwd("/scratch/bell/mathur20/MQU/ch3/revise/hetperkb/byindiv/")

inds <- read.table("/scratch/bell/mathur20/MQU/ch3/revise/lists/best66.samples.list", header=F)
inds <- as.array(inds[,1])
for (chr in c(1:33))
{
	a <- read.table(paste("chr",chr,"/",inds[1],".best66.hwe",sep=""), header=T)
	a <- a[,c(1,2)]
	for (i in 1:length(inds))
	{
		het <- read.table(paste("chr",chr,"/",inds[i],".best66.het",sep=""), header=T)
		a <- cbind(a,het)
		colnames(a)[i+2] <- as.character(inds[i])
	}
	write.table(a,paste("chr",chr,"/","allinds.het",sep=""),quote=F,row.names=F)
}


hetsite = function(x){
	p <- length(which(x !=0))
    return(p)
}

means <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(means) <- c("CHROM","Mean_AZ","SD_AZ","Mean_WTX","SD_WTX","Mean_ETX","SD_ETX","Mean_MX","SD_MX")

pvals <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(pvals) <- c("AZ-WTX","AZ-ETX","AZ-MX","WTX-ETX","WTX-MX","ETX-MX")


cols <- c(rep(col_az,28),rep(col_mx,4),rep(col_wtx,31),rep(col_etx,3))

count <- 0
for (chr in c(1:28,30:33))
{
	count <- count + 1
	bins <- read.table(paste("../bins/bins100kb/chr",chr,"_bin100kb.txt",sep=""),header=T)
	het <- read.table(paste("chr",chr,"/allinds.het",sep=""),header=T)
	colnames(het)[1] <- "CHROM"
	df <- data.frame(matrix(ncol = 69, nrow = 0))
	colnames(df) <- c("CHROM","START","END",colnames(het)[3:68])

	for (i in 1:nrow(bins))
	{
		a <- het[which(het$POS > bins$START[i] & het$POS < bins$END[i]),]
		hets <- sapply(a[3:68],hetsite)
		df[i,] <- c(bins$CHROM[i],bins$START[i],bins$END[i],hets[c(1:66)])
	}
	df <- na.omit(df)
	write.table(df,paste("chr",chr,"/","allind.hetper100kb",sep=""),quote=F,row.names=F)
	
	# Draw barplot
	for (i in 4:69)
	{
		pdf(file = paste("../plots/byindiv/chr",chr,".",colnames(df)[i],".het100kb.pdf",sep = ""), width = 9, height = 6)
		barplot(df[,i]/1e3,col=cols[i-3], border=cols[i-3],xlab="Position (Mb)", 
  		ylab="Het/kb", ylim=c(0,2), main=colnames(df)[i], yaxt = "n")
  		axis(2,at=seq(0,2,0.5))
  		#abline(h=mean(df[,i]/1e3),lwd=2,col=cols[i-3])
  		dev.off()
	}
	
  	par(mfrow=c(4,1))
  	barplot(df$AZ/1e3,col=col_az, border=col_az,xlab="Position (Mb)", 
  		ylab="Het/kb", ylim=c(0,2), main="AZ", yaxt = "n")
  	axis(2,at=seq(0,2,0.5))
  	abline(h=mean(df$AZ/1e3),lwd=2,col=col_az)
  	
  	barplot(df$WTX/1e3,col=col_wtx,border=col_wtx,xlab="Position (Mb)", 
  		ylab="Het/kb",ylim=c(0,2), main="WTX", yaxt = "n")
  	axis(2,at=seq(0,2,0.5))
  	abline(h=mean(df$WTX/1e3),lwd=2,col=col_wtx)

  	barplot(df$ETX/1e3,col=col_etx,border=col_etx, xlab="Position (Mb)", 
  		ylab="Het/kb",ylim=c(0,2), main="CTX", yaxt = "n")
  	axis(2,at=seq(0,2,0.5))
  	abline(h=mean(df$ETX/1e3),lwd=2,col=col_etx)

  	barplot(df$MX/1e3,col=col_mx,border=col_mx,xlab="Position (Mb)", 
  		ylab="Het/kb",ylim=c(0,2), main="MX", yaxt = "n")
  	axis(2,at=seq(0,2,0.5))
  	abline(h=mean(df$MX/1e3),lwd=2,col=col_mx)

  	dev.off()

  	#Draw histogram
  	#pdf(file = paste("../plots/chr",chr,".het100kb.hist.pdf",sep = ""), width = 12, height = 12)
  	#par(mfrow=c(4,1))
  	#h1 <- hist(df$AZ/1e3, plot=F)
  	#h2 <- hist(df$WTX/1e3, plot=F)
  	#h3 <- hist(df$ETX/1e3, plot=F)
  	#h4 <- hist(df$MX/1e3, plot=F)
  	#ylim=max(h1$counts,h2$counts,h3$counts,h4$counts)+50
  	#xlim=max(h1$breaks,h2$breaks,h3$breaks,h4$breaks)+0.5

  	#hist(df$AZ/1e3,col=col_az, ylab="No. of windows", xlab="Het/kb", main="AZ",ylim=c(0,ylim),xlim=c(0,xlim))
  	#abline(v=mean(df$AZ/1e3), lty=2, lwd=4,col=col_az)
  	#hist(df$WTX/1e3,col=col_wtx,ylab="No. of windows", xlab="Het/kb", main="WTX",ylim=c(0,ylim),xlim=c(0,xlim))
  	#abline(v=mean(df$WTX/1e3), lty=2, lwd=4,col=col_wtx)
  	#hist(df$ETX/1e3,col=col_etx,ylab="No. of windows", xlab="Het/kb", main="CTX",ylim=c(0,ylim),xlim=c(0,xlim))
  	#abline(v=mean(df$ETX/1e3), lty=2, lwd=4,col=col_etx)
  	#hist(df$MX/1e3,col=col_mx,ylab="No. of windows", xlab="Het/kb", main="MX",ylim=c(0,ylim),xlim=c(0,xlim))
  	#abline(v=mean(df$MX/1e3), lty=2, lwd=4,col=col_mx)
  	#dev.off()

  	#Get means and compare 
  	#means[count,] <- c(paste("chr",chr,sep=""),mean(df$AZ/1e3),sd(df$AZ/1e3),
  	#	mean(df$WTX/1e3),sd(df$WTX/1e3),mean(df$ETX/1e3),sd(df$ETX/1e3),
  	#	mean(df$MX/1e3),sd(df$MX/1e3))

  	#df2 <- rbind(data.frame(Site="AZ",value=df$AZ),
  	#	data.frame(Site="WTX",value=df$WTX),
  	#	data.frame(Site="ETX",value=df$ETX),
  	#	data.frame(Site="MX",value=df$MX))

  	#test <- pairwise.wilcox.test(df2$value, df2$Site,
    #                 p.adjust.method = "fdr")
  	#pvals[count,] <- c(paste("chr",chr,sep=""),as.numeric(na.omit(test$p.value[1:9])))



}












