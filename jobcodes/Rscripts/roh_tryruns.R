library(sensitivity)

setwd("/scratch/bell/mathur20/MQU/ch3/revise/roh/")

##### ( For try runs 1) #####
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

# Get files and calculate Froh for each individual for each set of parameters

filenames <- list.files("tryruns1", pattern="*indiv", full.names=TRUE)

total <- 960796.788
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

# run sensitivity analysis for each individual for each set of parameters

popmeans <- cbind(rohs,params[1:nrow(rohs),])
write.table(popmeans,"popmeans.tryruns.txt",quote=F,row.names=F)

lm_az <- lm(m_az ~ snp + den + gap + het + wsnp + whet + wmiss + wthresh, data=popmeans) 
lm_mx <- lm(m_mx ~ snp + den + gap + het + wsnp + whet + wmiss + wthresh, data=popmeans) 
lm_wtx <- lm(m_wtx ~ snp + den + gap + het + wsnp + whet + wmiss + wthresh, data=popmeans) 
lm_etx <- lm(m_etx ~ snp + den + gap + het + wsnp + whet + wmiss + wthresh, data=popmeans) 

summary(lm_sum)

X <- indiv1[,c(2:9)]
y <- indiv1[,1]
x <- src(X,y, nboot = 100))

##################

###### FOR TRY RUNS 2 ########

# Get Parameters

columns <- c("snp","het","wsnp","whet","wmiss")
params <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(params) = columns

for (snp in c(50,75,100))
{
	for (het in c(1,2,3))
	{
		for (wsnp in c(10,20,30,40,50))
		{
			for (whet in c(1,2,3))
			{
				for (wmiss in c(1,2,5))
				{
					p <- c(snp,het,wsnp,whet,wmiss)
					params <- rbind(params,p)
				}
			}
		}
	}
}

colnames(params) <- columns

# Get files and calculate Froh for each individual for each set of parameters

filenames <- list.files("tryruns2", pattern="*indiv", full.names=TRUE)

total <- 960796.788
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
write.table(popmeans,"popmeans.tryruns2.txt",quote=F,row.names=F)

##################

###### FOR TRY RUNS 3 ########

# Get Parameters

columns <- c("snp","het","wmiss")
params <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(params) = columns

for (snp in c(25,50,60,70,75,80))
{
	for (het in c(0,1,2,3))
	{
			for (wmiss in c(1,2,3))
			{
				p <- c(snp,het,wmiss)
				params <- rbind(params,p)
			}
	}
}

colnames(params) <- columns

# Get files and calculate Froh for each individual for each set of parameters

filenames <- list.files("tryruns3", pattern="*indiv", full.names=TRUE)

total <- 960796.788
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
write.table(popmeans,"popmeans.tryruns3.txt",quote=F,row.names=F)

###### FOR TRY RUNS 4 ########

# Get Parameters

columns <- c("snp")
params <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(params) = columns

for (snp in c(70,80,90,100,110,120,130,140,150))
{
	p <- c(snp)
	params <- rbind(params,p)
}

colnames(params) <- columns

# Get files and calculate Froh for each individual for each set of parameters

filenames <- list.files("tryruns4", pattern="*indiv", full.names=TRUE)

total <- 960796.788
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
write.table(popmeans,"popmeans.tryruns4.txt",quote=F,row.names=F)