###########################################################################
###                          Samarth Mathur, PhD                     	  ###
###                        The Ohio State University                 	  ###
###                                                                     ###
###########################################################################
###########################################################################
###                   mutation_age.R           		                     ###
###########################################################################
#### PREREQUISITES #####
# load packages
library(ggplot2)
library(ggsignif)
library(RcppCNPy)
library(tidyr)
library(reshape2)
library(ggpubr)
library(cowplot)
library(plotrix)
library(wesanderson)
library(plyr)
library(dplyr)
library(DataCombine)
library(DescTools)
library(rcartocolor)
library(car)
library(detectRUNS)
library(tidyverse)
library(randomcoloR)
library(IdeoViz)

# Load colors
#display_carto_all(colorblind_friendly = TRUE)

col_az <- carto_pal(12,"Vivid")[11]
col_wtx <- carto_pal(7,"Purp")[5]
col_etx <-carto_pal(7,"Peach")[5]
col_mx <- "#3598FF"
colss <-c(rep(col_az,28),rep(col_mx,4),rep(col_wtx,31),rep(col_etx,3))
cols_l <- adjustcolor(cols, alpha.f = 0.5)
col_azwtx <- "#90CCC2"
cols <- c(col_az,col_wtx,col_etx, col_mx,col_azwtx)
#### PREREQUISITES: END #####

setwd("~/Documents/Thesis_Research/Final Results/Ch3/age/")

del_az <- read.csv("az.del.age.withFreq", header = T)
del_wtx <- read.csv("wtx.del.age.withFreq", header = T)

tol_az <- read.csv("az.weakdel.age.withFreq", header = T)
tol_wtx <- read.csv("wtx.weakdel.age.withFreq", header = T)

del_az <- del_az[-which(del_az$PostMode > 1.5*1e5),] # only using mutations in the last 150,000 years
del_wtx <- del_wtx[-which(del_wtx$PostMode > 1.5*1e5),]

tol_az <- tol_az[-which(tol_az$PostMode >  1.5*1e5),]
tol_wtx <- tol_wtx[-which(tol_wtx$PostMode >  1.5*1e5),]


del_az <- del_az[which(del_az$Clock == "J"),] # Using only joint clock # N= 12,224
del_wtx <- del_wtx[which(del_wtx$Clock == "J"),] # N = 11,406

tol_az <- tol_az[which(tol_az$Clock == "J"),] # N = 22,452
tol_wtx <- tol_wtx[which(tol_wtx$Clock == "J"),] # N = 22,318

del_df <- rbind(data.frame(Site="AZ", value = del_az$PostMode),
                data.frame(Site="WTX", value = del_wtx$PostMode))

tol_df <- rbind(data.frame(Site="AZ", value = tol_az$PostMode),
                data.frame(Site="WTX", value = tol_wtx$PostMode))


age_df <- rbind(data.frame(Type="Deleterious",Site=del_df$Site,Age=del_df$value),
                data.frame(Type="Tolerated",Site=tol_df$Site,Age=tol_df$value))

ggplot(age_df,aes(x=Type,y=Age,fill=Site))+
  geom_boxplot()

# Compare mean age of a mutation

mean(del_az$PostMode) # 40194.22
mean(del_wtx$PostMode) # 41947.83

wilcox.test(del_az$PostMode,del_wtx$PostMode) # W = 68398677, p-value = 0.0121

mean(tol_az$PostMode) # 40649.9
mean(tol_wtx$PostMode) # 43255.52

wilcox.test(tol_az$PostMode,tol_wtx$PostMode) # W = 261795028, p-value = 3.277e-14

del_azVar <- del_az[,c(1,2,10,11)]
del_wtxVar <- del_wtx[,c(1,2,10,11)]
for (i in 1:10)
{
  del_azVar <- cbind(del_azVar,del_az$PostMode*rnorm(12224,1.5,0.5))
  del_wtxVar <- cbind(del_wtxVar,del_wtx$PostMode*rnorm(11406,1.5,0.5))
}

snps <- NULL
for (i in 5:14)
{
  a <- length(which(del_azVar[,i] < 1e5))
  t <- length(which(del_wtxVar[,i] < 1e5))
  d <- data.frame(Run=i-4,AZ=a,WTX=t)
  snps <- rbind(snps,d)
  
}

pre <- NULL
for (i in 5:14)
{
  a <- length(which(del_azVar[,i] < 25000))
  t <- length(which(del_wtxVar[,i] < 25000))
  d <- data.frame(Run=i-4,AZ=a,WTX=t)
  pre <- rbind(pre,d)
  
}


# Plot no. of mutations by bin age
bin_del_az <- NULL
bin_del_wtx <- NULL
for (i in 5:14)
{
  a <- c(length((del_azVar$Position[which(del_azVar[,i]<100000 & del_azVar[,i]>=50000)])),
                  length((del_azVar$Position[which(del_azVar[,i]<50000 & del_azVar[,i]>=25000)])),
                  length((del_azVar$Position[which(del_azVar[,i]<25000 & del_azVar[,i]>=15000)])),
                  length((del_azVar$Position[which(del_azVar[,i]<15000 & del_azVar[,i]>=10000)])),
                  length((del_azVar$Position[which(del_azVar[,i]<10000 & del_azVar[,i]>=5000)])),
                  length((del_azVar$Position[which(del_azVar[,i]<5000 & del_azVar[,i]>=0)])))
  
  b <- c(length((del_wtxVar$Position[which(del_wtxVar[,i]<100000 & del_wtxVar[,i]>=50000)])),
         length((del_wtxVar$Position[which(del_wtxVar[,i]<50000 & del_wtxVar[,i]>=25000)])),
         length((del_wtxVar$Position[which(del_wtxVar[,i]<25000 & del_wtxVar[,i]>=15000)])),
         length((del_wtxVar$Position[which(del_wtxVar[,i]<15000 & del_wtxVar[,i]>=10000)])),
         length((del_wtxVar$Position[which(del_wtxVar[,i]<10000 & del_wtxVar[,i]>=5000)])),
         length((del_wtxVar$Position[which(del_wtxVar[,i]<5000 & del_wtxVar[,i]>=0)])))
  bin_del_az <- c(bin_del_az,a)
  bin_del_wtx <- c(bin_del_wtx,b)
  
}

time <- rep(c("50-100","25-50","15-25","10-15","5-10","0-5"),5)

snpsaz <- data.frame(Site=rep("AZ",30),Time=time,SNPs=bin_del_az)
snpswtx <- data.frame(Site=rep("WTX",30),Time=time,SNPs=bin_del_wtx)
age <- rbind(snpswtx,snpsaz)

pal1 <- c(wes_palette("GrandBudapest2"),wes_palette("GrandBudapest1"))
pal2 <- rev(pal1[1:6])

age$Time <- factor(age$Time, levels = c("50-100","25-50","15-25","10-15","5-10","0-5"))
age$Time <- factor(age$Time, levels = c("0-5","5-10","10-15","15-25","25-50","50-100"))

ageMean <- ddply(age, c("Site","Time"), summarise, MeanSNPs=round(mean(SNPs)),sdSNPs=round(sd(SNPs)))

ggplot(ageMean, aes(x=MeanSNPs, y=Site,fill=Time)) + 
  geom_bar(stat="identity", color = "black")+
  #geom_errorbar(aes(xmin=MeanSNPs-sdSNPs, xmax=MeanSNPs+sdSNPs), width=.2,position="identity")+
  theme_classic(base_size = 25)+
  scale_x_continuous(breaks=scales::pretty_breaks(n = 10))+
  xlab("No. of deleterious SNPs")+
  labs(fill = "Time (kya)")+
  theme(axis.title.y = element_blank())+
  scale_fill_manual(values=pal2)+
  theme(legend.title=element_blank())
        #legend.position = "none")
my_comparisons <- list( c("AZ", "WTX"))
#par(mar=c(5,10,2,10))

# PLot frequency by bottleneck age (25kya)

pre_bot_az <- NULL
pre_bot_wtx <- NULL
post_bot_az <- NULL
post_bot_wtx <- NULL
for (i in 5:14)
{
  a <- c(length((del_azVar$Position[which(del_azVar[,i]>= 25000 & del_azVar[,i]< 100000)])))
  b <- c(length((del_wtxVar$Position[which(del_wtxVar[,i]>= 25000 & del_wtxVar[,i]< 100000)])))
  
  c <- c(length((del_azVar$Position[which(del_azVar[,i] <= 15000)])))
  d <- c(length((del_wtxVar$Position[which(del_wtxVar[,i] <= 15000)])))
  
  pre_bot_az <- c(pre_bot_az,a)
  pre_bot_wtx <- c(pre_bot_wtx,b)
  post_bot_az <- c(post_bot_az,c)
  post_bot_wtx <- c(post_bot_wtx,d)
         
  
}

df_snps <- rbind(data.frame(Site="AZ",Type="Pre-Bottleneck",snps=pre_bot_az),
                 data.frame(Site="WTX",Type="Pre-Bottleneck",snps=pre_bot_wtx),
                 data.frame(Site="AZ",Type="Post-Bottleneck",snps=post_bot_az),
                 data.frame(Site="WTX",Type="Post-Bottleneck",snps=post_bot_wtx))
                 
my_comparisons <- list( c("AZ", "WTX"))

# Plot no. of snps
df_snps$Site <- factor(df_snps$Site , levels=c("AZ", "WTX"))
df_snps$Type <- factor(df_snps$Type , levels=c("Pre-Bottleneck", "Post-Bottleneck"))

snpsMean <- ddply(df_snps, c("Site","Type"),summarise, MeanSNPs=round(mean(snps)),sdSNPs=round(sd(snps)))

ggplot(snpsMean, aes(fill=Site, y=MeanSNPs, x=Type)) + 
  geom_bar(position="dodge", stat="identity")+
  geom_errorbar(aes(ymin=MeanSNPs-sdSNPs, ymax=MeanSNPs+sdSNPs), width=.2, position=position_dodge(.9))+
  theme_classic(base_size = 15)+ylim(0,4500)+
  scale_fill_manual(values=cols[1:2])+
  labs(x="Time", y ="No. of deleterious SNPs")+
  theme(axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=15))

wilcox.test(pre_bot_az,pre_bot_wtx) # W = 100, p-value = 0.0001817
wilcox.test(post_bot_az,post_bot_wtx) # W = 100, p-value = 0.0001786

# Models for age
# All
df_az <- df_freq[which(df_freq$Site == "AZ"),]
df_wtx <- df_freq[which(df_freq$Site == "WTX"),]

df_freq <- rbind(data.frame(Site="AZ",Type="Pre-Bottleneck",snps=length(pre_bot_az$PostMode),freq=pre_bot_az$freq_B,age=pre_bot_az$PostMode),
                 data.frame(Site="AZ",Type="Post-Bottleneck",snps=length(post_bot_az$PostMode),freq=post_bot_az$freq_B,age=post_bot_az$PostMode),
                 data.frame(Site="WTX",Type="Pre-Bottleneck",snps=length(pre_bot_wtx$PostMode),freq=pre_bot_wtx$freq_B,age=pre_bot_wtx$PostMode),
                 data.frame(Site="WTX",Type="Post-Bottleneck",snps=length(post_bot_wtx$PostMode),freq=post_bot_wtx$freq_B,age=post_bot_wtx$PostMode))

mean(pre_bot_az$freq_B) # 0.1968249
mean(pre_bot_wtx$freq_B) # 0.2130443

mean(post_bot_az$freq_B) # 0.04703992
mean(post_bot_wtx$freq_B) # 0.05428657

ggplot(df_freq,aes(fill=Site, y=freq, x=age))+ylim(0,1)+
  geom_point(color=adjustcolor(col_az,alpha.f = 0.2), size=0.5, shape=16)+ 
  geom_point(color=adjustcolor(col_wtx,alpha.f = 0.2),  size=0.5, shape=16)+ 
  geom_smooth(method=loess, span=0.80, aes(x=age, y=freq), data = df_az, color=col_az, fill="#69b3a2", se=F, size=2.5)+
  geom_smooth(method=loess, span=0.80, aes(x=age, y=freq), data = df_wtx, color=col_wtx, fill="#69b3a2", se=F, size=2.5)+
  geom_smooth(method=loess, span=0.80, aes(x=age, y=freq), data = df_az, color="black", fill="#69b3a2", se=F, size=0.5)+
  geom_smooth(method=loess, span=0.80, aes(x=age, y=freq), data = df_wtx, color="black", fill="#69b3a2", se=F, size=0.5)+
  theme_classic(base_size = 15)+
  labs(x="Age of Mutation (Ya)", y="Allele frequency (MAF)")

# Pre
pre_df <- df_freq[which(df_freq$Type == "Pre-Bottleneck"),]
pre_df_az <- pre_df[which(pre_df$Site == "AZ"),]
pre_df_wtx <- pre_df[which(pre_df$Site == "WTX"),]

ggplot(pre_df,aes(fill=Site, y=freq, x=age))+ylim(0,0.1)+
  geom_point(color=adjustcolor(col_az,alpha.f = 0.4), size=0.5, shape=16)+ 
  geom_point(color=adjustcolor(col_wtx,alpha.f = 0.4),  size=0.5, shape=16)+ 
  geom_smooth(method=loess, span=0.80, aes(x=age, y=freq), data = pre_df_az, color=col_az, fill="#69b3a2", se=F, size=2.5)+
  geom_smooth(method=loess, span=0.80, aes(x=age, y=freq), data = pre_df_wtx, color=col_wtx, fill="#69b3a2", se=F, size=2.5)+
  geom_smooth(method=loess, span=0.80, aes(x=age, y=freq), data = pre_df_az, color="black", fill=col_az, se=T, size=0.5)+
  geom_smooth(method=loess, span=0.80, aes(x=age, y=freq), data = pre_df_wtx, color="black", fill=col_wtx, se=T, size=0.5)+
  theme_classic(base_size = 15)+
  labs(x="Age of Mutation (Ya)", y="Allele frequency (MAF)")

summary(loess(pre_df_az$freq ~pre_df_az$age, span=0.5))

# Post
post_df <- df_freq[which(df_freq$Type == "Post-Bottleneck"),]
post_df_az <- post_df[which(post_df$Site == "AZ"),]
post_df_wtx <- post_df[which(post_df$Site == "WTX"),]

ggplot(post_df,aes(fill=Site, y=freq, x=age))+ylim(0,0.1)+
  geom_point(color=adjustcolor(col_az,alpha.f = 0.4), size=0.5, shape=16)+ 
  geom_point(color=adjustcolor(col_wtx,alpha.f = 0.4),  size=0.5, shape=16)+ 
  geom_smooth(method=loess, span=0.80, aes(x=age, y=freq), data = post_df_az, color=col_az, fill="#69b3a2", se=F, size=2.5)+
  geom_smooth(method=loess, span=0.80, aes(x=age, y=freq), data = post_df_wtx, color=col_wtx, fill="#69b3a2", se=F, size=2.5)+
  geom_smooth(method=loess, span=0.80, aes(x=age, y=freq), data = post_df_az, color="black", fill=col_az, se=F, size=0.5)+
  geom_smooth(method=loess, span=0.80, aes(x=age, y=freq), data = post_df_wtx, color="black", fill=col_wtx, se=F, size=0.5)+
  theme_classic(base_size = 15)+
  labs(x="Age of Mutation (Kya)", y="Allele frequency (MAF)")



# Shared vs Private

az_del_post <- del_az[which(del_az$PostMode <= 25000),] # 6114
wtx_del_post <- del_wtx[which(del_wtx$PostMode <= 25000),] # 5308

del_post_shared <- inner_join(az_del_post,wtx_del_post, by=c("Chromosome","Position")) # 2328

del_post_priv_az <- anti_join(az_del_post,del_post_shared,by=c("Chromosome","Position")) # 3786
del_post_priv_wtx <- anti_join(wtx_del_post,del_post_shared,by=c("Chromosome","Position")) # 2980

df <- rbind(data.frame(Site="AZ",Type="Shared",freq=del_post_shared$freq_B.x,age=del_post_shared$PostMode.x),
            data.frame(Site="AZ",Type="Private",freq=del_post_priv_az$freq_B,age=del_post_priv_az$PostMode),
            data.frame(Site="WTX",Type="Shared",freq=del_post_shared$freq_B.y,age=del_post_shared$PostMode.y),
            data.frame(Site="WTX",Type="Private",freq=del_post_priv_wtx$freq_B,age=del_post_priv_wtx$PostMode))

# Shared
df2 <- df[which(df$Type=="Shared"),]
ggplot(df2,aes(fill=Site, y=freq, x=age))+ylim(0.01,0.1)+
  geom_point(color=adjustcolor(col_az,alpha.f = 0.4), size=0.5, shape=16)+ 
  geom_point(color=adjustcolor(col_wtx,alpha.f = 0.4),  size=0.5, shape=16)+ 
  geom_smooth(method=loess, span=0.80, aes(x=PostMode.x, y=freq_B.x), data = del_post_shared, color=col_az, fill="#69b3a2", se=F, size=2.5)+
  geom_smooth(method=loess, span=0.80, aes(x=PostMode.y, y=freq_B.y), data = del_post_shared, color=col_wtx, fill="#69b3a2", se=F, size=2.5)+
  geom_smooth(method=loess, span=0.80, aes(x=PostMode.x, y=freq_B.x), data = del_post_shared, color="black", fill=col_az, se=T, size=0.5)+
  geom_smooth(method=loess, span=0.80, aes(x=PostMode.y, y=freq_B.y), data = del_post_shared, color="black", fill=col_wtx, se=T, size=0.5)+
  labs(x="Years before present", y="Derived Allele Frequency")+
  theme_classic(base_size = 18)+scale_x_reverse()+
  theme(legend.title=element_blank(),
        legend.position = "none")

# Private
df3 <- df[which(df$Type=="Private"),]

ggplot(df3,aes(fill=Site, y=freq, x=age))+ylim(0.01,0.1)+
  geom_point(color=adjustcolor(col_az,alpha.f = 0.4), size=0.5, shape=16)+ 
  geom_point(color=adjustcolor(col_wtx,alpha.f = 0.4),  size=0.5, shape=16)+ 
  geom_smooth(method=loess, span=0.80, aes(x=PostMode, y=freq_B), data = del_post_priv_az, color=col_az, fill="#69b3a2", se=F, size=2.5)+
  geom_smooth(method=loess, span=0.80, aes(x=PostMode, y=freq_B), data = del_post_priv_wtx, color=col_wtx, fill="#69b3a2", se=F, size=2.5)+
  geom_smooth(method=loess, span=0.80, aes(x=PostMode, y=freq_B), data = del_post_priv_az, color="black", fill=col_az, se=T, size=0.5)+
  geom_smooth(method=loess, span=0.80, aes(x=PostMode, y=freq_B), data = del_post_priv_wtx, color="black", fill=col_wtx, se=T, size=0.5)+
  labs(x="Years before present", y="Derived Allele Frequency")+
  theme_classic(base_size = 18)+scale_x_reverse()+
  theme(legend.title=element_blank(),
        legend.position = "none")


## Plot shared/private density

shared <- read.csv("shared.del.age.withFreq",header = T)
priv.az <- read.csv("priv.az.del.age.withFreq",header = T)
priv.wtx <- read.csv("priv.wtx.del.age.withFreq",header = T)

shared <- shared[-which(shared$PostMode > 1.5*1e5),] # only using mutations in the last 100,000 years
priv.az <- priv.az[-which(priv.az$PostMode > 1.5*1e5),] 
priv.wtx <- priv.wtx[-which(priv.wtx$PostMode > 1.5*1e5),] 

shared <- shared[which(shared$Clock == "J"),] # Using only joint clock # N= 7147
priv.az <- priv.az[which(priv.az$Clock == "J"),] # N = 4212
priv.wtx <- priv.wtx[which(priv.wtx$Clock == "J"),] # N = 3399

age_df <- rbind(data.frame(Site="Shared", Age = shared$PostMode),
                data.frame(Site="Private AZ", Age = priv.az$PostMode),
                data.frame(Site="Private WTX", Age = priv.wtx$PostMode))

age_df$Site <- factor(age_df$Site, levels=c("Shared","Private AZ","Private WTX"))

# PLot 10 runs
a <- age_df
for (i in 1:9)
{
  a <- cbind(a,age_df$Age*rnorm(nrow(age_df),1,0.4))
}

snps <- NULL
for (i in 2:11)
{
  sh <- length(which(a[,i] < 1e5))
  paz <- length(which(del_wtxVar[,i] < 1e5))
  ptx <- data.frame(Run=i-4,AZ=a,WTX=t)
  snps <- rbind(snps,d)
  
}

age_df50k <- age_df[which(age_df$Age <=50000),]
age_df150k <- age_df[which(age_df$value <=150000),]

shared_age <- age_df50k[which(age_df50k$Site=="Shared"),]
privaz_age <- age_df50k[which(age_df50k$Site=="Private AZ"),]
privwtx_age <- age_df50k[which(age_df50k$Site=="Private WTX"),]

for (i in 1:9)
{
  shared_age <- cbind(shared_age,shared_age$Age*rnorm(nrow(shared_age),1,0.4))
  privaz_age <- cbind(privaz_age,privaz_age$Age*rnorm(nrow(privaz_age),1,0.4))
  privwtx_age <- cbind(privwtx_age,privwtx_age$Age*rnorm(nrow(privwtx_age),1,0))
}

p0 <- ggplot(age_df50k) +
  xlim(0,50000) + xlab(expression("Age of mutation (Ya)"))+
  ylab("No. of deleterious mutations")+ ylim(0,200)+
  theme_classic(base_size = 10)+
  scale_fill_manual(values=cols[c(1,2,5)])+
  theme(legend.title=element_blank(),
        legend.position = "none")

for (i in 2:11)
{
  s <- shared_age[,c(1,i)]
  paz <- privaz_age[,c(1,i)]
  ptx <- privwtx_age[,c(1,i)]

  colnames(s)[2] <- "value"
  colnames(paz)[2] <- "value"
  colnames(ptx)[2] <- "value"
  
  name <- paste("p",i,sep="")
  pi <- p0 + geom_histogram(alpha=0.9, bins = 100, aes(value), data= s, fill=cols[5], col="black")+
    geom_histogram(alpha=0.9, bins = 100, aes(value), data= paz, fill=col_wtx,col="black")+
    geom_histogram(alpha=0.9, bins = 100, aes(value), data= ptx, fill=col_az, col="black")
  assign(name,pi)
}

grid.arrange(p2,p3,p4,p5,p6,p7,p8,p9,p10,p11, nrow = 5)

p0 <- ggplot(age_df150k) +
  geom_histogram(alpha=0.5, bins = 100, aes(value), data= shared_age, fill=cols[5], col="black")+
  geom_histogram(alpha=0.5, bins = 100, aes(value), data= privwtx_age, fill=col_wtx,col="black")+
  geom_histogram(alpha=0.5, bins = 100, aes(value), data= privaz_age, fill=col_az, col="black")+
  xlim(0,150000) + xlab(expression("Age of mutation (Ya)"))+
  ylab("No. of deleterious mutations")+ ylim(0,450)+
  theme_classic(base_size = 20)+
  #scale_fill_manual(values=cols[c(1,2,5)])+
  theme(axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title=element_blank(),
        legend.position = "none")

p0

age_df <- cbind(age_df,age_df$value*rnorm(14758,1.5,0.5),age_df$value*rnorm(14758,1.5,0.5))
colnames(age_df) <- c("Site","Mean","Min","Max")



p1 <- p0+  geom_histogram(alpha=0.7, bins = 100, aes(Mean), col="black") 

p2 <- p0 + geom_histogram(alpha=0.5, bins = 100, aes(Min))
p3 <- p0 + geom_histogram(alpha=0.9, bins = 100, aes(Max))

grid.arrange(p2,p1, p3, nrow = 1)

# Get genotypes for shared/private

setwd("~/Documents/Thesis_Research/Final Results/Ch3/revise/age/")

hetero <- read.table("Heterozygotes.txt",header = T)
homo <- read.table("AltHomos.txt",header = T)

priv.Total <- c(rep(dim(priv.az)[1],28),rep(dim(priv.wtx)[1],31))

homo_prop <- cbind(homo,homo$PrivDel/priv.Total,homo$SharedDel/dim(shared)[1])
colnames(homo_prop)[5:6] <- c("PrivProp","SharedProp")

my_comparisons <- list( c("AZ", "WTX"))
                        

df.shpriv <- rbind(data.frame(Site="AZ",Type="Shared",SNPs=homo_prop$SharedDel[which(homo_prop$Site == "AZ")]),
                   data.frame(Site="AZ",Type="Private",SNPs=homo_prop$PrivDel[which(homo_prop$Site == "AZ")]),
                   data.frame(Site="WTX",Type="Shared",SNPs=homo_prop$SharedDel[which(homo_prop$Site == "WTX")]),
                   data.frame(Site="WTX",Type="Private",SNPs=homo_prop$PrivDel[which(homo_prop$Site == "WTX")]))

ggplot(df.shpriv,aes(x=Type,y=SNPs,fill=Site))+
  geom_boxplot(alpha = 0.5,outlier.shape = NA, position=position_dodge())+
  geom_point(position=position_jitterdodge(),size=2,aes(color=Site),color="black",pch=21)+
  labs(x="Site", y = expression("No. of derived homozygous SNPs"))+
  theme_classic(base_size = 22)+
  scale_fill_manual(values=cols[1:2])+
  scale_color_manual(values=cols[1:2])+
  theme(axis.title.y=element_text(size=17),
    axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position = "none")


wilcox.test(df.shpriv$SNPs[which(df.shpriv$Type=="Shared" & df.shpriv$Site=="AZ")],
            df.shpriv$SNPs[which(df.shpriv$Type=="Shared" & df.shpriv$Site=="WTX")],conf.int = T) # W = 53.5, p-value = 7.997e-09

wilcox.test(df.shpriv$SNPs[which(df.shpriv$Type=="Private" & df.shpriv$Site=="AZ")],
            df.shpriv$SNPs[which(df.shpriv$Type=="Private" & df.shpriv$Site=="WTX")],conf.int = T) # W = 241.5, p-value = 0.003558
