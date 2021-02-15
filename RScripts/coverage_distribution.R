# Create frequency histograms for each individual sample based on samtools depth output",

setwd("/scratch/snyder/m/mathur20/MQU/ch3/align/MQU/final/stats/")

samples <- c("E8946",
"E8947",
"E8948",
"E8949",
"E8954",
"E9030",
"E9031",
"E9032",
"E9033",
"E9034",
"E9035",
"E9036",
"E9037",
"E9038",
"E9039",
"E9040",
"E9041",
"E9042",
"E9043",
"E9044",
"E9045",
"E9046",
"E9047",
"E9048",
"E9049",
"E9050",
"E9051",
"E9052",
"E9053",
"E9054",
"E9055",
"E9056",
"E9057",
"E9058",
"E9059")

for (i in samples)
{
	f_name <- paste(i, "_mqu_final_depth.txt", sep = "")
	f <- read.table(f_name, header= F)
	#name <- paste(i, sep="")
	#assign(name,f)
	setwd ("coverage_plots/")
	plot_name <- paste(i,"_coverage.pdf",sep="")
	pdf(plot_name, width=9, height=5)
	par(mar=c(5,5,5,5))
	#a <- get(name)
	hist(f$V3, breaks=100,
	main=paste(i,sep=""), 
	ylab="Bases",
    xlab="Depth of coverage", cex.lab=1)
    mtext(paste("Mean=",mean(f$V3), "SD=",signif(sd(f$V3),digits=5), sep=" "), side=3)
	dev.off()
	setwd("../")
	rm(f) 
}