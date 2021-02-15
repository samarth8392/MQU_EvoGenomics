library("optparse")

option_list = list(
  make_option(c("-f","--file"), type="character",default=NULL,help="Plink hom file",metavar="character"),
  make_option(c("-o","--out"), type="character",default=NULL,help="Output file",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

map <- read.table(opt$file, header=T)

pos=map[,1]

rate=map[,2]*1e9

N=length(pos)

dm=c(0,1e-6*rate[-N]*(pos[-1]-pos[-N]))

m=cumsum(dm)

#m=m+map[1,3]

output <- cbind(pos,rate,m)


#colnames(output) <- c("Position(bp)","Rate(cM/Mb)","Map(cM)")
write.table(file=opt$out, output, quote=F)