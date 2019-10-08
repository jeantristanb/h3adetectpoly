#!/usr/bin/R
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
File="stat_test/out.depth"
Begin=45403049
Fin=45403083
Out="CmpDepthIndividuals.pdf"
}else{
File=args[1]
Begin=as.integer(args[2])
Fin=as.integer(args[3])
Out=args[4]
}


pdf(out)
data<-read.table(File, header=T)
for(Cmt in 3:ncol(data)){
plot(data[,2], data[,Cmt], xlab='Positions along genomes', ylab='depth', main=basenames(names(data)[Cmt]), cex.main=0.5)
lines(c(Begin, Fin),c(min(data[,Cmt]), min(data[,Cmt])), lwd=4, col='red')
}
dev.off()

