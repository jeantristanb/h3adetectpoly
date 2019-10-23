#!/usr/bin/env Rscript
library("optparse")
#source('FctExtractAlleleV2.r')
Mode<-function(x){
tbx<-table(x)
return(as.numeric(names(tbx)[which.max(tbx)][1]))
}
GetChisq<-function(xI){
MinX<-max(min(xI),0)
x<-xI-min(xI)+1
ValTest=0:(max(x)+1)
MeanX<-mean(x)/length(x)
if(MeanX>1)return(GetChisq(xI[xI!=min(xI)]))
ResBinom<-data.frame(Val=ValTest,Bin=dbinom(ValTest,size=length(x), p=MeanX))
TabValI<-merge(ResBinom,as.data.frame(table(x)),all=T, by=1)
TabValI[is.na(TabValI$Freq),'Freq']<-0
TabValI[,1]<-TabValI[,1]+MinX-1
return(TabValI)
}



#x<-Data2$lenSeq[Data2$lenSeq>Min]
GetAllele<-function(x){
if(length(x)==0)return(list(resume=data.frame(A1=NA,A2=NA, NbA1=0, NbA2=0, Type=0), summres=NA, Chisq1Law=NA, Chisq2Law=NA))
if(length(unique(x))==1)return(list(resume=data.frame(A1=unique(x),A2=unique(x), NbA1=length(x)/2, NbA2=length(x)/2, Type=1),summres=NA ,Chisq1Law=NA, Chisq2Law=NA))
MeanX<-mean(range(x))
A1Dist=x[x<=MeanX]
A2Dist=x[x>MeanX]
if(length(A2Dist)<=2)return(GetAllele(A1Dist))
if(length(A1Dist)<=2)return(GetAllele(A2Dist))
library(MASS)
library(vcd)## loading vcd package 
gf<-goodfit(x,type= "poisson",method= "MinChisq")
gf1<-goodfit(A1Dist,type= "poisson",method= "MinChisq")
gf2<-goodfit(A2Dist,type= "poisson",method= "MinChisq")


gfmerge<-merge(merge(cbind(Observed=gf$observed,Count=gf$count,Nall=gf$fitted),cbind(Count=gf1$count,N1=gf1$fitted),all=T), cbind(Count=gf2$count,N2=gf2$fitted),all=T) 
gfmerge$Nall[is.na(gfmerge$Nall)]<-0;gfmerge$N1[is.na(gfmerge$N1)]<-0;gfmerge$N2[is.na(gfmerge$N2)]<-0
gfmerge$Pall<-gfmerge$Nall/sum(gfmerge$Nall)
gfmerge$Pall2law<-(gfmerge$N1+gfmerge$N2)/sum((gfmerge$N1+gfmerge$N2))
ResAll<-chisq.test(gfmerge$Observed,p=gfmerge$Pall,simulate.p.value=T)
ResAll2<-chisq.test(gfmerge$Observed,p=gfmerge$Pall2law,simulate.p.value=T)

if(ResAll$p.value>ResAll2$p.value)return(list(resume=data.frame(A1=Mode(x),A2=Mode(x), NbA1=length(x)/2, NbA2=length(x)/2, Type=2), summres=gfmerge,Chisq1Law=ResAll, Chisq2law=ResAll2))
else return(list(resume=data.frame(A1=Mode(A1Dist)[1],A2=Mode(A2Dist)[1], NbA1=length(A1Dist), NbA2=length(A2Dist), Type=3),  summres=gfmerge,Chisq1Law=ResAll, Chisq2law=ResAll2))
}





option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-e", "--header"), type="character", default="temp", 
              help="header name for analyse[default= %default]", metavar="character"),
  make_option(c("-l", "--lheader"), type="character", default="NbRepetI,NbRepetNewAl", 
              help="header name for analyse[default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
FileI=opt$file
Head=opt$header
Out=opt$out

#FileI="ATM0C.md.recal_bowtielocal.stat"
#Head='Tmp'
listheadernbrepet=strsplit(opt$lheader, split=",")[[1]]
listheadernbrepet<-listheadernbrepet[nchar(listheadernbrepet)>0]

Data<-read.table(FileI, header=T)
header1="NbRepetI"
svg(paste(Out,'.svg',sep=''),width = 7*2, height = 7)
par(mfrow=c(1,2))
Cmt=1
for(header1 in listheadernbrepet){
NbAllele<-GetAllele(Data[,header1])
ResumeNbAllele=NbAllele[['resume']]

tbval=table(Data[,header1])
if(nrow(Data)>0){
if(ResumeNbAllele$Type==0){
plot(as.integer(names(tbval)),unlist(tbval), type='h', ylim=range(0,max(tbval)), yaxt='n', lwd=4, xlab="Repetition Number", ylab='Observation number', main=paste('Inital :',Head), sub='No allele found')
axis(2,at=seq(0, max(tbval),2), label=seq(0, max(tbval),2))
}else if(ResumeNbAllele$Type==1){
Sub=paste('No test Allele done. Homoz. Allele : ',ResumeNbAllele$A1)
plot(as.integer(names(tbval)),unlist(tbval), type='h', ylim=range(0,max(tbval)), yaxt='n', lwd=4, xlab="Repetition Number", ylab='Observation number', main=paste('Inital :',Head), sub=Sub)
text(ResumeNbAllele$A1, ResumeNbAllele$NbA1+ ResumeNbAllele$NbA1, '*', cex=5, col='red')
axis(2,at=seq(0, max(tbval),2), label=seq(0, max(tbval),2))
}else if(ResumeNbAllele$Type==2){
Sub=paste('Test Allele done, Homoz Allele : ',ResumeNbAllele$A1)
DistAll<-NbAllele$summres
xran<-range(Data[,header1],  DistAll[,'Count'], na.rm=T)
yran<-range(0,tbval, DistAll[,c('N1','N2','Nall')],na.rm=T);yran[2]=yran[2]+1
plot(as.integer(names(tbval)),unlist(tbval), type='h', yaxt='n', lwd=4, xlab="Repetition Number", ylab='Observation number', main=paste('Inital :',Head), sub=Sub, xlim=xran, ylim=yran)
lines(DistAll[,1],DistAll[,'Nall'], lwd=2)
lines(DistAll[,1],DistAll[,'N1'], lwd=1, lty=2, col='red')
lines(DistAll[,1],DistAll[,'N2'], lwd=1, lty=3, col='blue')
axis(2,at=seq(0, max(tbval),2), label=seq(0, max(tbval),2))
text(ResumeNbAllele$A1, tbval[names(tbval)==ResumeNbAllele$A1], '*', cex=5, col='red')
legend('topright', legend=c('Model All', 'Model A1', 'Model A2'), lty=c(1,2,3), col=c('black', 'red', 'blue'), bty='n')
}else if(ResumeNbAllele$Type==3){
Sub=paste('Test Allele done, two models Allele : ',ResumeNbAllele$A1, ResumeNbAllele$A2)
DistAll<-NbAllele$summres
xran<-range(Data[,header1],  DistAll[,'Count'], na.rm=T)
yran<-range(0,tbval, DistAll[,c('N1','N2','Nall')],na.rm=T);yran[2]=yran[2]+1
plot(as.integer(names(tbval)),unlist(tbval), type='h', yaxt='n', lwd=4, xlab="Repetition Number", ylab='Observation number', main=paste('Inital :',Head), sub=Sub, xlim=xran, ylim=yran)
lines(DistAll[,1],DistAll[,'Nall'], lwd=1)
lines(DistAll[,1],DistAll[,'N1'], lwd=2, lty=2, col='red')
lines(DistAll[,1],DistAll[,'N2'], lwd=2, lty=3, col='blue')
axis(2,at=seq(0, max(tbval),2), label=seq(0, max(tbval),2))
text(ResumeNbAllele$A1, tbval[names(tbval)==ResumeNbAllele$A1], '*', cex=5, col='red')
text(ResumeNbAllele$A2, tbval[names(tbval)==ResumeNbAllele$A2], '*', cex=5, col='blue')
legend('topright', legend=c('Model All', 'Model A1', 'Model A2'), lty=c(1,2,3), col=c('black', 'red', 'blue'),bty='n')
}
}else{
plot(1:10,1:10,type='n')
text(c(5), 6, 'Not PE aligned', cex=2)
}
names(ResumeNbAllele)<-paste(names(ResumeNbAllele), header1, sep="_")
if(Cmt==1)AllResume<-ResumeNbAllele
else AllResume<-cbind(AllResume,ResumeNbAllele)
Cmt<-Cmt+1
}
dev.off()
AllResume<-cbind(Name=Head, AllResume)
write.table(AllResume, file=paste(Out,'_resume.stat',sep=''), col.names=T, row.names=F, quote=F)
