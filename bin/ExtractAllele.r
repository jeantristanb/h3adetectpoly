source('FctExtractAllele.r')
Min<-10
Type='bwa'
Type2='bwalocal'
listFile=dir(paste('stat_test/align/',Type,'/',sep=''), pattern="*.out", full.names=T)

pdf('CmpVariousAlign.pdf')
CmtAll=1
for(File1 in listFile){
#ANG0J.md.recal_alignsim_detailseq.out
#File1="stat_test/align/bwa/ANG0J.md.recal_alignsim_detailseq.out"
par(mfrow=c(2,2))
File2<-paste('stat_test/align/',Type2,'/',basename(File1),sep='')
Data1<-read.table(File1, header=T)
Data2<-read.table(File2, header=T)
NbAlleleBWA<-GetAllele(Data1$NbRepRef)
NbAlleleBWALoc<-GetAllele(Data2$lenSeq[Data2$lenSeq>Min])
NbAlleleBWALoc<-GetAllele(Data2$lenSeq[Data2$lenSeq>Min])
hist(Data1$NbRepRef, main=paste(Type, basename(File1)), xlab="nb repet", 50, sub=paste(unlist(NbAlleleBWA),collapse=","))
hist(Data2$lenSeq[Data2$lenSeq>Min], main=Type2, xlab="nb repet",50,  sub=paste(unlist(NbAlleleBWALoc),collapse=","))
Intersect=merge(Data1,Data2, by=1, all=T)
if(length(unique(Data1$NbRepetNewAl))==1){
barplot(table(Data1$NbRepetNewAl), main=paste('real :',Type), xlab="nb repet")
}else{
hist(Data1$NbRepetNewAl, main=paste('real :',Type), xlab="nb repet",50)
}
if(length(unique(Data2$NbRepetNewAl))==1){
barplot(table(Data2$NbRepetNewAl), main=paste('real :',Type2), xlab="nb repet")
}else{
hist(Data2$NbRepetNewAl, main=paste('real :',Type2), xlab="nb repet", 50)
}
Intersect$File<-basename(File1)
if(CmtAll==1)AllRes<-Intersect else AllRes<-rbind(AllRes,Intersect)
CmtAll<-CmtAll+1
}
dev.off()
