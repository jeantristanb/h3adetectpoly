
Mode = function(x){
    ta = table(x)
    if(length(ta)==1)return(unique(x))
    tam = max(ta)
    if (all(ta == tam))
         mod = min(ta)
    else
         if(is.numeric(x))
    mod = as.numeric(names(ta)[ta == tam])
    else
         mod = names(ta)[ta == tam]
    return(mod)
}

Mode<-function(x){
tbx<-table(x)
return(as.numeric(names(tbx)[which.max(tbx)][1]))
}
#x<-Data2$lenSeq[Data2$lenSeq>Min]
GetAllele<-function(x){
if(length(x)==0)return(list(resume=data.frame(A1=NA,A2=NA, NbA1=0, NbA2=0, Type=0), DistAll=NA, DistA1=NA, DistA2=NA))
getchisq<-function(x){
MinX<-max(min(x),0)
x<-x-min(x)+1
ValTest=0:(max(x)+1)
MeanX<-mean(x)/length(x)
print(MeanX)
ResBinom<-data.frame(Val=ValTest,Bin=dbinom(ValTest,size=length(x), p=MeanX))
TabValI<-merge(ResBinom,as.data.frame(table(x)),all=T, by=1)
TabValI[is.na(TabValI$Freq),'Freq']<-0
#TabValI<-TabValI[TabValI$Freq>0,]
TabValI[,1]<-TabValI[,1]+MinX-1
return(TabValI)
}
if(length(unique(x))==1)return(list(resume=data.frame(A1=unique(x),A2=unique(x), NbA1=length(x)/2, NbA2=length(x)/2, Type=1), DistAll=NA, DistA1=NA, DistA2=NA))
MeanX<-mean(range(x))
A1Dist=x[x<=MeanX]
A2Dist=x[x>=MeanX]
if(length(A2Dist)<=2)return(GetAllele(A1Dist))
if(length(A1Dist)<=2)return(GetAllele(A2Dist))
ResAllI<-getchisq(x)
ResAll<-ResAllI[ResAllI$Freq>0, ]
ResA1<-getchisq(A1Dist)
ResA2<-getchisq(A2Dist)
ResModel2<-rbind(ResA1[ResA1$Freq>0,],ResA2[ResA1$Freq>0,])
ResModel2$Bin<-ResModel2$Bin/sum(ResModel2$Bin)
ResAll$Bin<-ResAll$Bin/sum(ResAll$Bin)
ResModel2$Freq[ResModel2[,1]==MeanX]<-ResModel2$Freq[ResModel2[,1]==MeanX]/2
ResModel2<-aggregate(.~Val,ResModel2, sum)
ChisqMod=fisher.test(rbind(ResModel2$Freq, round(ResModel2$Bin*sum(ResModel2$Freq))), simulate.p.value=T)
ChisqAll=fisher.test(rbind(ResAll$Freq, round(ResAll$Bin*sum(ResModel2$Freq))), simulate.p.value=T)

if(ChisqAll$p.value>ChisqMod$p.value)return(list(resume=data.frame(A1=Mode(x),A2=Mode(x), NbA1=length(x)/2, NbA2=length(x)/2, Type=2), DistAll=ResAllI, DistA1=ResA1, DistA2=ResA2))
else return(list(resume=data.frame(A1=Mode(A1Dist)[1],A2=Mode(A2Dist)[1], NbA1=length(A1Dist), NbA2=length(A2Dist), Type=3), DistAll=ResAllI, DistA1=ResA1, DistA2=ResA2))
}
