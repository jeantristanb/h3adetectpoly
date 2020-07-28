#!/usr/bin/env python3
import sys
import argparse

def GetChroFasta(File, chro, fileseq=None):
 read=open(File) 
 Chaine=""
 Save=False
 for line in read :
   if line[0]==">" :  
      Chro=line.replace('\n','').replace('>','').split(' ')[0]#.replace('chr','').replace('Chr','')
      if Save==True and Chro!=chro :
        if fileseq!=None :
          ecrire=open(fileseq,'w')
          ecrire.write('>'+chro+'\n'+Chaine+'\n')
          ecrire.close()
        return Chaine
      if Chro==chro :
         Save=True
   elif Save :
        Chaine+=line.replace('\n','')
 return Chaine

def CheckPolyElement(sequence, repet, Check):
   Final=''
   try :
     nbpolyelement=int(len(sequence)/float(len(repet)))
   except :
      Out1='Probleme with '+repet +' in sequence :' + sequence +' (reference genome) repetetion number is not possible : '+str(len(sequence)/float(len(repet)))
      if Check :
          Final+=Out1+'\n'
      else :
        sys.exit(1) 
   if repet*nbpolyelement!=sequence :
      Out2='Problem with '+repet +' in sequence :' + sequence+' (reference genome) with '+str(len(sequence)/float(len(repet)))+' repetitions different of sequence generated '+ repet*nbpolyelement
      if Check :
        Final+=Out2
      else :
         print('Problem with '+repet +' in sequence :' + sequence+' (reference genome) with '+str(len(sequence)/float(len(repet)))+' repetitions different of sequence generated '+ repet*nbpolyelement)
         sys.exit(1) 
   #print('nb repetition observed : '+str(nbpolyelement))
   return (nbpolyelement,Final)
def CheckAroundPolyelem(SeqChro, repet,posbegin, posend):
   Chaine=''
   Before=SeqChro[(posbegin-len(repet)-1):(posbegin-1)]
   After=SeqChro[posend:(posend+len(repet))]
   if Before==repet :
     Chaine+="Warning : bases before element contain repetition"
   if After==repet :
     Chaine+="Warning : bases after element contain repetition"
   return Chaine
     
def parseArguments():
    parser = argparse.ArgumentParser(description='extract sequence of interrest in fasta file')
    parser.add_argument('--fasta',type=str,required=True, help="fasta file")
    parser.add_argument('--repet',type=str,required=False, help="fasta file")
    parser.add_argument('--chro',type=str,required=True,help="chro")
    parser.add_argument('--posbegin', type=int,help="pos begin for chro", required=True)
    parser.add_argument('--posend', type=int,help="pos end of sequence of chro", required=True)
    parser.add_argument('--around', type=int,help="how many position around", required=False, default=0)
    parser.add_argument('--around2', type=int,help="how many position around", required=False, default=-1)
    parser.add_argument('--minnbrepet', type=int,help="how many position around", required=False, default=0)
    parser.add_argument('--maxnbrepet', type=int,help="how many position around", required=False, default=100)
    parser.add_argument('--out', type=str,help="out header",default="outseq", required=False)
    parser.add_argument('--check',action="store_true")
    args = parser.parse_args()
    return args

## first
args=parseArguments()


OutHead=args.out
chro=args.chro
#chro=chro.replace('chr','').replace('Chr','')
posend=args.posend
posbegin=args.posbegin
around=args.around
around2=args.around2
if around2<0 :
  around2=around
repet=args.repet.upper()
LogFile=open(OutHead+'.log','w')

if args.check :
   Type='Check'
else :
   Type= ' generate sequences '
LogFile.write('#Report to '+Type+'\n')
LogFile.write('##Params:\n')
res=str(args).replace('Namespace(','').replace(')','').split(',')
LogFile.write('###'+'\n###'.join(res)+'\n')

   


if around<0 :
   sys.exit("around negative : "+str(around))
if posend < posbegin or posend <=0 or posbegin<=0 :
   sys.exit('error in posbegin '+ str(posbegin)+' or posend '+str(posend)+': posend < posbegin or posend/posbegin<0 ')
### keep sequence
SeqChro=GetChroFasta(args.fasta, chro)
## check sequence
LogFile.write('## Check sequences and parameters: ')
LenChro=len(SeqChro)
if LenChro==0 :
   Message='not found chro '+chro+' in fasta file '+args.fasta
   if args.check :
     LogFile.write(Message)
     sys.exit(1)
   else :
     LogFile.write(Message)
     sys.exit(1)

if LenChro<posend :
   Message='pos end '+ str(posend)+' is higher than len of chro '+chro+' :'+str(len(SeqChro))+' in fasta file '+args.fasta
   if args.check :
     LogFile.write(Message)
     sys.exit(1)
   else :
     LogFile.write(Message)
     sys.exit(1)
## keep polyelement
PolyElementSeq=SeqChro[(posbegin-1):(posend)]
## check around poly element
FinalAroundPoly=CheckAroundPolyelem(SeqChro, repet,posbegin, posend)

(NbRepetition,FinalInfo)=CheckPolyElement(PolyElementSeq, repet,sys.exit)
LogFile.write('###Nb times of repetition '+str(NbRepetition)+'\n')
LogFile.write('###Error Message in Check element poly :\n'+FinalInfo+'\n')
LogFile.write('###Warning Message around element poly :\n'+FinalAroundPoly+'\n')

if args.check :
   sys.exit(0)


####### write sequences

   
PolyElement=">"+str(chro)+'_'+str(posbegin)+'_'+str(posend)+'\n'+PolyElementSeq+'\n'
Write=open(OutHead+".repet.fasta", 'w')
Write.write(PolyElement)
Write.close()


## around 
posbegaround=max(posbegin-around,0)
posendaround=min(posend+around,LenChro)
newposbeg=around+(posbegaround-(posbegin-around))+1
newposend=newposbeg+len(PolyElementSeq)-1
SeqInterest=">"+'chr'+str(chro)+'_'+str(posbegaround)+'_'+str(posendaround)+' '+str(newposbeg)+'_'+str(newposend)+'\n'
SeqInterest+=SeqChro[(posbegaround-1):(posendaround)]+'\n'

Write=open(OutHead+".around.fasta", 'w')
Write.write(SeqInterest)
Write.close()

newposbeg=around+(posbegaround-(posbegin-around))+1
seqbefore=SeqChro[(posbegaround-1):(posbegin-1)]
seqafter=SeqChro[posend:posendaround]

Write=open(OutHead+".nopolyt.fasta", 'w')
SeqInterest=">"+'chr'+str(chro)+'_'+str(posbegaround)+'_'+str(posendaround)+'\n'
SeqInterest+=seqbefore+seqafter+'\n'
Write.write(SeqInterest)
Write.close()

Write=open(OutHead+".begin.fasta", 'w')
SeqInterest=">"+'chr'+str(chro)+'_'+str(posbegaround)+'_'+str(posendaround)+'\n'
SeqInterest+=seqbefore+'\n'
Write.write(SeqInterest)
Write.close()

Write=open(OutHead+".end.fasta", 'w')
SeqInterest=">"+'chr'+str(chro)+'_'+str(posbegaround)+'_'+str(posendaround)+'\n'
SeqInterest+=seqafter+'\n'
Write.write(SeqInterest)
Write.close()



Write=open(OutHead+".insertvar.fasta", 'w')
for CmtRep in range(args.minnbrepet,args.maxnbrepet+1):
   NewPolyElem=repet*CmtRep
   newposend=newposbeg+len(NewPolyElem)-1
   SeqSimul=">"'Nb'+str(CmtRep)+'_'+str(newposbeg)+'_'+str(newposend)+'\n'
   SeqSimul+= seqbefore+NewPolyElem+seqafter+'\n'
   Write.write(SeqSimul)
Write.close()
