#!/usr/bin/env python3

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys
import argparse


def DoAlignment(SeqToCheck, Polyele, MinRep, MaxRep, LenAround):
   ListScore=[]
   SaveMinAl={}
   Polyele=Polyele.upper()
   for CmtPoly in range(MinRep,MaxRep) :
      SeqPolyele=CmtPoly*Polyele
      ##
      alignments = pairwise2.align.localms(SeqToCheck, SeqPolyele, 1, -2, -5, -0.5)
      ScorAl=[x[2] for x in alignments]
      ScorMin=min(ScorAl)
      SaveMinAl[CmtPoly]=alignments[ScorAl.index(ScorMin)]
      ListScore.append(ScorMin)
   KeyMaxScore=MinRep+ListScore.index(max(ListScore))
   if SaveMinAl[KeyMaxScore][4]>=(len(SeqToCheck)-LenAround) or SaveMinAl[KeyMaxScore][3]<=LenAround :
      return("NA", "NA", SaveMinAl)
   return(max(ListScore),MinRep+ListScore.index(max(ListScore)), SaveMinAl)




def intersect(a, b):
     return list(set(a) & set(b))

ListeFlag=set(['I', 'M', 'D', 'N', 'S', 'H'])
def GetInfoFlag(FlagSam, detail=False):
   Chaine=""
   Len=0
   listflag=[]
   listnbflag=[]
   for x in FlagSam :
      if x in ListeFlag :
          if x!='I' and x!='S' and x!='H':
              Len+=int(Chaine)
          listnbflag.append(int(Chaine))
          listflag.append(x)
          Chaine=""
      else :
          Chaine+=x
   if detail :
      return (Len,listflag,listnbflag)
   return Len



def ReadUnMapped(flag) :
   if flag & 4 :
      return True
   return False

def IsSecondPE(flag):
  if flag & 64 :
     return True
  return False


def parseArguments():
    parser = argparse.ArgumentParser(description='extract sequence of interrest in fasta file')
    parser.add_argument('--fasta',type=str,required=False, help="fasta file")
    parser.add_argument('--sam_begin',type=str,required=True, help="fasta file")
    parser.add_argument('--sam_end',type=str,required=True, help="fasta file")
    parser.add_argument('--out', type=str,help="out header",default="outseq", required=False)
    parser.add_argument('--oldpos_poly', type=int,help="out header", required=True)
    parser.add_argument('--minnbrepet', type=int,help="how many position around", required=False, default=0)
    parser.add_argument('--maxnbrepet', type=int,help="how many position around", required=False, default=100)
    parser.add_argument('--repet',type=str,required=False, help="fasta file")
    parser.add_argument('--nb_bp_threshold', type=int,help="out header", required=False, default=5)
    args = parser.parse_args()
    return args


args=parseArguments()
readsam1=open(args.sam_begin)
readsam2=open(args.sam_end)
around=args.nb_bp_threshold
oldpospoly=args.oldpos_poly

dicseqbegin={}
for line in readsam1 :
    if line[0]=='@':
      continue
    splline=line.split()
    FlagSam=int(splline[1])
    if ReadUnMapped(int(splline[1])) :
       continue
    splline=line.split('\t')
    BegInSeq=int(splline[3])
    Cigar=splline[5]
    (lenSeq,listflag,listnbflag)=GetInfoFlag(Cigar,True)
    EndInSeq=BegInSeq+lenSeq
    Balise=EndInSeq>=(oldpospoly-around) and listflag[len(listflag)-1] in ['S','H', 'I'] 
    PeType="Pe1"
    if IsSecondPE(FlagSam) :
       PeType="Pe2"
    #print('sam1',splline[0], PeType, FlagSam, splline)
    if Balise:
      dicseqbegin[splline[0]+PeType]=splline

dicseqend={}
for line in readsam2 :
    if line[0]=='@':
      continue
    splline=line.split()
    FlagSam=int(splline[1])
    if ReadUnMapped(int(splline[1])) :
       continue
    splline=line.split('\t')
    BegInSeq=int(splline[3])
    Cigar=splline[5]
    (lenSeq,listflag,listnbflag)=GetInfoFlag(Cigar,True)
    Balise=BegInSeq<=around and listflag[0] in ['S','H','I'] #and (listflag.count('S')+listflag.count('H'))==1 
    PeType="Pe1"
    if IsSecondPE(FlagSam) :
       PeType="Pe2"
    #print('sam2',splline[0], PeType, FlagSam, splline)
    if Balise :
      dicseqend[splline[0]+PeType]=splline


Common=intersect(dicseqend.keys(), dicseqbegin.keys())
Write=open(args.out, 'w')
Head=['SeqName','FlagSamBeg', 'PosInRefBeg','CigarBeg', 'FlagSamEnd', 'PosInRefEnd','CigarEnd', 'NbRep', 'ScoreNewAl','NbRepetNewAl','Seq']
Write.write("\t".join(Head)+'\n')
for PeName in Common :
    pebegin=dicseqbegin[PeName]
    peend=dicseqend[PeName]
    NbRep=(len(peend[9])-(GetInfoFlag(pebegin[5])+GetInfoFlag(peend[5])))/len(args.repet)
    if NbRep >= args.minnbrepet and NbRep<=args.maxnbrepet :
       NewAlignment=DoAlignment(peend[9], args.repet,args.minnbrepet,args.maxnbrepet, around)
       if NewAlignment[0]!="NA" :
         res=[PeName,pebegin[1],pebegin[3],pebegin[5],peend[1],peend[3],peend[5],NbRep, NewAlignment[0],NewAlignment[1],peend[9]]
         Write.write( "\t".join([str(x) for x in res])+'\n')
Write.close()

