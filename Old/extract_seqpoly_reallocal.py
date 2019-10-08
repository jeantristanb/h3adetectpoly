#!/usr/bin/env python3
import sys
import argparse

ListeFlag=set(['I', 'M', 'D', 'N', 'S'])
def GetInfoFlag(FlagSam):
   Chaine=""
   Longueur=0
   for x in FlagSam :
      if x in ListeFlag :
          if x!='I' and x!='S':
              Longueur+=int(Chaine)
          Chaine=""
      else :
          Chaine+=x
   return Longueur



def ReadUnMapped(flag) :
   if flag & 8 :
      return True
   return False


def parseArguments():
    parser = argparse.ArgumentParser(description='extract sequence of interrest in fasta file')
    parser.add_argument('--fasta',type=str,required=False, help="fasta file")
    parser.add_argument('--sam',type=str,required=True, help="fasta file")
    parser.add_argument('--out', type=str,help="out header",default="outseq", required=False)
    parser.add_argument('--oldpos_poly', type=int,help="out header", required=True)
    parser.add_argument('--around', type=int,help="out header", required=False, default=5)
    args = parser.parse_args()
    return args


args=parseArguments()
readsam=open(args.sam)
Write=open(args.out+"_detailseq.out", 'w')
around=args.around
oldpospoly=args.oldpos_poly

for line in readsam :
    if line[0]=='@':
      continue
    splline=line.split('\t')
    if ReadUnMapped(int(splline[1])) :
       continue
    ###  positon in sequence 
    Cigar=splline[5]
    lenseq=GetInfoFlag(Cigar)
    BegInSeq=int(splline[3])
    EndInSeq=BegInSeq+lenseq
    print(BegInSeq,EndInSeq, Cigar)
    BaliseBegin=BegInSeq>=(oldpospoly-around) and BegInSeq<=(oldpospoly+around) and 'S' in Cigar
    BaliseEnd=EndInSeq>=(oldpospoly-around) and EndInSeq<=(oldpospoly+around) and 'S' in Cigar
    if BaliseBegin or BaliseEnd :
       Write.write( "\t".join([str(x) for x in [args.sam,oldpospoly, BegInSeq, EndInSeq, splline[5], splline[0], splline[1], splline[9]]])+'\n')

Write.close()
