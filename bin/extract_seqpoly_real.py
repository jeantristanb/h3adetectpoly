#!/usr/bin/env python3
import sys
import argparse
from Bio import pairwise2
from Bio.pairwise2 import format_alignment




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
    parser.add_argument('--sam',type=str,required=True, help="fasta file")
    parser.add_argument('--out', type=str,help="out header",default="outseq", required=False)
    parser.add_argument('--minnbrepet', type=int,help="how many position around", required=False, default=0)
    parser.add_argument('--maxnbrepet', type=int,help="how many position around", required=False, default=100)
    parser.add_argument('--repet',type=str,required=False, help="fasta file")
    parser.add_argument('--nb_bp_threshold', type=int,help="out header",default=5, required=False)
    args = parser.parse_args()
    return args


args=parseArguments()
readsam=open(args.sam)
Write=open(args.out, 'w')
around=args.nb_bp_threshold

Head=['SeqName','NbRepRef', 'BegPolyRef', 'EndPolyRef', 'BegSeqRef', 'EndSeqRef', 'Cigar', 'FlagSam', 'ScoreNewAl','NbRepetNewAl','Seq']
Write.write("\t".join(Head)+'\n')
for line in readsam :
    if line[0]=='@':
      continue
    splline=line.split('\t')
    FlagSam=int(splline[1])
    if ReadUnMapped(FlagSam) :
       continue
    Seq=splline[2].split('_')
    NbRep=int(Seq[0].replace('Nb',''))
    BegPoly=int(Seq[1])
    EndPoly=int(Seq[2])
    ###  positon in sequence 
    Cigar=splline[5]

    BegInSeq=int(splline[3])
    EndInSeq=BegInSeq+len(splline[9])
    PeType="Pe1"
    if IsSecondPE(FlagSam) :
       PeType="Pe2"
    splline[0]=splline[0]+PeType
    if BegInSeq<(BegPoly-around) and EndInSeq> (EndPoly+around) :
       NewAlignment=DoAlignment(splline[9], args.repet,args.minnbrepet,args.maxnbrepet, around)
       if NewAlignment[0]!="NA" :
         Write.write( "\t".join([str(x) for x in [splline[0],NbRep, BegPoly, EndPoly, BegInSeq, EndInSeq, Cigar, FlagSam, NewAlignment[0], NewAlignment[1],splline[9]]])+'\n')

Write.close()
