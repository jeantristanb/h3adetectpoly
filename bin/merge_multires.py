#!/usr/bin/env python3
import sys
import argparse

def parseArguments():
    parser = argparse.ArgumentParser(description='extract sequence of interrest in fasta file')
    parser.add_argument('--listfile',type=str,required=True, help="list of file used for choice best")
    parser.add_argument('--out', type=str,help="out header",default="outseq", required=False)
    parser.add_argument('--head_polycount', type=str,help="header of nb repet of poly to used", required=False, default='NbRepetNewAl')
    parser.add_argument('--head_pe', type=str,help="header of pe to merge", required=False, default='SeqName')
    args = parser.parse_args()
    return args


args=parseArguments()
DicInfoPe={}
ListFile=args.listfile.split(',')
HeadPe=args.head_pe
HeadPCount=args.head_polycount
for File in ListFile :
    Read=open(File)
    Cmt=0
    for line in Read :
      SplL=line.split()
      if Cmt==0 :
         Head=SplL
         PosHead=Head.index(HeadPe)
         PosPCount=Head.index(HeadPCount)
      else :
         if SplL[PosHead] not in DicInfoPe :
            DicInfoPe[SplL[PosHead]]=[int(SplL[PosPCount]),1]
         else :
            DicInfoPe[SplL[PosHead]][1]+=1
            if DicInfoPe[SplL[PosHead]][0]!=int(SplL[PosPCount]) :
              print("warning "+SplL[PosHead]+" different in count "+str(SplL[PosPCount])+" "+SplL[PosPCount])
      Cmt+=1 
    Read.close() 

Write=open(args.out,'w')
Write.write(HeadPe+"\t"+HeadPCount+"\tNbMethod\n")
for Keys in DicInfoPe.keys():
    Write.write(Keys+'\t'+"\t".join([str(x) for x in DicInfoPe[Keys]])+'\n')
Write.close()

