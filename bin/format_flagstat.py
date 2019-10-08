#!/usr/bin/env python3
import sys
import os
import argparse

ListFile=sys.argv[1::]
print("\t".join(["FileName","NbPeI", "NbSeqMap", "NbPeWellMap"]))
for File in ListFile:
   FileB=os.path.basename(File) 
   ReadF=open(File) 
   AllRead=ReadF.readlines()
   nbSeqTot=int(AllRead[0].split('+')[0])/2
   nbSeqMap=int(AllRead[4].split('+')[0])
   nbPEPropMap=int(AllRead[8].split('+')[0])/2
   print("\t".join([str(x) for x in [FileB, nbSeqTot,nbSeqMap,nbPEPropMap]]))
