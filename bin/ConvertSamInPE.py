#!/usr/bin/python3
# -*- coding: utf-8 -*-
#Version 0.8
import fileinput
import sys
### http://blog.nextgenetics.net/?e=18
##returne vrail si il y a necessite de revese
def intersect(a, b):
     return list(set(a) & set(b))

def difference(a, b):
     return list(set(a) - set(b))

def NeedRev(flag) :
   if flag & 16 or flag & 32:
      return True
   else:
      return False


def NeedPosPair(flag) :
   if flag & 128:
      return 2
   else:
      return 1

def ReadPropPaired(flag) :
    if flag & 3 :
       True
    return False

def ReadUnMapped(flag) :
   if flag & 8 :
      return True
   return False
    
def Rev(x) :
   if x=='A':
      return 'T'
   if x=='T':
      return 'A'
   if x=='C':
      return 'G'
   if x=='G':
      return 'C'
   return x


def RevSeq(Sequence) :
    Seq=""
    for x in Sequence :
      Seq=Rev(x)+Seq
    return Seq


#ListDicPE1[NameSeq]=(FlagSam, Seq, Quality) 
#PrintFormatIllumina(WriteR2,Pe, ListDicPE2[Pe][0], ListDicPE2[Pe][1], ListDicPE2[Pe][2], ListDicPE22Pe][3] )
def PrintFormatIllumina(WriteSortie,NamePE,FlagSam, Seq, Quality):
   WriteSortie.write("@"+NamePE+"\n")
   if len(Seq)!=len(Quality) :
      return
   if NeedRev(int(FlagSam)) :
       WriteSortie.write(RevSeq(Seq)+"\n")
       WriteSortie.write("+\n")
       ### bug detecter 29 avril
       WriteSortie.write(Quality[::-1]+"\n")
   else :
       WriteSortie.write(Seq+"\n")
       WriteSortie.write("+\n")
       WriteSortie.write(Quality+"\n")


if len(sys.argv)!=4 :
   print("Exe SortieR1 SortieR2 ErrorRead")
   sys.exit()


try :
    WriteR1=open(sys.argv[1], 'w')
except :
   print("file "+ sys.argv[1] + " Can't open")

try :
    WriteMonoEnd=open(sys.argv[3], 'w')
except :
   print("file "+ sys.argv[3] + " Can't open")

try :
    WriteR2=open(sys.argv[2], 'w')
except :
   print("file "+ sys.argv[2] + " Can't open")

ListDicPE1={}
ListDicPE2={}

for line in sys.stdin:
    splitligne=line.split()
    FlagSam=int(splitligne[1])
    NameSeq=splitligne[0]
    Seq=splitligne[9]
    Quality=splitligne[10]
    if NeedPosPair(FlagSam)==1 :
       ListDicPE1[NameSeq]=(FlagSam, Seq, Quality) 
    else :
       ListDicPE2[NameSeq]=(FlagSam, Seq, Quality) 

PeCommun=intersect(ListDicPE1.keys(), ListDicPE2.keys())
for Pe in PeCommun :
    PrintFormatIllumina(WriteR1,Pe, ListDicPE1[Pe][0], ListDicPE1[Pe][1], ListDicPE1[Pe][2])
    PrintFormatIllumina(WriteR2,Pe, ListDicPE2[Pe][0], ListDicPE2[Pe][1], ListDicPE2[Pe][2])


PeDiff1=difference(ListDicPE1.keys(), ListDicPE2.keys())
for Pe in PeDiff1:
    PrintFormatIllumina(WriteMonoEnd,Pe, ListDicPE1[Pe][0], ListDicPE1[Pe][1], ListDicPE1[Pe][2])

PeDiff2=difference(ListDicPE2.keys(), ListDicPE1.keys())
for Pe in PeDiff2:
    PrintFormatIllumina(WriteMonoEnd,Pe, ListDicPE2[Pe][0], ListDicPE2[Pe][1], ListDicPE2[Pe][2])

WriteR1.close() 
WriteR2.close() 
WriteMonoEnd.close()

