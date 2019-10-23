#!/usr/bin/env python3
# (c) University of the Witwatersand, Johannesburg on behalf of the H3ABioinformatics Network Consortium
# 2016-2018
# Licensed under the Creative Commons Attribution 4.0 International Licence. 
# See the "LICENSE" file for details

import sys

ListeFile=sys.argv[1::]
Cmt=0
for File in ListeFile :
    Read=open(File)
    line=Read.readline().replace('\n','')
    if Cmt==0:
      print(line) 
    line=Read.readline().replace('\n','')
    if line :
      print(line) 
    Cmt=Cmt+1
    Read.close()
