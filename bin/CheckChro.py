#!/usr/bin/env python3
# (c) University of the Witwatersand, Johannesburg on behalf of the H3ABioinformatics Network Consortium
# 2016-2018
# Licensed under the Creative Commons Attribution 4.0 International Licence. 
# See the "LICENSE" file for details

import sys

chro=sys.argv[1]
read=sys.stdin
Balise=False
#@SQ	SN:3	LN:198022430
ListChro=[]
for line in read :
   if line[0:3]=='@SQ' :
      Chro=line.split('SN:')[1].split()[0]   
      if Chro==chro :
         Balise=True
      ListChro.append(Chro)
   
if Balise==False :
   print('Chro not found in bam file '+chro)
   print(' '.join(ListChro))
   sys.exit(1)


