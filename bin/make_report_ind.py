#!/usr/bin/env python3

# (c) University of the Witwatersand, Johannesburg on behalf of the H3ABioinformatics Network Consortium
# 2016-2018
# Licensed under the MIT Licence
# See the "LICENSE" file for details

import glob
import sys
import os
from  subprocess import CalledProcessError
import subprocess
import re

def check_output(x,shell=False):
   ans=subprocess.check_output(x,shell=shell)
   ans=str(ans,'ascii') # encoding option only frmo python 3.6
   return ans


def parseArguments():
    parser = argparse.ArgumentParser(description='build report for each individuals')
    parser.add_argument('--order',type=str,required=False, help="fasta file")
    parser.add_argument('--indname',type=str,required=False, help="fasta file")
    parser.add_argument('--figures',type=str,required=True, help="fasta file")
    parser.add_argument('--filestat',type=str,required=True, help="fasta file")
    parser.add_argument('--out', type=str,help="out header",default="outseq", required=False)
    args = parser.parse_args()
    return args



try:
   kpsewhich=check_output("which kpsewhich",shell=True)
   if kpsewhich:
      kpsewhich=check_output("kpsewhich datetime.sty",shell=True)
except CalledProcessError:
   kpsewhich=""

fancy="""
*-usepackage{fancyhdr}
*-usepackage[yyyymmdd,hhmmss]{datetime}
*-pagestyle{fancy}
*-rfoot{Completed on *-today*- at *-currenttime}
*-cfoot{}
*-lfoot{Page *-thepage}
"""

try:
   kpsewhich=check_output("which kpsewhich",shell=True)
   if kpsewhich:
      kpsewhich=check_output("kpsewhich datetime.sty",shell=True)
except CalledProcessError:
   kpsewhich=""

fancy="""
*-usepackage{fancyhdr}
*-usepackage[yyyymmdd,hhmmss]{datetime}
*-pagestyle{fancy}
*-rfoot{Completed on *-today*- at *-currenttime}
*-cfoot{}
*-lfoot{Page *-thepage}
"""

dateheader=""
if len(kpsewhich)>1:
   dfmt = kpsewhich.rstrip()
   if os.access(dfmt,os.R_OK):
      with open(dfmt) as f:
         for line in f:
            m=re.search("ProvidesPackage.datetime..(..../../..)",line)
            if m and m.group(1) >= "2010/09/21":
               dateheader=fancy

if len(sys.argv)<=1:
   sys.argv = "make_assoc_report.py ${params.pheno} $texf".split()


Individu  = args.indname
out    = args.out


