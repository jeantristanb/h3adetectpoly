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
template='''
*-documentclass[11pt]{article}

*-usepackage[paper=a4paper,left=2cm,right=2cm,top=2cm,bottom=2cm]{geometry}
*-usepackage{graphicx}
*-usepackage{subfig}
*-usepackage{listings}
*-usepackage{longtable}
*-usepackage{array}
*-usepackage{booktabs}
*-usepackage{float}
*-usepackage{dcolumn}
*-floatstyle{ruled}
*-restylefloat{figure}
*-restylefloat{table}
*-newcommand{*-lefttblcol}{*-raggedright*-hspace{0pt}}
*-newcommand{*-righttblcol}{*-raggedleft*-hspace{0pt}}
*-newcommand{*-centretblcol}{*-centering*-hspace{0pt}}

*-newcolumntype{P}[1]{>{*-lefttblcol}p{#1}}
*-newcolumntype{Q}[1]{>{*-righttblcol}p{#1}}
*-newcolumntype{R}[1]{>{*-centretblcol}p{#1}}
*-lstset{
basicstyle=*-small*-ttfamily,
columns=flexible,
breaklines=true
}
*-usepackage{url}
*-title{Association Testing  %(base)s : %(pheno)s}
*-date{%(date)s}
'''+dateheader+('''
*-author{SBIMB : detection of Poly element }

*-newcommand{*-ourfig}[3]{*-begin{figure}[ht]*-begin{center}*-includegraphics[scale=0.6]{#3} *-end{center} *-caption{#2 [File is *-protect*-url{#3}]}  *-label{#1}*-end{figure}}
*-begin{document}

*-maketitle

*-section{Introduction}

This report gives a brief overview of the run of the detection of poly element in genomes
*-begin{itemize}
*-item You were testing for the following individus %s 
*-end{itemize}
'''%Individu)

EOL = chr(10)


images = "${workflow.container}"
if images=="[:]":
   pdict["dockerimages"] = ": locally installed binaries used"
else:
   images = getImages("${workflow.container}")
   pdict["dockerimages"] = ": the docker images used are found in "+images


template = template % pdict



template=template.replace("*-",chr(92)).replace("##",chr(36)).replace("@.@",chr(10))

g=open(out,"w")
g.write(template)
g.close()

os.system("pdflatex %s >& /dev/null"%out)
os.system("pdflatex %s"%out)




