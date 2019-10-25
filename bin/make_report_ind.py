#!/usr/bin/env python3

# (c) University of the Witwatersand, Johannesburg on behalf of the H3ABioinformatics Network Consortium
# 2016-2018
# Licensed under the MIT Licence
# See the "LICENSE" file for details

import glob
import argparse
import sys
import os
from  subprocess import CalledProcessError
import subprocess
import re

def ltx(x):
   x=x.replace('_','\\_')
   return(x)
   
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

args=parseArguments()
Individu=args.indname


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
*-usepackage{grffile}
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
*-title{Poly detection of element : %(Individu)s}
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

listtype=args.order.replace(' ','').split(',')
liststat=args.filestat.replace(' ','').split(',')
listfigure=args.figures.replace(' ','').split(',')

ListRes={}
ListPosHead=[]
for Cmt in  range(len(liststat)) :
    Read=open(liststat[Cmt]) 
    ListHead=Read.readline().replace('\n','').split()
    ListInfo=Read.readline().replace('\n','').split()
    CmtHead=0
    for Head in ListHead :
      if Head!="Name" and "Type" not in Head:
        if Head not in ListRes : 
            ListRes[Head]=['-']*len(liststat)
            ListPosHead.append(Head)
        ListRes[Head][Cmt]=ListInfo[CmtHead]
      CmtHead+=1


### 
TableRes="\n*-begin{table}[ht]\n*-centering\n*-resizebox{*-textwidth}{!}{\n*-begin{tabular}{ |"+"c|"*(len(ListRes)+1)+"| }"
TableRes+="\n*-hline\nType\nAlignement & "+" & ".join([x.replace('_',' ') for x in ListPosHead]).replace('NbRepet', '')+"  *-*- \n*-hline\n"
for Cmt in range(len(liststat)):
   resind=[listtype[Cmt]]
   resind+=[ListRes[Head][Cmt] for Head in ListPosHead]
   TableRes+=" & ".join(resind)+" *-*- \n "
TableRes+="*-hline\n*-end{tabular}\n}\n*-caption{resume for each alignment nb PE found for each Allele}\n*-end{table} \n"

FigRes=""
for CmtFig in range(len(listfigure)):
   Figure=listfigure[CmtFig]
   #NewFigure=Figure.replace('.','-').replace('-svg', '.png')
   #print('convert '+ Figure+' '+NewFigure)
   #os.system('convert '+ Figure+' '+NewFigure)
   FigRes+="\subsection{"+listtype[CmtFig]+"}\n"
   FigRes+="File used for figure "+ltx(Figure)
   FigRes+="\n*-begin{figure}[h]\n*-includegraphics[width=*-linewidth]{"+listfigure[CmtFig]+"}\n*-caption{ "+ listtype[CmtFig] + " }\n*-label{fig:"+listtype[CmtFig]+"}\n*-end{figure}\n"
  

template+="\n*-section{Resume results}\n"+TableRes+"\n"
template+="\n*-section{Figures distribution}\n"+FigRes+"\n"

template=template.replace("*-",chr(92)).replace("##",chr(36)).replace("@.@",chr(10))
template+="\end{document}\n"
hashd = { 'Individu':Individu.replace("_","-"), "date":check_output("date").strip()}
template=template%hashd

   


out=args.out
g=open(out,"w")
g.write(template)
g.close()

os.system("pdflatex %s >& /dev/null"%out)
os.system("pdflatex %s"%out)




