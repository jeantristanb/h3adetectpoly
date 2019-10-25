# Detection of number repeted element in whole genome sequencing data

## Background : 
Objectif is to detected in  specific position of genome repetion number of polyelement using WGS paired end data. Using data already aligned with software (as bwa with option mem), pipeline extract all Paired end should be contains element, realigned with different approaches and after estimate poly element in genome. After we estimated if individuals is homozygote, heterozygote and size have allele

## 2.2 Pre-requisites
  * python3 library : fileinput, BioPython, argparse
  * R library : optparse, vcd
  * nextflow, samtools, bwa, bowtie2


### Parameters 
* *Input and Output* :
 * `files_bam` : files contains each bam file to analyse, bai must be in same folder and obtain by replace bam extension with bai
 * `reffasta_file` : reference fasta filo 
 * `out` : out pattern
 * `out_dir` : out direction 

* *Poly element information* :
 * `chro` : chromosome on reference where poly is found 
 * `pos_begin` : pos begin on poly element
 * `pos_end` : pos end on poly element
 * `polyrep` : what element is repeated : for instance if in sequence is acacac give as paramater ac 

* *PE information* :
  * `pe_len` : len of PE in intial data set [default : 150 ]

* *Parameters of analysis* :
 * `around` : how many position must be analysed around poly element and extract of reference sequences [default : 0]
 * `reppoly_min` : minimum of repetition of poly element T inserted in reference sequences used for _repetition_ procedure default [default : 0]
 * `reppoly_max=100 : maximum of repetition of poly element T inserted in reference sequences used for _repetition_ procedure default [default : 100]
 * `nb_bp_threshold` : minimum position keept in sequences realigned around element 

* *alignment option* :
 * bowtie2 : 
  * `opt_bowtie` : what parameter for local alignment : [default : -N 1 -L 20 --rdg 5,5 --rfg 5,5]

* *binaries software *:
 * `bin_samtools` : binary for samtools [default : samtools]
 * `bin_bowtie2`  : bowtie2 software (alignment) [default : bowtie2 ]
 * `bin_bowtie2_build`  : bowtie2 build reference software (alignment) [default : bowtie2-build ]
 * `bin_bwa` : [default : bwa]
* *memory and cpu allocation* :
 * `cpus_req_al` : cpu request for alignment 
 * `mem_req_samtools` : memory request for samtools 
 * `cpus_req_samtools` : cpu request for samtools

## `detectpoly.nf` : count of polyelement number in genome
main

### What's done
1. Extract Paired end of interrest (_Peinterest_) :
  1. around of region of interest : see args `pos_begin`, `pos_end` and `around`
  2. PE not well aligned using FlagSam : 77,141,73,133,89,121,165,181,101,117,153,185,69,137
  3. sequences will be convert in paired end (fastq file)

2. Doing a new alignment with reference fasta before and after poly element (procedure call _local_) :
  1. two new reference sequences with selection `around` positions around element repeted, deleted repeted element and build one fasta with sequences before alignment and one fasta with sequences after alignment
  2. align two times Pe with local procedure first on fasta before alignment and after alignment: `bwa mem` and `bowtie2` (option --local and see argument `opt_bowtie`)
  3. using both alignment with breakpoint to defined size of poly element with realignment using `bio python`

3. Doing a new alignment with n insertion of poly element in reference sequences and used best alignment to defined size element
  1. New reference fasta files contains n sequences with poly T element and `around` position around poly T repeted between `reppoly_min` and `reppoly_max` times
  2. align two times Pe with local procedure first on fasta before alignment and after alignment: `bwa mem` and `bowtie2` (global)
  3. using both alignment with breakpoint to defined size of poly element with realignment using `bio python`

4. Merge elements 

### Output :
  * stats folder
  * resume statistics 
  *

## `checkdata.nf` : Check, computed stat on bam file, fasta file in function of poly-tc and 
### What's done 
1. computed with flagstat percentage of sequence well aligned, PE well aligned
2. computed depth around and in poly-sequence
3. check poly element in reference sequence

### What output :
  1.
  3. using both alignment with breakpoint to defined size of poly element with realignment using `bio python`

3. Doing a new alignment with n insertion of poly element in reference sequences and used best alignment to defined size element

## `checkdata.nf` : Check, computed stat on bam file, fasta file in function of poly-tc and 
### What's done 
1. computed with flagstat percentage of sequence well aligned, PE well aligned
2. computed depth around and in poly-sequence
3. check poly element in reference sequence

### What output :

## `checkdata.nf` : Check, computed stat on bam file, fasta file in function of poly-tc and 
### What's done 
1. computed with flagstat percentage of sequence well aligned, PE well aligned
2. computed depth around and in poly-sequence
3. check poly element in reference sequence

### What output :
1. file `$out.depth` : contains depth for each individuals and position around  poly element
2. file `$out.statalignhall` : contains  for each bam files PE Initial, Sequence aligned and PE well aligned. : * if all PE is well aligned, your data have been too much cleaned *
3. file `$out.statfasta.log` : contains analyse of fasta file, polyrep etc...  : 
  1. if "Error Message" is not empty that mean error in some parameters and must be checker, otherwise algorithms didn't run
  2. if "Warning Message" is not empty that mean need to be checked

### 
