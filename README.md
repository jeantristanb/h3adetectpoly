# Len of sequences repeted in genomes

## Background : 

Todo

## 2.2 Pre-requisites
python, nextflow, samtools  



### Parameters 
* `files_bam` : files contains each bam file to analyse, bai must be in same folder and obtain by replace bam extension with bai
* `bin_samtools` : binary for samtools [default : samtools]
* `out_dir` : out direction 
* `mem_req_samtools` : memory request for samtools 
* `cpus_req_samtools` : cpu request for samtools
* `reffasta_file` : reference fasta filo 
s_begin', 'pos_end', 'polyrep', 'around'
* `out` : out pattern
* `chro` : chromosome on reference where poly is found 
* `pos_begin` : pos begin on poly element
* `pos_end` : pos end on poly element
* `polyrep` : what element is repeated : for instance if in sequence is acacac give as paramater ac 
* `around` : how many position must be analysed around poly element


## nf : Check, computed stat on bam file, fasta file in function of poly-tc and 
### What's done 
1. computed with flagstat percentage of sequence well aligned, PE well aligned
2. computed depth around and in poly-sequence
3. check poly element in reference sequence

### What output :
1. file `$out.depth` : contains depth for each individuals and position around  poly element
2. file `$out.statalignhall` : contains  for each bam files PE Initial, Sequence aligned and PE well aligned. : * if all PE is well aligned, your data have been too much cleaned *
3. file `$out.statfasta.log` : contains analyse of fasta file, polyrep etc...  : 
  * if "Error Message" is not empty that mean error in some parameters and must be checker, otherwise algorithms didn't run
  * if "Warning Message" is not empty that mean need to be checked

### 
