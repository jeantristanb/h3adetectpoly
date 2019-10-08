#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *      Jean-Tristan Brandenburg
 *
 * Description Nextflow Pipeline for :
 *  * Do statistics on bam file and check if not well aligned sequences is present
 *  * Computed 
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths


def helps = [ 'help' : 'help' ]

allowed_params =['files_bam', 'bin_samtools', 'out_dir', 'mem_req_samtools', 'cpus_req_samtools', 'reffasta_file', 'out']
allowed_params+=['chro', 'pos_begin', 'pos_end', 'polyrep', 'around']


params.bin_samtools='samtools'
params.files_bam=''
params.mem_req_samtools='10GB'
params.cpus_req_samtools=1

params.out_dir='out'
params.out='out'

params.reffasta_file=''
params.chro=''
params.pos_begin=-1
params.pos_end=-1
params.polyrep=''
params.around=0
params.maxForks_samtools=25


/*check if file_bam exist*/

if(params.files_bam==""){
error('files_bam paramaters not found ')
}


list_bamfile=Channel.fromPath(file(params.files_bam).readLines())

process ComputeStatBam{
   memory params.mem_req_samtools
   cpus params.cpus_req_samtools
   maxForks params.maxForks_samtools
   input :
      file(bamfile) from list_bamfile
   output:
     file(statfile) into statbam
   script :
      statfile=bamfile.baseName+".stat"
   """
   samtools flagstat -@ ${params.cpus_req_samtools} $bamfile > $statfile
   """
}
allstatbam=statbam.collect()
process FormatStatBam{
   input :
      file(stat) from allstatbam
   publishDir "${params.out_dir}/stat", overwrite:true, mode:'copy'
   output :
       file("$out")
   script :
      out=params.out+'.statalignhall'
      AllFile=stat.join(' ')
   """
   format_flagstat.py $AllFile > $out
   """
}


if(params.reffasta_file==''){
error('no reference fasta file give see reffasta_file')
}
if(params.pos_begin<0 || params.pos_end<0 || params.around<0){
error('negative values of pos_end, pos_begin or around')
}
if(params.polyrep==''){
error('element repeted not done')
}
fastafile_ch_check=Channel.fromPath(params.reffasta_file)
process CheckRefSeq{
   input :
     file(filefasta) from fastafile_ch_check
   publishDir "${params.out_dir}/stat", overwrite:true, mode:'copy'
   output :
     file("${out}.log")
   script :
   out=params.out+'.statfasta'
   """
   extract_seq.py  --fasta $filefasta --chro ${params.chro} --posbegin ${params.pos_begin} --posend ${params.pos_end} --rep ${params.polyrep}  --around ${params.around} --check --out  $out
   """
}

list_bamfile_check=Channel.fromPath(file(params.files_bam).readLines())
process CheckChro{
   maxForks params.maxForks_samtools
   input :
     file(bam) from list_bamfile_check
   script :
      """
      samtools view -H $bam |CheckChro.py ${params.chro}
      """
}


if(1==2){
list_bamfile_bai=Channel.fromPath(file(params.files_bam).readLines())
process doBai{
   memory params.mem_req_samtools
   cpus params.cpus_req_samtools
   input :
     file(bam) from list_bamfile_bai
   output :
     file('${bam}.bai') into list_bai_file
   script :
     """
     samtools index -@ ${params.cpus_req_samtools}  $bam
     """
}
}

list_bamfile_depth=Channel.fromPath(file(params.files_bam).readLines()).collect()
listresbai=file(params.files_bam).readLines()
i=0
while (i < listresbai.size()){
listresbai[i]=listresbai[i].replace('.bam','.bai')
i++
}
list_baifile_depth=Channel.fromPath(listresbai).collect()

process ComputeDepth{
   input :
      file(allbam) from list_bamfile_depth
      file(allbai) from list_baifile_depth
   publishDir "${params.out_dir}/stat/", overwrite:true, mode:'copy'
   output :
      file(out)
   script :
     out=params.out+".depth"
     outpdf=params.out+"_depth.pdf"
     range=params.chro+':'+(params.pos_begin-params.around)+'-'+(params.pos_end+params.around)
     listfile=allbam.join(' ') 
     Head='Chro Pos '+allbam.join('\t').replace('bam','')
     """
     echo $Head > $out
     samtools depth -aa -r $range $listfile >> $out
     plot_couv.r $out  ${params.pos_begin} ${params.pos_end} $outpdf
     """
}



