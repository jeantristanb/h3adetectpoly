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

allowed_params =['files_bam', 'bin_samtools', 'out_dir', 'mem_req_samtools', 'cpus_req_samtools', 'reffasta_file', 'out', 'bin_bowtie2', 'bin_bowtie_build']
allowed_params+=['chro', 'pos_begin', 'pos_end', 'polyrep', 'around', 'opt_bowtie']


params.bin_samtools='samtools'
params.bin_bowtie2='bowtie2'
params.bin_bowtie2_build='bowtie2-build'
params.bin_bwa='bwa'

// files contain list of bam to analyseendfas  
params.files_bam=''
params.mem_req_samtools='10GB'
params.cpus_req_samtools=1
params.cpus_req_al=1
params.opt_bowtie="-N 1 -L 20 --rdg 5,5 --rfg 5,5"
params.pe_len=150

params.out_dir='out'
params.out='out'

params.reffasta_file=''
// chro of poly element
params.chro=''
// begin in fasta file of poly element
params.pos_begin=-1
params.pos_end=-1
params.polyrep=''
//nb 
params.around=0
params.around_al2=-1
//nb base pair to check around poly element
params.nb_bp_threshold=5
// Nb repetition of poly element to do
params.reppoly_min=0
params.reppoly_max=100

params.big_time='100h'
params.maxForks_samtools=25


/*check if file_bam exist*/

if(params.files_bam==""){
error('files_bam paramaters not found ')
}
around_al2=params.around_al2
if(around_al2==1)around_al2=params.around

list_bamfile=Channel.fromPath(file(params.files_bam).readLines())

if(params.reffasta_file==''){
error('no reference fasta file give see reffasta_file')
}
if(params.pos_begin<0 || params.pos_end<0 || params.around<0){
error('negative values of pos_end, pos_begin or around')
}
if(params.polyrep==''){
error('element repeted not done')
}
fastafile_ch_sim=Channel.fromPath(params.reffasta_file)
process GenerateSeq{
   input :
     file(filefasta) from fastafile_ch_sim
   publishDir "${params.out_dir}/Seq", overwrite:true, mode:'copy'
   output :
     file("${out}.log")
     file("${out}.insertvar.fasta") into (fasta_sim_var_bowtie, fasta_sim_var_bwa)
     set file("${out}.begin.fasta"), file("${out}.end.fasta") into (fasta_sim_var_bowtie2,fasta_sim_var_bwa2)
   script :
   out=params.out+'_seq'
   """
   extract_seq.py  --fasta $filefasta --chro ${params.chro} --posbegin ${params.pos_begin} --posend ${params.pos_end} --rep ${params.polyrep}  --around ${params.around} --out  $out --minnbrepet ${params.reppoly_min} --maxnbrepet ${params.reppoly_max} --around2 ${params.around_al2}
   """
}

list_bamfile_notmapped=Channel.fromPath(file(params.files_bam).readLines())
process ExtractPENotMapped{
   maxForks params.maxForks_samtools
   time params.big_time
   input :
     file(bam) from list_bamfile_notmapped
   output :
     set val(basenam), file(fileout), file(bam) into sam_not_mapped
   script:
      basenam=bam.baseName
      fileout=bam.baseName+"_umapped.sam"
      """
      ${params.bin_samtools} view  $bam |awk \'{if(\$2==77||\$2==141||\$2==73||\$2==133||\$2==89||\$2==121||\$2==165||\$2==181||\$2==101||\$2==117||\$2==153||\$2==185||\$2==69||\$2==137)print \$0}\'> $fileout
 
      """
}

list_bamfile_region=Channel.fromPath(file(params.files_bam).readLines())
listresbai=file(params.files_bam).readLines()
i=0
while (i < listresbai.size()){
listresbai[i]=listresbai[i].replace('.bam','.bai')
i++
}
list_baifile_region=Channel.fromPath(listresbai)
process ExtractRegion{
   maxForks params.maxForks_samtools
   time params.big_time
   input :
     file(bam) from list_bamfile_region
     file(bai) from list_baifile_region
   output :
     set val(basenam), file(fileout) into sam_region
   script:
      basenam=bam.baseName
      fileout=bam.baseName+"_region.sam"
      Chro=params.chro
      PosBeginAround=params.pos_begin - params.around
      PosEndAround=params.pos_end + params.around
      """
      ${params.bin_samtools} view $bam $Chro:$PosBeginAround-$PosEndAround |awk '{if(\$2!=77&&\$2!=141&&\$2!=73&&\$2!=133&&\$2!=89&&\$2!=121&&\$2!=165&&\$2!=181&&\$2!=101&&\$2!=117&&\$2!=153&&\$2!=185&&\$2!=69&&\$2!=137)print \$0}'>> $fileout
      """
}
// merge
sam_all=sam_not_mapped.join(sam_region)
process MergeSam{
   time params.big_time
   input :
    set val(basenam), file(sam_unmapped), file(bam), file(sam_region) from sam_all
   output :
     file(oubamsort)
     set val(basenam),file(pe1file), file(pe2file), file(uefile) into (pe_ch_bowtie, pe_ch_bwa,pe_ch_bowtie2, pe_ch_bwa2)
   script :
    ousam=bam.baseName+"_sub.sam"
    oubam=bam.baseName+"_tmp.bam"
    oubamsort=bam.baseName+".sort.bam"
    pe1file=bam.baseName+"_pe1.fastq"
    pe2file=bam.baseName+"_pe2.fastq"
    uefile=bam.baseName+"_ue.fastq"
    """
    ${params.bin_samtools} view -H $bam > $ousam
    cat ${sam_unmapped} >> $ousam
    cat ${sam_region} >> ${ousam}
    ${params.bin_samtools} view -S -b ${ousam} >>  ${oubam}
    ${params.bin_samtools} sort -n -o $oubamsort ${oubam}
    ${params.bin_samtools} view $oubamsort | ConvertSamInPE.py $pe1file $pe2file $uefile
    """
}

pe_fas_bowtie_ch=pe_ch_bowtie.combine(fasta_sim_var_bowtie)
process AlignSeqBowtie{
   cpus params.cpus_req_al
   time params.big_time
   input :
    set val(bambase),file(pe1file), file(pe2file), file(uefile), file(ref)  from pe_fas_bowtie_ch
    publishDir "${params.out_dir}/align/bowtie", overwrite:true, mode:'copy'
   output :
    set val(bambase),file(out_bowtie) into bowtie_out_ch
   script :
    refindexhead=ref.baseName
    out_bowtie=bambase+'_alignsim.sam'
    """
    ${params.bin_bowtie2_build} -f $ref $refindexhead
    ${params.bin_bowtie2} -p ${params.cpus_req_al} -x $refindexhead -1 $pe1file -2 $pe2file ${params.opt_bowtie} > ${out_bowtie}
    """
}


pe_fas_bwa_ch=pe_ch_bwa.combine(fasta_sim_var_bwa)
process AlignSeqBWA{
   cpus params.cpus_req_al
   time params.big_time
   input :
    set val(bambase),file(pe1file), file(pe2file), file(uefile), file(ref)  from pe_fas_bwa_ch
    publishDir "${params.out_dir}/align/bwa", overwrite:true, mode:'copy'
   output :
    set val(bambase), file(out_bwa) into bwa_out_ch
   script :
    refindexhead=ref.baseName
    out_bwa=bambase+'_alignsim.sam'
    """
    ${params.bin_bwa}  index $ref
    bwa mem $ref $pe1file $pe2file > ${out_bwa}
    """
}


pe_fas_bowtie2_ch=pe_ch_bowtie2.combine(fasta_sim_var_bowtie2)
process AlignSeqBowtieLocal{
   cpus params.cpus_req_al
   time params.big_time
   input :
    set val(bambase),file(pe1file), file(pe2file), file(uefile), file(refbeg), file(refend)  from pe_fas_bowtie2_ch
    publishDir "${params.out_dir}/align/bowtielocal", overwrite:true, mode:'copy'
  output:
    set val(bambase),file(out_bowtielocalbeg), file(out_bowtielocalend) into bowtielocal_out_ch
  script :
    refindexheadbeg=refbeg.baseName
    out_bowtielocalbeg=bambase+'beg_alignsim.sam'
    refindexheadend=refend.baseName
    out_bowtielocalend=bambase+'end_alignsim.sam'
    """
    ${params.bin_bowtie2_build} -f $refbeg $refindexheadbeg
    ${params.bin_bowtie2} --local -p ${params.cpus_req_al} -x $refindexheadbeg -1 $pe1file -2 $pe2file -L 10 > ${out_bowtielocalbeg}
    ${params.bin_bowtie2_build} -f $refend $refindexheadend
    ${params.bin_bowtie2} --local -p ${params.cpus_req_al} -x $refindexheadend -1 $pe1file -2 $pe2file -L 10 > ${out_bowtielocalend}
    """
}

pe_fas_bwa2_ch=pe_ch_bwa2.combine(fasta_sim_var_bwa2)
process AlignSeqBWALocal{
   cpus params.cpus_req_al
   time params.big_time
   input :
    set val(bambase),file(pe1file), file(pe2file), file(uefile), file(refbeg), file(refend)  from pe_fas_bwa2_ch
    publishDir "${params.out_dir}/align/bwalocal", overwrite:true, mode:'copy'
  output:
    set val(bambase),file(out_bwalocalbeg), file(out_bwalocalend) into bwalocal_out_ch
  script :
    refindexheadbeg=refbeg.baseName
    out_bwalocalbeg=bambase+'beg_alignsim.sam'
    refindexheadend=refend.baseName
    out_bwalocalend=bambase+'end_alignsim.sam'
    pelenlim=params.pe_len-10
    """
    ${params.bin_bwa}  index $refbeg
    bwa mem -U 0 -k 10 -w $pelenlim -H  -p $refbeg $pe1file $pe2file > ${out_bwalocalbeg}
    ${params.bin_bwa}  index $refend
    bwa mem -U 0 -k 10 -w $pelenlim -H  -p $refend $pe1file $pe2file > ${out_bwalocalend}
    """
}
process ComputeStatBWALocal{
   input :
     set val(bambase), file(bwalocalbeg), file(bwalocalend) from bwalocal_out_ch
  publishDir "${params.out_dir}/stats/bwalocal", overwrite:true, mode:'copy'
   output :
      file(out) into stat_bwa_local
  script :
     out=bambase+"_bwalocal.stat"
     """
     extract_seqpoly_reallocalv2.py --sam_begin $bwalocalbeg --sam_end $bwalocalend --out $out --oldpos_poly ${params.around_al2} --nb_bp_threshold ${params.nb_bp_threshold}  --repet ${params.polyrep} --minnbrepet ${params.reppoly_min} --maxnbrepet ${params.reppoly_max}
     """
}

bwalocal_col_ch=bwalocal_out_ch.collect()

process ComputeStatBowtieLocal{
   input :
       set val(bambase), file(bowtielocalbeg), file(bowtielocalend) from bowtielocal_out_ch
        publishDir "${params.out_dir}/stats/bowtielocal", overwrite:true, mode:'copy'
   output :
       set val(bambase),file(out) into stat_bowtie_local
   script :
       out=bambase+"_bowtielocal.stat"
       """
       extract_seqpoly_reallocalv2.py --sam_begin $bowtielocalbeg --sam_end $bowtielocalend --out $out --oldpos_poly ${params.around_al2} --nb_bp_threshold ${params.nb_bp_threshold}  --repet ${params.polyrep} --minnbrepet ${params.reppoly_min} --maxnbrepet ${params.reppoly_max}
       """
}

process ComputeStatBWA{
   input :
     set val(bambase), file(bwaout) from bwa_out_ch
  publishDir "${params.out_dir}/stats/bwa", overwrite:true, mode:'copy'
   output :
      set val(bambase),file(out) into stat_bwa
  script :
     out=bambase+"_bwa.stat"
     """
     extract_seqpoly_real.py --sam $bwaout --out $out --repet ${params.polyrep} --minnbrepet  ${params.reppoly_min} --maxnbrepet ${params.reppoly_max} --nb_bp_threshold ${params.nb_bp_threshold}
     """
}


process ComputeStatBowtie{
   input :
     set val(bambase), file(bowtieout) from bowtie_out_ch
  publishDir "${params.out_dir}/stats/bowtie", overwrite:true, mode:'copy'
   output :
      set val(bambase),file(out) into stat_bowtie
  script :
     out=bambase+"_bowtie.stat"
     """
     extract_seqpoly_real.py --sam $bowtieout --out $out --repet ${params.polyrep} --minnbrepet  ${params.reppoly_min} --maxnbrepet ${params.reppoly_max} --nb_bp_threshold ${params.nb_bp_threshold}
     """
}




