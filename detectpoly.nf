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
params.mem_req_samtools='20GB'
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
params.around_depth=500
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
   memory params.mem_req_samtools
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
      set val(bambase),file(filedistbwalocal) into fig_bwalocal
      set val(bambase),file(outresume) into (stat_bwalocal, stat_bwalocal2)
      set val(bambase),file(out) into stat_bwalocal_merge
  script :
     out=bambase+"_bwalocal.stat"
     outresume=bambase+"_bwalocal_resume.stat"
     filedistbwalocal=bambase+"_bwalocal.pdf"
     """
     extract_seqpoly_reallocalv2.py --sam_begin $bwalocalbeg --sam_end $bwalocalend --out $out --oldpos_poly ${params.around_al2} --nb_bp_threshold ${params.nb_bp_threshold}  --repet ${params.polyrep} --minnbrepet ${params.reppoly_min} --maxnbrepet ${params.reppoly_max}
     ExtractAllele_V2.r --file $out --out $bambase"_bwalocal" --header $bambase
     """
}



process ComputeStatBowtieLocal{
   input :
       set val(bambase), file(bowtielocalbeg), file(bowtielocalend) from bowtielocal_out_ch
        publishDir "${params.out_dir}/stats/bowtielocal", overwrite:true, mode:'copy'
   output :
      set val(bambase),file(filedistbowtielocal) into fig_bowtielocal
      set val(bambase), file(outresume) into (stat_bowtielocal, stat_bowtielocal2)
      set val(bambase),file(out) into stat_bowtielocal_merge
   script :
       out=bambase+"_bowtielocal.stat"
      outresume=bambase+"_bowtielocal_resume.stat"
      filedistbowtielocal=bambase+"_bowtielocal.pdf"
       """
       extract_seqpoly_reallocalv2.py --sam_begin $bowtielocalbeg --sam_end $bowtielocalend --out $out --oldpos_poly ${params.around_al2} --nb_bp_threshold ${params.nb_bp_threshold}  --repet ${params.polyrep} --minnbrepet ${params.reppoly_min} --maxnbrepet ${params.reppoly_max}
     ExtractAllele_V2.r --file $out --out $bambase"_bowtielocal" --header $bambase
       """
}

process ComputeStatBWA{
   input :
     set val(bambase), file(bwaout) from bwa_out_ch
  publishDir "${params.out_dir}/stats/bwa", overwrite:true, mode:'copy'
   output :
      set val(bambase),file(filedistbwa) into fig_bwa
      set val(bambase), file(outresume) into (stat_bwa, stat_bwa2)
      set val(bambase),file(out) into stat_bwa_merge
  script :
     out=bambase+"_bwa.stat"
     outresume=bambase+"_bwa_resume.stat"
     filedistbwa=bambase+"_bwa.pdf"
     """
     extract_seqpoly_real.py --sam $bwaout --out $out --repet ${params.polyrep} --minnbrepet  ${params.reppoly_min} --maxnbrepet ${params.reppoly_max} --nb_bp_threshold ${params.nb_bp_threshold}
     ExtractAllele_V2.r --file $out --out $bambase"_bwa" --header $bambase
     """
}


process ComputeStatBowtie{
   input :
     set val(bambase), file(bowtieout) from bowtie_out_ch
  publishDir "${params.out_dir}/stats/bowtie", overwrite:true, mode:'copy'
   output :
      set val(bambase),file(filedistbowtie) into fig_bowtie
      set val(bambase), file(outresume) into (stat_bowtie, stat_bowtie2)
      set val(bambase),file(out) into stat_bowtie_merge
  script :
     out=bambase+"_bowtie.stat"
     outresume=bambase+"_bowtie_resume.stat"
     filedistbowtie=bambase+"_bowtie.pdf"
     """
     extract_seqpoly_real.py --sam $bowtieout --out $out --repet ${params.polyrep} --minnbrepet  ${params.reppoly_min} --maxnbrepet ${params.reppoly_max} --nb_bp_threshold ${params.nb_bp_threshold}
     ExtractAllele_V2.r --file $out --out $bambase"_bowtie" --header $bambase
     """
}
join_res=stat_bowtie_merge.join(stat_bowtielocal_merge).join(stat_bwa_merge).join(stat_bwalocal_merge)
process ComputeStatAll{
   input :
      set val(bambase), file(bowtie), file(bowtieloc), file(bwa), file(bwalocal) from join_res
  publishDir "${params.out_dir}/stats/all", overwrite:true, mode:'copy'
   output :
      set val(bambase),file(filedistall) into fig_all
      set val(bambase),file(outresume) into (stat_all, stat_all2)
      set val(bambase),file(out) into stat_all_merge
   script :
     out=bambase+"_all.stat"
     outresume=bambase+"_all_resume.stat"
     filedistall=bambase+"_all.pdf"
      """
     merge_multires.py --listfile $bwalocal,$bwa,$bowtie,$bowtieloc --out $out
     ExtractAllele_V2.r --file $out --out $bambase"_all" --header $bambase --lheader NbRepetNewAl
      """
}

stat_all_col=stat_all.flatMap{n->n[1]}.collect()
process MergeStatAll{
    input :
        file(statmerge) from stat_all_col
    publishDir "${params.out_dir}/stats", overwrite:true, mode:'copy'
    output :
      file(out) into merge_all
    script :
      mergestat=statmerge.join(" ")
      out="resume.all"
      """
      mergefile.py $mergestat > $out
      """
}

stat_bwa_col=stat_bwa.flatMap{n->n[1]}.collect()
process MergeStatBWA{
    input :
        file(statmerge) from stat_bwa_col
    publishDir "${params.out_dir}/stats", overwrite:true, mode:'copy'
    output :
      file(out) into merge_bwa
    script :
      mergestat=statmerge.join(" ")
      out="resume.bwa"
      """
      mergefile.py $mergestat > $out
      """
}

stat_bowtielocal_col=stat_bowtielocal.flatMap{n->n[1]}.collect()
process MergeStatBowtieLocal{
    input :
        file(statmerge) from stat_bowtielocal_col
    publishDir "${params.out_dir}/stats", overwrite:true, mode:'copy'
    output :
      file(out) into merge_bowtieloc
    script :
      mergestat=statmerge.join(" ")
      out="resume.bowtielocal"
      """
      mergefile.py $mergestat > $out
      """
}


stat_bwalocal_col=stat_bwalocal.flatMap{n->n[1]}.collect()
process MergeStatBWALocal{
    input :
        file(statmerge) from stat_bwalocal_col
    publishDir "${params.out_dir}/stats", overwrite:true, mode:'copy'
    output :
      file(out) into merge_BWALoc
    script :
      mergestat=statmerge.join(" ")
      out="resume.bwalocal"
      """
      mergefile.py $mergestat > $out
      """
}



stat_bowtie_col=stat_bowtie.flatMap{n->n[1]}.collect()
process MergeStatBowtie{
    input :
        file(statmerge) from stat_bowtie_col
    publishDir "${params.out_dir}/stats", overwrite:true, mode:'copy'
    output :
      file(out) into merge_bowtie
    script :
      mergestat=statmerge.join(" ")
      out="resume.bowtie"
      """
      mergefile.py $mergestat > $out
      """
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
   publishDir "${params.out_dir}/depth/", overwrite:true, mode:'copy'
   output :
      file(out) into depth_all
      file(outpdf) 
   script :
     out=params.out+".depth"
     outpdf=params.out+"_depth.pdf"
     range=params.chro+':'+(params.pos_begin-params.around_depth)+'-'+(params.pos_end+params.around_depth)
     listfile=allbam.join(' ')
     Head='Chro Pos '+allbam.join('\\t').replace('bam','')
     """
     echo -e \"$Head\" > $out
     samtools depth -aa -r $range $listfile >> $out
     plot_couv.r $out  ${params.pos_begin} ${params.pos_end} $outpdf
     """
}

fig_merge_process=fig_bwa.join(fig_bwalocal).join(fig_bowtie).join(fig_bowtielocal).join(fig_all)
stat_merge_process=stat_bwa2.join(stat_bwalocal2).join(stat_bowtie2).join(stat_bowtielocal2).join(stat_all2)
report_ind=fig_merge_process.join(stat_merge_process).combine(depth_all)

process DoReportInd{
   input :
       set val(bambase), file(figbwa), file(figbwalocal), file(figbowtie), file(figbowtielocal), file(figall),  file(statbwa),file(statbwalocal),file(statbowtie),file(statbowtielocal) ,file(statall), file(depth_all) from report_ind
  publishDir "${params.out_dir}/report_ind", overwrite:true, mode:'copy'
   output :
     file("${bambase}.tex")
     file("${bambase}.pdf")
   script :
     outpdfdepth=bambase+"_depth.pdf"
     """
     plot_couv_ind.r ${depth_all} ${params.pos_begin} ${params.pos_end} $outpdfdepth $bambase"."
     make_report_ind.py --out $bambase".tex" --order "all,bwa,bwalocal,bowtie2,bowtie2local" --indname $bambase --figures $figall,$figbwa,$figbwalocal,$figbowtie,$figbowtielocal --filestat $statall,$statbwa,$statbwalocal,$statbowtie,$statbowtielocal --couv  $outpdfdepth
     """

}


