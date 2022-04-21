

##1, QC
study_name=2022-03-25_AAV_SEQ
MinIONQC_program=/drive2/wli/software2/biosoft/seq_ana/minion_qc/MinIONQC.R
seq_dir=/home/Research_Archive/RawData/Gene_Therapy/Nanopore/$study_name/
seq_sum_f=$seq_dir/sequencing_summary_FAR99886_eebaec55.txt #sequencing_summary_FAS00280_baea760c.txt #sequencing_summary_FAR52309_4f2c1716.txt #sequencing_summary_FAQ39952_ff966d87.txt

out_dir=/home/Research_Archive/ProcessedData/Nanopore_seq/$study_name/MinIONQC/
mkdir -p $out_dir
cd $out_dir

/bin/Rscript $MinIONQC_program  -i $seq_sum_f -o $out_dir -s TRUE


##2, demultiplexing:
study_name=2022-03-25_AAV_SEQ

mkdir -p /home/Research_Archive/ProcessedData/Nanopore_seq/$study_name/fastq/
cd /home/Research_Archive/ProcessedData/Nanopore_seq/$study_name/fastq/

src_dir=/home/Research_Archive/RawData/Gene_Therapy/Nanopore/$study_name/fastq_pass
# cat $src_dir/*/*.fastq | ~/anaconda3/bin/qcat -b ./ --trim  -k PBC001 >demultiplex.stats.txt 2>&1  & 

for barcode in barcode07 barcode08 barcode09 barcode10 barcode11 barcode12; do
  echo $barcode
  cat $src_dir/$barcode/*fastq >$barcode.fastq &
done
wait

##rename samples
mkdir renamed
cd renamed

cp ../barcode07.fastq SC-AAV8-Stuffer-U7-hUsh2aEX13-46-01.fastq #1-Rep2-NewHelper.fastq #ID23.fastq
cp ../barcode08.fastq SC-AAV8-Stuffer-U7-hUsh2aEX13-48-01.fastq #2-Rep2-pHelper.fastq #SCA3_3.fastq
cp ../barcode09.fastq SC-AAV2-Stuffer-U7-hUsh2aEX13-46-01.fastq #3-Rep2-PADHelper.fastq #SCA3_6.fastq
cp ../barcode10.fastq SC-AAV2-Stuffer-U7-hUsh2aEX13-48-01.fastq #4-Rep2-Aldx80Helper.fastq #SCA3_12.fastq
cp ../barcode11.fastq SC-AAV9-SCA3-E10x2-CAGDel7-ID43-01.fastq #5-Rep3-pHelper.fastq #SS-AAV9-WT.CMV-Reelin-R3p6v2-V5-01.fastq
cp ../barcode12.fastq SC-AAV8-Stuffer-U7-hUsh2aEX13-46-02.fastq #6-Rep3-NewHelper.fastq
#cp ../barcode09.fastq Dual-Transfection.fastq #SS-AAV9-WT.CMV-Reelin-R3p6v2-V5-03-04.fastq

# no barcode11 ln ../barcode11.fastq 
# no barcode12 ln ../barcode12.fastq 

##2021-03-16 vi ~/.bash_profile # change PATH=$PATH:/drive2/wli/software2/biosoft/seq_ana/ggsashimi # :wq

##main settings:
study_name=2022-03-25_AAV_SEQ

samples=(SC-AAV8-Stuffer-U7-hUsh2aEX13-46-01 SC-AAV8-Stuffer-U7-hUsh2aEX13-48-01 SC-AAV2-Stuffer-U7-hUsh2aEX13-46-01 SC-AAV2-Stuffer-U7-hUsh2aEX13-48-01 SC-AAV9-SCA3-E10x2-CAGDel7-ID43-01 SC-AAV8-Stuffer-U7-hUsh2aEX13-46-02)
helpers=(pHelper-KanR pHelper-KanR pHelper-KanR pHelper-KanR pHelper-KanR pHelper-KanR)
RepCaps=(pAAV-Rep2Cap8 pAAV-Rep2Cap8 pAAV-RC2 pAAV-RC2 pAAV-Rep2Cap9-KanR pAAV-Rep2Cap8)
GOIs=(Stuffer-U7hUsh2aEX13-46_shift Stuffer-U7hUsh2aEX13-48_shift Stuffer-U7hUsh2aEX13-46_shift Stuffer-U7hUsh2aEX13-48_shift scAAV-SCA3-E10x2_CAGdel7 Stuffer-U7hUsh2aEX13-46_shift)
GOIfolders=(Ush2a Ush2a Ush2a Ush2a SCA3 Ush2a)
colors=(orange green red pink red3 black)
barcodes=(NBD07 NBD08 NBD09 NBD10 NBD11 NBD12)
aav_sample_ids=(`seq 1 ${#samples[*]}`)

all_AAV_samples=()
for i in ${aav_sample_ids[*]}; do
  indx=$(($i-1));
  sample=${samples[indx]};
  all_AAV_samples[${#all_AAV_samples[@]}]="$sample";
done
echo ${all_AAV_samples[*]}

# minimap2 folder: export PATH="$PATH:"`pwd`
# samtools folder: export PATH="$PATH:"`pwd`
# Use Build_minimap2_index.R(/drive2/wli/analyze/dep_seq/othProject/Nanopore_seq/) to build an index folder

GOI_fa_dir=/drive2/wli/analyze/summary/GeneTherapy/NanoporeSequencing/VirusGenomeSequences
index_folder=/home/Research_Archive/ProcessedData/Nanopore_seq/target_index
sample_info_f=/drive2/wli/analyze/summary/GeneTherapy/NanoporeSequencing/$study_name/$study_name-Samples.xlsx

root_dir=/drive2/bwang #wli/
pipeline_root=/drive2/wli/analyze/Pipeline/ #this is root directory of the pipeline/scripts (where the setGlobalVars.sh is located)
study_name=2022-03-25_AAV_SEQ #2021-03-10_AAV_SEQ_POOL1
mapping_data_root=/home/Research_Archive/ProcessedData/Nanopore_seq/$study_name/
project_root=/drive2/wli/analyze/dep_seq/othProject/Nanopore_seq/2022-03-25_AAV_SEQ/ #2021-03-10_AAV_SEQ_POOL1/ #project folder (all the input and output files related to this project)
MinIONQC_out_dir=$mapping_data_root/MinIONQC/$study_name/


##3, calculate total read number, and gzip fastq files
cd /home/Research_Archive/ProcessedData/Nanopore_seq/$study_name/fastq/renamed

out_readNum_f=readnum.txt
echo -e "Sample\tRaw read" >$out_readNum_f
p=0
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  sample1=${samples[indx]}
  rowNum=`wc -l $sample1.fastq | sed s/[[:blank:]].*//`
  rawReadNum=$(( $rowNum / 4 ))
  echo -e "$sample1\t$rawReadNum" >>$out_readNum_f
  gzip $sample1.fastq &
  p=$(($p+1)); echo $p; if [ $p -ge 8 ]; then  p=0 ;  wait; fi
done
wait


##4, mapping

##4.1 map to GOI-vector first
cd $mapping_data_root
mkdir minimap2
fq_dir=$mapping_data_root/fastq/renamed/ #the directory of the fastq files 


for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  fq1=$fq_dir/${samples[indx]}.fastq.gz
  out_dir=minimap2
  out_sam_f=$out_dir/${samples[indx]}.map2GOI.sam
  minimap2_index_f=$GOI_fa_dir/${GOIfolders[indx]}/${GOIs[indx]}.fa
  mkdir $out_dir
  minimap2 -ax map-ont -t 6 $minimap2_index_f $fq1 >$out_sam_f
done

##4.2, find unmapped reads and output fastq format
study_name=2022-03-25_AAV_SEQ #2021-03-10_AAV_SEQ_POOL1
cd $pipeline_root/RNAseq

p=0
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  perl 01.Cal_geno_junc_readsnum_frSAM.pl -n $study_name -s "${samples[indx]}"  -d "$mapping_data_root/minimap2/" -o "$project_root/01.ReadTable/map2GOI/" \
    -e ".map2GOI.sam" -t minimap2 -q 0 -m 30 -u 0 -U $mapping_data_root/fastq/${samples[indx]}.unmapped2GOI.fq &
  p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
done
wait

cd $pipeline_root/RNAseq/
perl 01a.cal_readnum.pl -o $project_root/01.ReadTable/map2GOI/


# ##5.1 map to all other GOI sequences (performed in same sequencing batch), 5.1 and 5.2 can be skipped if the GOI is the only GOI in the sequencing batch
cd $mapping_data_root
mkdir minimap2
fq_dir=$mapping_data_root/fastq/ #the directory of the fastq files 
 
 
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  fq1=$fq_dir/${samples[indx]}.unmapped2GOI.fq
  out_dir=minimap2
  out_sam_f=$out_dir/${samples[indx]}.map2otherGOI.sam
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  minimap2_index_f=$project_root/othGOI_fa/$index_name.othGOI.fa
  mkdir $out_dir
  minimap2 -ax map-ont -t 6 $minimap2_index_f $fq1 >$out_sam_f
done

# ##5.2, find unmapped reads and output fastq format
study_name=2022-03-25_AAV_SEQ #2021-03-10_AAV_SEQ_POOL1
cd $pipeline_root/RNAseq
 
p=0
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  perl 01.Cal_geno_junc_readsnum_frSAM.pl -n $study_name -s "${samples[indx]}"  -d "$mapping_data_root/minimap2/" -o "$project_root/01.ReadTable/map2othGOI/" \
    -e ".map2otherGOI.sam" -t minimap2 -q 0 -m 30 -u 0 -U $mapping_data_root/fastq/${samples[indx]}.unmapped2GOI.fq &
  p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
done
wait
# 
cd $pipeline_root/RNAseq/
perl 01a.cal_readnum.pl -o $project_root/01.ReadTable/map2othGOI/


##6, map to the rest of vectors and genomes
cd $mapping_data_root
mkdir minimap2
fq_dir=$mapping_data_root/fastq/ #the directory of the fastq files 


for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  fq1=$fq_dir/${samples[indx]}.unmapped2GOI.fq
  out_dir=minimap2
  out_sam_f=$out_dir/${samples[indx]}.map2oths.sam
  echo out_sam_f
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  minimap2_index_f=$index_folder/$index_name/$index_name.mmi
  mkdir $out_dir
  minimap2 -ax map-ont -t 6 $minimap2_index_f $fq1 >$out_sam_f
done

##7, merge bam files
cd $mapping_data_root
out_dir=minimap2
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  sam1=$out_dir/${samples[indx]}.map2GOI
  samtools view -b -o $sam1.bam $sam1.sam &
  sam2=$out_dir/${samples[indx]}.map2oths
  samtools view -b -o $sam2.bam $sam2.sam &
  wait
  samtools merge -f -@ 6 $out_dir/${samples[indx]}.bam $sam1.bam $sam2.bam
done
####2021/02/17: AAV-GTX-temp.bam had an issue to  merge !!! :Aborted (core dumped)
##resolved after installed a new version of samtools: /home/wli/soft/biosoft/seq_map/samtools/samtools-1.11/install/bin/

##3/14/2021

##8, find uniquely mapped reads with good quality from sam or bam file calculate read mapping to different plasmids or chromosomes
study_name=2022-03-25_AAV_SEQ #2021-03-10_AAV_SEQ_POOL1
cd $project_root/..
p=0
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  sample=${samples[indx]}
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  GOI_info_f=$index_folder/$index_name/$index_name.feature.txt

  perl 02.count_AAV_read_fromBam.pl -n $study_name -s "$sample"  -b "$mapping_data_root/minimap2/$sample.bam" \
    -p "$study_name/02_Bam2ReadCounts/"  \
    -t minimap2 -u 10 -m 30  -v "$GOI_info_f" \
    -S "$mapping_data_root/minimap2/$sample.highQuality.sam" &
  p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
done
wait

##combine mapping summary
cd $pipeline_root/RNAseq/
Rscript 02a.combinePE.readstats.R -in_dir "$project_root/02_Bam2ReadCounts/" -out_f $project_root/02_Bam2ReadCounts/All.read.counts.txt \
  -inf_pattern .read.counts.txt.log 




##9, create bam index (and bigwig) files for visualization:
study_name=2022-03-25_AAV_SEQ #2021-03-10_AAV_SEQ_POOL1
geno=hg19
PROG_DIR=$pipeline_root/Software/UCSC
fext='.highQuality.sam' #string after sample name in input sam file name

function sam2bw {
 #totalReadNum=`samtools view -c -q 10 $sample$fext`
 #echo "# $sample for file $sample$fext, TotalReadNum=$totalReadNum "
 echo "# $sample, convert to bam "
 samtools view -S -b -T $genofasta -o $sample.highQuality.bam  $sample$fext 
 #echo "# $sample filter original bam file (keep uniquely mapped reads only) "
 #samtools view -b -q 10 -o $sample.filtered.bam $sample.bam
 
 echo "# $sample, sort bam file"
 samtools sort -o $sample.sorted.bam $sample.highQuality.bam
 
 echo "# $sample, create index file for bam"
 samtools index $sample.sorted.bam
 #rm $sample.bam
 
 echo "# $sample, convert to bedgraph"
 genomeCoverageBed -bg -split -ibam $sample.sorted.bam -g $genofasta.sizes > $sample.bedgraph
 echo "# $sample, convert to bw"
 $PROG_DIR/bedGraphToBigWig $sample.bedgraph $genofasta.sizes $sample.bw 
}


p=0
cd $mapping_data_root/minimap2/
for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  sample=${samples[indx]}
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  genofasta=$index_folder/$index_name/$index_name.fa
  sam2bw  &
  p=$(($p+1)); echo $p; if [ $p -ge 6 ]; then  p=0 ;  wait; fi
done
wait

cd $mapping_data_root/minimap2/
# rm *.bedgraph
# rm *.sam
# rm *.map2GOI.bam
# rm *.map2oths.bam
# rm *.highQuality.bam


##10, analyze read mapping results
cd $project_root/..

for i in `seq 1 ${#samples[*]}`; do
  indx=$(($i-1));
  sample=${samples[indx]}
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  ref_sizes_f=$index_folder/$index_name/$index_name.fa.sizes
  GOI_info_f=$index_folder/$index_name/$index_name.feature.txt

  Rscript 03.ana_readMapping.R -study_name $study_name -samples "$sample" \
    -read_in_dir $study_name/02_Bam2ReadCounts/  -out_dir $study_name/03_ana_readMap/$sample. \
    -ref_sizes_f $ref_sizes_f -GOI_info_f $GOI_info_f
done



##11, draw read coverage in different chromosomes:
study_name=2022-03-25_AAV_SEQ #2021-03-10_AAV_SEQ_POOL1
cd $project_root/..

Rscript 04GeneReadCovPlot.R -study_name $study_name -sample_info_f $sample_info_f \
   -mapping_data_root $mapping_data_root -index_folder $index_folder


##12, call peaks of human genome mapping

##or 12.b use a custom R script to cluster reads
cd $project_root/..
for i in ${aav_sample_ids[*]}; do
  indx=$(($i-1));
  sample=${samples[indx]}
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  ref_sizes_f=$index_folder/$index_name/$index_name.fa.sizes
  Rscript 07clusterReads.R -study_name $study_name -samples "$sample" \
    -read_in_dir $study_name/02_Bam2ReadCounts/  -out_dir $study_name/07_readCluster/$sample. \
    -ref_sizes_f $ref_sizes_f 
done

all_AAV_samples=()
for i in ${aav_sample_ids[*]}; do
  indx=$(($i-1));
  sample=${samples[indx]};
  all_AAV_samples[${#all_AAV_samples[@]}]="$sample";
done
echo ${all_AAV_samples[*]}

##cluster all samples:
  i=${aav_sample_ids[0]}
  indx=$(($i-1));
  sample=${samples[indx]}
  index_name=${helpers[indx]}.${RepCaps[indx]}.${GOIs[indx]}
  ref_sizes_f=$index_folder/$index_name/$index_name.fa.sizes
  Rscript 07clusterReads.R -study_name $study_name -samples "${all_AAV_samples[*]}" \
    -read_in_dir $study_name/02_Bam2ReadCounts/  -out_dir $study_name/07_readCluster/AllSamples. \
    -ref_sizes_f $ref_sizes_f 

##12.c annotate peaks or clusters (bed format file)
cd $project_root/..
Rscript 08annotate_cluster.R -study_name $study_name -samples "AllSamples ${all_AAV_samples[*]}"




##13, generate sequencing report
study_name=2022-03-25_AAV_SEQ #2021-03-10_AAV_SEQ_POOL1

cd $project_root/..
/bin/Rscript 05prepare_report.R -study_name $study_name -samples "${samples[*]}" -barcodes "${barcodes[*]}" -project_root $project_root \
  -mapping_data_root $mapping_data_root -MinIONQC_out_dir "$MinIONQC_out_dir" -sample_info_f $sample_info_f \
  -AAV_samples "AllSamples ${all_AAV_samples[*]}"

cd $project_root/05.Report
for sample in All ${samples[*]}; do
  echo $sample
  /bin/Rscript -e "out_sample='$sample'; rmarkdown::render('NanoporeSeq_report.Rmd', output_format='html_document', output_file = 'NanoporeSeq_report.$sample.html') " 
  /bin/Rscript -e "out_sample='$sample'; rmarkdown::render('NanoporeSeq_report.Rmd', output_format='word_document', output_file = 'NanoporeSeq_report.$sample.docx') " 
done
