#!/bin/bash
#PBS -q standard
#PBS -l select=1:ncpus=7:mem=42gb:pcmem=6gb
#PBS -l pvmem=168gb
#PBS -N bATAC
#PBS -W group_list=cusanovichlab
#PBS -l place=pack:shared
#PBS -l walltime=24:00:00
### cput ("total cputime") = walltime * ncpus
##PBS -l cput=16:00:00 

date
module load bowtie2/2.2.9
module load samtools/1.4
module load bedtools/2.26.0
module load R/3.6.1
#############################
##EDIT THESE AS NEEDED:
#############################
#genome=mm
genome=hs
#genome_dir=/rsgrps/cusanovichlab/references/bowtie/mm10/mm10
genome_dir=/rsgrps/cusanovichlab/references/bowtie/hg19/hg19
#tss_bed=/rsgrps/cusanovichlab/references/annotations/mm10/mm10.tss.liftedover.bed.gz
tss_bed=/rsgrps/cusanovichlab/references/annotations/hg19/hg19.tss.bed.gz
fastq1=fastq1
fastq2=fastq2
base=base
OUTDIR=OUTDIR
#############################

mkdir $OUTDIR/fastqs_trimmed
mkdir $OUTDIR/bams
mkdir $OUTDIR/beds
mkdir $OUTDIR/reports
mkdir $OUTDIR/macs

echo running
# do not forget to change Trimmomatic path and adapter path, download picard.jar, install macs2, Path for R scripts
java -jar /rsgrps/cusanovichlab/sbin/Trimmomatic-0.36/trimmomatic-0.36.jar PE $OUTDIR/fastqs/$fastq1 $OUTDIR/fastqs/$fastq2 $OUTDIR/fastqs_trimmed/${base}_R1.fastq.paired.trimmed.gz $OUTDIR/fastqs_trimmed/${base}_R1.fastq.unpaired.trimmed.gz $OUTDIR/fastqs_trimmed/${base}_R2.fastq.paired.trimmed.gz $OUTDIR/fastqs_trimmed/${base}_R2.fastq.unpaired.trimmed.gz ILLUMINACLIP:/home/u23/haozhang1/sciatac_pipeline/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10:1:true TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:20 2> $OUTDIR/reports/${base}_trimmomatic.log;
bowtie2 -p 8 -X 2000 -3 1 -x $genome_dir -1 $OUTDIR/fastqs_trimmed/${base}_R1.fastq.paired.trimmed.gz -2 $OUTDIR/fastqs_trimmed/${base}_R2.fastq.paired.trimmed.gz 2> $OUTDIR/reports/${base}.bowtie2.log | samtools view -bS - > $OUTDIR/bams/${base}.bam;
samtools view -h -f3 -F12 -q10 $OUTDIR/bams/${base}.bam | grep -v '[0-9]'$'\t'chrM | grep -v '[0-9]'$'\t'chrUn | grep -v '[0-9]'$'\t'chr.*_random$'\t' | samtools view -Su - | samtools sort -@ 8 - -o $OUTDIR/bams/${base}.q10.sort.bam;
# if get error for "Exception in thread "main" java.lang.OutOfMemoryError: Java heap space", try to increase Xmx to 4 or 8 G
java -Xmx1G -jar /rsgrps/cusanovichlab/sbin/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$OUTDIR/bams/${base}.q10.sort.bam O=$OUTDIR/bams/${base}.nodups.bam M=$OUTDIR/reports/${base}.dedup.log.txt;
java -Xmx1G -jar /rsgrps/cusanovichlab/sbin/picard.jar CollectInsertSizeMetrics INPUT=$OUTDIR/bams/${base}.nodups.bam OUTPUT=$OUTDIR/reports/${base}.insertsizes.txt H=$OUTDIR/reports/${base}.insertsizes.pdf;
bedtools bamtobed -i $OUTDIR/bams/${base}.nodups.bam > $OUTDIR/beds/${base}.nodups.bed;
macs2 callpeak -t $OUTDIR/beds/${base}.nodups.bed -f BED -g $genome --nomodel --keep-dup all --extsize 200 --shift -100 --call-summits -n $OUTDIR/macs/${base}.nodups;
sort -k 8gr,8gr $OUTDIR/macs/${base}.nodups_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | sort -k1,1 -k2,2n -k3,3n | gzip -c > $OUTDIR/macs/${base}.nodups.narrowPeak.gz;
rm $OUTDIR/macs/${base}.nodups_peaks.narrowPeak;
zcat $tss_bed \
		| awk -vOFS="\\t" -vEXT=2000 '{{ print $1,$2-EXT,$3+EXT,$4,0,$6}}' \
		| bedtools coverage -sorted -d -a stdin -b $OUTDIR/beds/${base}.nodups.bed \
		| awk '{{ if ($8 > 0) print $0 }}' \
		| gzip > $OUTDIR/reports/${base}.nodups.tssfile.temp.gz;
Rscript /rsgrps/cusanovichlab/sbin/bulk_atac/aggregate_per_base_tss_region_counts.R $OUTDIR/reports/${base}.nodups.tssfile.temp.gz $OUTDIR/reports/${base}.nodups.tsscoverage.tsv;
rm $OUTDIR/reports/${base}.nodups.tssfile.temp.gz;
zcat $OUTDIR/macs/${base}.nodups.narrowPeak.gz | bedtools merge -i stdin > $OUTDIR/macs/${base}.nodups.narrowPeak_merged.bed
Rscript /rsgrps/cusanovichlab/sbin/bulk_atac/tss_enrichment_plot.R --bowtie_log=$OUTDIR/reports/${base}.bowtie2.log --picard_log=$OUTDIR/reports/${base}.dedup.log.txt --picard_inserts_file=$OUTDIR/reports/${base}.insertsizes.txt --bam_file=$OUTDIR/bams/${base}.bam --dedup_file=$OUTDIR/bams/${base}.nodups.bam --bed_file=$OUTDIR/beds/${base}.nodups.bed --peak_file=$OUTDIR/macs/${base}.nodups.narrowPeak_merged.bed --per_base_tss_region_coverage_file=$OUTDIR/reports/${base}.nodups.tsscoverage.tsv --plot=$OUTDIR/reports/${base}.nodups.tss_enrichment.pdf;
date
