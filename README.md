# Bulk ATAC-seq pipeline
## This pipeline aims to preprocess bulk ATAC-seq data, call peaks, and create a count matrix for differential analysis.
### Linux packages needed to install.
```
bowtie2/2.4.1
samtools/1.10
Trimmomatic-0.36
picard.jar
bedtools
macs2
```
### R packages needed to install.
```
dplyr
reader
```
### Before starting to run the pipeline, define the PATH or value for the variables below.
```
OUTDIR=/path/to/your/output/dir
Rscript=/path/to/your/output/dir
# Genome size parameters for MACS2, see https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html for details
genome=mm
genome_dir=/dir/to/bowtie/genome/reference
tss_bed=/bed/mm10.v23.tss.bed.gz
fastq1=/dir/to/fastq/read1
fastq2=/dir/to/fastq2/read2
base=file_name_prefix
```
### Create the following folders.
```
mkdir $OUTDIR/fastqs_trimmed
mkdir $OUTDIR/bams
mkdir $OUTDIR/beds
mkdir $OUTDIR/reports
mkdir $OUTDIR/macs
```
### Date preprocessing
#### Remove sequencing adaptor and low-quality reads.
```
java -jar Trimmomatic-0.36/trimmomatic-0.36.jar \
PE $fastq1 $fastq2 \
$OUTDIR/fastqs_trimmed/${base}_R1.fastq.paired.trimmed.gz $OUTDIR/fastqs_trimmed/${base}_R1.fastq.unpaired.trimmed.gz \
$OUTDIR/fastqs_trimmed/${base}_R2.fastq.paired.trimmed.gz $OUTDIR/fastqs_trimmed/${base}_R2.fastq.unpaired.trimmed.gz \
ILLUMINACLIP:Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10:1:true \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:20 2> $OUTDIR/reports/${base}_trimmomatic.log;
```
#### Mapping reads and then filtering the mapped reads to only include the reads that were properly paired and confidently mapped to the assembled nuclear chromosomes (MAPQ ≥ 10).
```
bowtie2 -p 8 -X 2000 -3 1 -x $genome_dir \
-1 $OUTDIR/fastqs_trimmed/${base}_R1.fastq.paired.trimmed.gz \
-2 $OUTDIR/fastqs_trimmed/${base}_R2.fastq.paired.trimmed.gz \
2> $OUTDIR/reports/${base}.bowtie2.log | samtools view -bS - > $OUTDIR/bams/${base}.bam;

samtools view -h -f3 -F12 -q10 $OUTDIR/bams/${base}.bam \
| grep -v '[0-9]'$'\t'chrM | grep -v '[0-9]'$'\t'chrUn | grep -v '[0-9]'$'\t'chr.*_random$'\t' \
| samtools view -Su - | samtools sort -@ 8 - -o $OUTDIR/bams/${base}.q10.sort.bam;
```
#### Deduplication and generate fragment size matrics.
```
java -Xmx1G -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$OUTDIR/bams/${base}.q10.sort.bam O=$OUTDIR/bams/${base}.nodups.bam M=$OUTDIR/reports/${base}.dedup.log.txt;

java -Xmx1G -jar picard.jar CollectInsertSizeMetrics INPUT=$OUTDIR/bams/${base}.nodups.bam OUTPUT=$OUTDIR/reports/${base}.insertsizes.txt H=$OUTDIR/reports/${base}.insertsizes.pdf;
```
#### Call peaks with MACS2.
```
bedtools bamtobed -i $OUTDIR/bams/${base}.nodups.bam > $OUTDIR/beds/${base}.nodups.bed;

macs2 callpeak -t $OUTDIR/beds/${base}.nodups.bed -f BED -g $genome --nomodel --keep-dup all --extsize 200 --shift -100 --call-summits -n $OUTDIR/macs/${base}.nodups;

sort -k 8gr,8gr $OUTDIR/macs/${base}.nodups_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | sort -k1,1 -k2,2n -k3,3n | gzip -c > $OUTDIR/macs/${base}.nodups.narrowPeak.gz;

rm $OUTDIR/macs/${base}.nodups_peaks.narrowPeak;

zcat $OUTDIR/macs/${base}.nodups.narrowPeak.gz | bedtools merge -i stdin > $OUTDIR/macs/${base}.nodups.narrowPeak_merged.bed
```
#### Evaluate per base read counts around TSSs
```
zcat $tss_bed \
| awk -vOFS="\\t" -vEXT=2000 '{{ print $1,$2-EXT,$3+EXT,$4,0,$6}}' \
| bedtools coverage -sorted -d -a stdin -b $OUTDIR/beds/${base}.nodups.bed \
| awk '{{ if ($8 > 0) print $0 }}' \
| gzip > $OUTDIR/reports/${base}.nodups.tssfile.temp.gz;

$Rscript /bin/aggregate_per_base_tss_region_counts.R $OUTDIR/reports/${base}.nodups.tssfile.temp.gz $OUTDIR/reports/${base}.nodups.tsscoverage.tsv;

rm $OUTDIR/reports/${base}.nodups.tssfile.temp.gz;
```
#### Generate QC report
```
$Rscript /bin/tss_enrichment_plot.R --bowtie_log=$OUTDIR/reports/${base}.bowtie2.log \
--picard_log=$OUTDIR/reports/${base}.dedup.log.txt --picard_inserts_file=$OUTDIR/reports/${base}.insertsizes.txt \
--bam_file=$OUTDIR/bams/${base}.bam --dedup_file=$OUTDIR/bams/${base}.nodups.bam \
--bed_file=$OUTDIR/beds/${base}.nodups.bed \
--peak_file=$OUTDIR/macs/${base}.nodups.narrowPeak_merged.bed \
--per_base_tss_region_coverage_file=$OUTDIR/reports/${base}.nodups.tsscoverage.tsv \
--plot=$OUTDIR/reports/${base}.nodups.tss_enrichment.pdf;
```
### Call peaks on aggregated replicates and summarize reads for each peak across all sample
#### Set variables
```
mkdir $OUTDIR/deeptools
OUTDIR=$OUTDIR/deeptools
bams=$OUTDIR/bams
sample=glo2ko
mm10_blacklist_sort=/bed/mm10-blacklist.v2.sorted.bed
```
#### Make directory
```
mkdir $OUTDIR/bams
mkdir $OUTDIR/peak_ref
mkdir $OUTDIR/beds
mkdir $OUTDIR/macs
```
#### Merge replicates within each condition
```
for file in $bams/*1.nodups.bam
do
base=$(basename ${file} 1.nodups.bam)

echo "Merging $base for ${base}1.nodups.bam ${base}2.nodups.bam and ${base}3.nodups.bam..."
samtools merge $OUTDIR/bams/${base}.nodups.bam $bams/${base}1.nodups.bam $bams/${base}2.nodups.bam $bams/${base}3.nodups.bam
done
```
#### Call peaks on aggregated replicates
```
echo "Index bams..."
for file in $bams/*.nodups.bam; do samtools index $file; done

echo "Call peaks..."
for file in $OUTDIR/bams/*.nodups.bam
do
base=$(basename ${file} .nodups.bam)
samtools index $OUTDIR/bams/${base}.nodups.bam
echo "Converting ${base} bam to bed......"
bedtools bamtobed -i $OUTDIR/bams/${base}.nodups.bam > $OUTDIR/beds/${base}.nodups.bed
tot=$(cat $OUTDIR/beds/${base}.nodups.bed |wc -l)
echo "${base} total reads: ${tot}"
echo "Calling peaks on ${base} with ${genome} genome......"
macs2 callpeak -t $OUTDIR/beds/${base}.nodups.bed -f BED -g $genome --nomodel --keep-dup all --extsize 200 --shift -100 --call-summits -n $OUTDIR/macs/${base}.nodups;
sort -k 8gr,8gr $OUTDIR/macs/${base}.nodups_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' |\
sort -k1,1 -k2,2n -k3,3n | gzip -c > $OUTDIR/macs/${base}.nodups.narrowPeak.gz;
zcat $OUTDIR/macs/${base}.nodups.narrowPeak.gz | bedtools merge -i stdin > $OUTDIR/macs/${base}.nodups.narrowPeak_merged.bed
echo "Removing blacklist for ${base}......"
bedtools intersect -v -a $OUTDIR/macs/${base}.nodups.narrowPeak_merged.bed -b $mm10_blacklist_sort > $OUTDIR/peak_ref/${base}.whitelist.bed
done
```
#### Merge condition-specific peaks
```
cat $OUTDIR/peak_ref/*.whitelist.bed |\
sort -k1,1 -k2,2n -k3,3n | bedtools merge -i stdin >$OUTDIR/peak_ref/${sample}.merged_whitelist.bed
```
#### Make count matrix
```
peaks=$OUTDIR/peak_ref/${sample}.merged_whitelist.bed

echo "Making read count matrix......"
multiBamSummary BED-file --BED $peaks \
--bamfiles $bams/*.nodups.bam \
-o $OUTDIR/${sample}.npz \
--outRawCounts $OUTDIR/${sample}_rawcounts.txt -p max/2
```
#### Run differential analysis using edgeR_glmQL.R in the bin folder.
