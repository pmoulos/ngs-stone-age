#!/bin/bash

# Directory variables
HOME_PATH=/media/raid/tmp/karolos_xs
FASTQ_PATH=$HOME_PATH/fastq
TOPHAT_OUT=$HOME_PATH/tophat_out

# Command tool variables
SAMTOOLS_COMMAND=/opt/ngstools/samtools-0.1.19/samtools
#TOPHAT_COMMAND=/opt/ngstools/tophat/tophat2
BOWTIE2_COMMAND=/opt/ngstools/bowtie2/bowtie2
BEDTOOLS_COMMAND=/opt/ngstools/bedtools2/bin/bedtools
PICARD_COMMAND="java -jar /opt/ngstools/picard/AddOrReplaceReadGroups.jar"

# Index
BOWTIE2_INDEX=$HOME_PATH/reference_genome/bmori2
#TRANSCRIPTOME=$HOME_PATH/reference_genome/genes

# Annotation
#GTF=$HOME_PATH/reference_genome/genes.gtf

# Organism
ORG=bmori2

# Number of cores to use
CORES=32

# Execute tophat2
for FILE in $FASTQ_PATH/*_R1.fastq
do
    SAMPLE=`basename $FILE | sed s/\_R1.fastq//`
    echo "==================== Mapping with tophat2 for $SAMPLE..."
    mkdir -p $TOPHAT_OUT/$SAMPLE
    $TOPHAT_COMMAND --no-coverage-search --num-threads $CORES --GTF $GTF --transcriptome-index $TRANSCRIPTOME --read-gap-length 3 --read-edit-dist 3  --mate-inner-dist 100 --mate-std-dev 120 --output-dir $TOPHAT_OUT/$SAMPLE $BOWTIE2_INDEX $FASTQ_PATH/$SAMPLE"_R1.fastq" $FASTQ_PATH/$SAMPLE"_R2.fastq"
    echo " "
    echo "==================== Name sorting unmapped reads to prepare for second round mapping for $SAMPLE..."
    $SAMTOOLS_COMMAND sort -n $TOPHAT_OUT/$SAMPLE/unmapped.bam $TOPHAT_OUT/$SAMPLE/unmapped_tmp
    echo "==================== Converting unmapped reads to fastq for $SAMPLE..."
    $BEDTOOLS_COMMAND bamtofastq -i $TOPHAT_OUT/$SAMPLE/unmapped_tmp.bam -fq $TOPHAT_OUT/$SAMPLE/unmapped_R1.fastq -fq2 $TOPHAT_OUT/$SAMPLE/unmapped_R2.fastq 2> /dev/null
    echo "==================== Trying now to map unmapped reads with bowtie2 for $SAMPLE..."
    $BOWTIE2_COMMAND --local --very-sensitive-local --maxins 1000 --dovetail -p $CORES -x $BOWTIE2_INDEX -1 $TOPHAT_OUT/$SAMPLE/unmapped_R1.fastq -2 $TOPHAT_OUT/$SAMPLE/unmapped_R2.fastq | $SAMTOOLS_COMMAND view -bhS -o $TOPHAT_OUT/$SAMPLE/unmapped_remap.bam -
    echo "==================== Merging all reads for $SAMPLE..."
    $SAMTOOLS_COMMAND merge -f $TOPHAT_OUT/$SAMPLE/tmp.bam $TOPHAT_OUT/$SAMPLE/accepted_hits.bam $TOPHAT_OUT/$SAMPLE/unmapped_remap.bam
    echo "==================== Name sorting all reads bam files for $SAMPLE..."
    $SAMTOOLS_COMMAND sort -n $TOPHAT_OUT/$SAMPLE/tmp.bam $TOPHAT_OUT/$SAMPLE/tmp.nsort
    echo "==================== Mate fixing all reads for $SAMPLE..."
    $SAMTOOLS_COMMAND fixmate $TOPHAT_OUT/$SAMPLE/tmp.nsort.bam $TOPHAT_OUT/$SAMPLE/tmp.fix.nsort.bam
    echo "==================== Coordinate sorting all mate fixed reads for $SAMPLE..."
    $SAMTOOLS_COMMAND sort $TOPHAT_OUT/$SAMPLE/tmp.fix.nsort.bam $TOPHAT_OUT/$SAMPLE/merged_tmp
    $PICARD_COMMAND I=$TOPHAT_OUT/$SAMPLE/merged_tmp.bam O=$TOPHAT_OUT/$SAMPLE/merged.bam LB=$ORG PL=illumina PU=illumina_merged SM=merged VALIDATION_STRINGENCY=SILENT
    echo "==================== Indexing all merged mate fixed reads for $SAMPLE..."
    $SAMTOOLS_COMMAND index $TOPHAT_OUT/$SAMPLE/merged.bam
    echo "==================== Exctracting final aligned reads for $SAMPLE..."
    $SAMTOOLS_COMMAND view -bh -F4 -o $TOPHAT_OUT/$SAMPLE/aligned.bam $TOPHAT_OUT/$SAMPLE/merged.bam 
    echo "==================== Indexing final aligned reads for $SAMPLE..."
    $SAMTOOLS_COMMAND index $TOPHAT_OUT/$SAMPLE/aligned.bam
    echo "==================== Removing intermediate garbage for $SAMPLE..."
    rm $TOPHAT_OUT/$SAMPLE/merged_tmp.bam $TOPHAT_OUT/$SAMPLE/unmapped_tmp.bam $TOPHAT_OUT/$SAMPLE/tmp.bam $TOPHAT_OUT/$SAMPLE/tmp.nsort.bam $TOPHAT_OUT/$SAMPLE/tmp.fix.nsort.bam $TOPHAT_OUT/$SAMPLE/unmapped*R*.fastq $TOPHAT_OUT/$SAMPLE/unmapped_remap.bam
    echo " "
done
