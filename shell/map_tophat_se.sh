#!/bin/bash

# Directory variables
HOME_PATH=/media/raid/tmp/kontaki
FASTQ_PATH=$HOME_PATH/fastq
TOPHAT_OUT=$HOME_PATH/tophat_out

# Command tool variables
SAMTOOLS_COMMAND=/opt/ngstools/samtools-0.1.19/samtools
TOPHAT_COMMAND=/opt/ngstools/tophat/tophat2
BOWTIE2_COMMAND=/opt/ngstools/bowtie2/bowtie2
BEDTOOLS_COMMAND=/opt/ngstools/bedtools2/bin/bedtools
PICARD_COMMAND="java -jar /opt/ngstools/picard/AddOrReplaceReadGroups.jar"

# Index
BOWTIE2_INDEX=/media/raid/resources/igenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome
TRANSCRIPTOME=/media/raid/resources/transcriptomes/mm10/Bowtie2Index/genes

# Annotation
GTF=/media/raid/resources/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf

# Organism
ORG=mm10

# Number of cores to use
CORES=32

# Execute tophat2
for FILE in $FASTQ_PATH/*.fastq
do
    SAMPLE=`basename $FILE | sed s/\.fastq//`
    echo "==================== Mapping with tophat2 for $SAMPLE..."
    mkdir -p $TOPHAT_OUT/$SAMPLE
    $TOPHAT_COMMAND --no-coverage-search --num-threads $CORES --GTF $GTF --transcriptome-index $TRANSCRIPTOME --output-dir $TOPHAT_OUT/$SAMPLE $BOWTIE2_INDEX $FASTQ_PATH/$SAMPLE".fastq"
    echo " "
    echo "==================== Name sorting unmapped reads to prepare for second round mapping for $SAMPLE..."
    $SAMTOOLS_COMMAND sort -n $TOPHAT_OUT/$SAMPLE/unmapped.bam $TOPHAT_OUT/$SAMPLE/unmapped_tmp
    echo "==================== Converting unmapped reads to fastq for $SAMPLE..."
    $BEDTOOLS_COMMAND bamtofastq -i $TOPHAT_OUT/$SAMPLE/unmapped_tmp.bam -fq $TOPHAT_OUT/$SAMPLE/unmapped.fastq 2> /dev/null
    echo "==================== Trying now to map unmapped reads with bowtie2 for $SAMPLE..."
    $BOWTIE2_COMMAND --local --very-sensitive-local -p $CORES -x $BOWTIE2_INDEX -U $TOPHAT_OUT/$SAMPLE/unmapped.fastq | $SAMTOOLS_COMMAND view -bhS -o $TOPHAT_OUT/$SAMPLE/unmapped_remap.bam -
    echo "==================== Merging all reads for $SAMPLE..."
    $SAMTOOLS_COMMAND merge -f $TOPHAT_OUT/$SAMPLE/tmp.bam $TOPHAT_OUT/$SAMPLE/accepted_hits.bam $TOPHAT_OUT/$SAMPLE/unmapped_remap.bam
    echo "==================== Coordinate sorting all reads for $SAMPLE..."
    $SAMTOOLS_COMMAND sort $TOPHAT_OUT/$SAMPLE/tmp.fix.nsort.bam $TOPHAT_OUT/$SAMPLE/merged_tmp
    $PICARD_COMMAND I=$TOPHAT_OUT/$SAMPLE/merged_tmp.bam O=$TOPHAT_OUT/$SAMPLE/merged.bam LB=$ORG PL=IONTORRENT PM=Ion Torrent Proton SM=merged VALIDATION_STRINGENCY=SILENT
    echo "==================== Indexing all merged reads for $SAMPLE..."
    $SAMTOOLS_COMMAND index $TOPHAT_OUT/$SAMPLE/merged.bam
    echo "==================== Exctracting final aligned reads for $SAMPLE..."
    $SAMTOOLS_COMMAND view -bh -F4 -o $TOPHAT_OUT/$SAMPLE/aligned.bam $TOPHAT_OUT/$SAMPLE/merged.bam 
    echo "==================== Indexing final aligned reads for $SAMPLE..."
    $SAMTOOLS_COMMAND index $TOPHAT_OUT/$SAMPLE/aligned.bam
    echo "==================== Removing intermediate garbage for $SAMPLE..."
    rm $TOPHAT_OUT/$SAMPLE/merged_tmp.bam $TOPHAT_OUT/$SAMPLE/unmapped_tmp.bam $TOPHAT_OUT/$SAMPLE/tmp.bam $TOPHAT_OUT/$SAMPLE/tmp.nsort.bam $TOPHAT_OUT/$SAMPLE/tmp.fix.nsort.bam $TOPHAT_OUT/$SAMPLE/unmapped*R*.fastq $TOPHAT_OUT/$SAMPLE/unmapped_remap.bam
    echo " "
done
