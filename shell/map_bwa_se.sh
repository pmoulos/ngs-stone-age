# Paths and vars
BWA_PATH=/opt/ngstools/bwa
FASTQ_PATH=/media/raid/tmp/yota_chipseq_hu/fastq
CLEAN_PATH=$FASTQ_PATH
SAMTOOLS=/opt/ngstools/samtools-0.1.19
THREADS=32
BWA_INDEX=/media/raid/resources/igenomes/Mus_musculus/UCSC/mm9/Sequence/BWAIndex/version0.6.0/genome.fa
BEDTOOLS_COMMAND="/opt/ngstools/bedtools2/bin/bedtools bamtobed"

printf "%s\t%s\t%s\t%s\t%s\t%s\n" "name" "total" "mapped" "qc_filtered" "unique" "clean" > $FASTQ_PATH/bwa_report.txt

echo "================================================================================"
echo "----> Mapping with bwa"
echo "================================================================================"
for FILE in $FASTQ_PATH/*.fastq.gz
do
    BASE=`basename $FILE | sed s/\.fastq\.gz//`
    echo "========== Processing $BASE..."   
    #echo "     ===== Unzipping..."
    #pigz -d -p $THREADS $FILE
    echo "     ===== Running bwa with defaults..."
    $BWA_PATH/bwa aln -t $THREADS $BWA_INDEX $FASTQ_PATH/$BASE".fastq.gz" > $FASTQ_PATH/$BASE".sai"
    $BWA_PATH/bwa samse $BWA_INDEX $FASTQ_PATH/$BASE".sai" $FASTQ_PATH/$BASE".fastq.gz" | $SAMTOOLS/samtools view -bS -o $FASTQ_PATH/$BASE".bam" -
    echo "     ===== Processing the output..."
    MAPPED=$FASTQ_PATH/$BASE".bam"
    printf "%s\t" $BASE >> $FASTQ_PATH/bwa_report.txt
    echo "     ===== Counting total reads..."
    printf "%d\t" `$SAMTOOLS/samtools view -c $MAPPED` >> $FASTQ_PATH/bwa_report.txt
    echo "     ===== Counting mapped reads..."
    printf "%d\t" `$SAMTOOLS/samtools view -c -F 4 $MAPPED` >> $FASTQ_PATH/bwa_report.txt
    echo "     ===== Counting mapped and high-quality reads..."
    printf "%d\t" `$SAMTOOLS/samtools view -c -F 4 -q 10 $MAPPED` >> $FASTQ_PATH/bwa_report.txt
    echo "     ===== Converting BAM file $MAPPED to BED file..."
    $BEDTOOLS_COMMAND -i $MAPPED > $FASTQ_PATH/$BASE".tmp.bed"
    echo "     ===== Cleaning the BED file..."
    printf "%d\t" `sort -k1,1 -k2g,2 -u $FASTQ_PATH/$BASE".tmp.bed" | wc -l` >> $FASTQ_PATH/bwa_report.txt
    echo "     ===== Creating the final BED file..."
    grep -vE 'chrM|rand' $FASTQ_PATH/$BASE".tmp.bed" | sort -k1,1 -k2g,2 -u | pigz -p $THREADS > $FASTQ_PATH/$BASE".bed.gz"
    rm $FASTQ_PATH/$BASE".tmp.bed" $FASTQ_PATH/$BASE".sai"
    printf "%d\n" `zcat $FASTQ_PATH/$BASE".bed.gz" | wc -l` >> $FASTQ_PATH/bwa_report.txt
    #echo "     ===== Zipping the original fastq..."
    #pigz -p $THREADS $FASTQ_PATH/$BASE".fq"
done

