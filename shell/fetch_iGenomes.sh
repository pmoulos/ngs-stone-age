#!/bin/bash

# We must be under sudo
if [ `id -u` -ne 0 ]; 
then
        echo "You must execute this script as root!"
        exit 1
fi

CURR_DIR=`pwd`
IGENOMES_HOME=/opt/iGenomes

if [ -d $IGENOMES_HOME ]
then
	echo "$IGENOMES_HOME has already been created..."
else
	mkdir -p $IGENOMES_HOME
fi

# Go there
cd $IGENOMES_HOME

# Fetch for human...
echo "==================== Fetching human iGenome... ===================="
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg18/Homo_sapiens_UCSC_hg18.tar.gz
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz

# Fetch for mouse...
echo "==================== Fetching mouse iGenome... ===================="
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm9/Mus_musculus_UCSC_mm9.tar.gz
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/Ensembl/NCBIM37/Mus_musculus_Ensembl_NCBIM37.tar.gz

echo "==================== Uncompressing human iGenome... ===================="
tar -zxvf $IGENOMES_HOME/Homo_sapiens_UCSC_hg18.tar.gz
tar -zxvf $IGENOMES_HOME/Homo_sapiens_UCSC_hg19.tar.gz
tar -zxvf $IGENOMES_HOME/Homo_sapiens_Ensembl_GRCh37.tar.gz

echo "==================== Uncompressing mouse iGenome... ===================="
tar -zxvf $IGENOMES_HOME/Mus_musculus_UCSC_mm9.tar.gz
tar -zxvf $IGENOMES_HOME/Mus_musculus_UCSC_mm10.tar.gz
tar -zxvf $IGENOMES_HOME/Mus_musculus_Ensembl_NCBIM37.tar.gz

chmod -R 777 $IGENOMES_HOME

cd $CUR_DIR

