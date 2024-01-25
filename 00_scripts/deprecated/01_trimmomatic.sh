#!/bin/bash
#date: 14-11-22
#Author:QR
#script to run trimmomatic

#required arguments: name of read1.fastq.gz and read2.fastq.gz
if [ $# -ne 2 ]; then
    echo "USAGE: $0 R1.fastq.gz R2.fastq.gz file "
    echo "Expecting fastq files from read1 and read2 as input"
    echo "Extension should be '*fastq.gz'"
    exit 1
else
    file1=$1
    file2=$2 
    echo "fastq file1 is : ${file1}"
    echo "fastq file2 is : ${file2}"
    echo " "
fi

#base=$(basename $file1)
base=$(basename ${file1%[_.]*[._]f**gz} )

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="99_log_files"
cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

#create folder if not existent:
mkdir $LOG_FOLDER 2>/dev/null
mkdir 02_trimmed  2>/dev/null

ADAPTERFILE="Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa"
NCPU=8

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip

java -jar -Xmx10G Trimmomatic-0.39/trimmomatic-0.39.jar PE \
        -threads 8 \
        -phred33 \
        "$file1" \
        "$file2" \
        02_trimmed/"$base"_R1.paired.fastq.gz \
        02_trimmed/"$base"_R1.single.fastq.gz \
        02_trimmed/"$base"_R2.paired.fastq.gz \
        02_trimmed/"$base"_R2.single.fastq.gz \
        ILLUMINACLIP:"$ADAPTERFILE":2:30:10:2:True \
	HEADCROP:9 \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:30:30 \
        MINLEN:100 2> $LOG_FOLDER/log.trimmomatic.pe."$base"."$TIMESTAMP"


