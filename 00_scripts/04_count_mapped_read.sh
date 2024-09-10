#!/bin/bash
if [ $# -lt 1 ] ; then
  echo ""
  echo "usage: 03_count_mapped_read.sh  [bam_file1] <bam_file2> ..|| *.bam"
  echo "counts the number of mapped reads in a bam file"
  echo ""
 exit 0
fi
     
bam="${!@}"
for i in "${bam[@]}"
do
    if [ ! -s  comptage_brute."${i%.bam}".txt ] 
    then
        #counting
        samtools view -c "$i" |awk -v var="$i" 'END {print var"\t"$1}' > comptage_brute."${i%.bam}".txt
        samtools view -F 0x4 "$i" | cut -f 1 | sort | uniq | wc -l |\
            awk -v var="$i" 'END {print var"\t"$1}' > comptage_F04."${i%.bam}".txt ;
    else
       echo "counting already done"
    fi 
done
