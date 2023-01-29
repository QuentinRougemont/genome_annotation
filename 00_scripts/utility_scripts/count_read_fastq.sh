#!/bin/bash
if [ $# -lt 1 ] ; then
  echo ""
  echo "usage: count_fastq.sh [fastq_file1] <fastq_file2> ..|| *.fastq"
  echo "counts the number of reads in a fastq file"
  echo ""
 exit 0
fi
     
fastq=${@};
for i in ${fastq[@]}
do
    lines=$(zcat $i |wc -l |cut -d " " -f 1)
    count=$(($lines / 4))
    echo -n -e "\t$i : "
    echo "$count"  
done


