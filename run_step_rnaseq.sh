#!/bin/bash 

#microscriot to run all rnaseq steps from read trimming to read mapping and count
#RNAseq=YES
haplotype=$1

cd $haplotype
genome=03_genome/"$haplotype".fa*

#check that the file in 03_genome exist:

#check that the raw_data for RNAseq are present:

	
echo "trimming read for RNAseq" 
paste <(ls 01_raw_data/*1.f**q.gz ) <(ls 01_raw_data/*2.f**q.gz ) > file1file2.tmp

echo "running trimmomatic" 
while IFS=$'\t' read -r -a read ; 
do 
    ./00_scripts/01_trimmomatic.sh ${read[0]} ${read2[1]}  
done <  file1file2.tmp 

if [ $? -eq 0 ]; then
   echo $trimmomatic complete
   rm file1file2.tmp 
else
  echo -e "\n#ERROR : Runnning trimmomatic failed. please check your input files"
  exit 1
fi

#launch gmap :
./00_scripts/02_gmap.sh $genome 

#launch gsnap - samtools and read count:
for read1 in $(ls 02_trimmed/*R1.paired.fastq.gz ) ; do
   ./00_scripts/03_gsnap.sh $genome $read1
done 


