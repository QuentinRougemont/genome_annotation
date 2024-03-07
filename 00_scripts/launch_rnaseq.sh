#!/bin/bash 

source ../config/config

echo rnaseqlist is $RNAseqlist 
#microscript to run all rnaseq steps from read trimming to read mapping and count
#RNAseq=YES
haplotype=$1

genome=03_genome/"$haplotype".fa*

#check that the file in 03_genome exist:

#check that the raw_data for RNAseq are present:

if [ -f Trimmomatic-0.39/trimmomatic-0.39.jar ] ; then 
    echo "found trimmomatic jar " ; 
else 
    echo "trimmmatic jar not found\nwill attempt to download" 
        wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
        unzip Trimmomatic-0.39.zip
fi

#launch gmap :
../00_scripts/02_gmap.sh $genome 

echo "trimming read for RNAseq" 
#paste <(ls 01_raw_data/*1.f**q.gz ) <(ls 01_raw_data/*2.f**q.gz ) > file1file2.tmp

ncol=$(awk '{print NF}'  $RNAseqlist |uniq)

if [[ $ncol = 2 ]] ; then
    #assuming PE:
    echo "running trimmomatic assuming reads are Paired-End" 
    while IFS=$'\t' read -r -a read ; 
    do 
        ../00_scripts/01_trimmomatic_PE.sh ${read[0]} ${read[1]}  
    done < $RNAseqlist #  file1file2.tmp 

    if [ $? -eq 0 ]; then
        echo $trimmomatic complete
    echo "counting the number of retained reads"        #rm file1file2.tmp 
    ../00_scripts/utility_scripts/count_read_fastq.sh 02_trimmed/*gz > read_count.txt

    else
        echo -e "\n#ERROR : Runnning trimmomatic failed. please check your input files"
        exit 1
    fi

    #launch gsnap - samtools and read count:
    for read1 in $(ls 02_trimmed/*R1.paired.fastq.gz ) ; do
        ../00_scripts/03_gsnap_PE.sh $genome $read1
    done 

else
    #assuming SE:
    echo "running trimmomatic" 
    while IFS=$'\t' read -r -a read ; 
    do 
        ../00_scripts/01_trimmomatic_SE.sh ${read[0]} 2>&1 |tee trimmo_"${read[0]}"_log  
    done < $RNAseqlist #  file1file2.tmp 
    
    if [ $? -eq 0 ]; then
        echo $trimmomatic complete
        #rm file1file2.tmp 
    else
        echo -e "\n#ERROR : Runnning trimmomatic failed. please check your input files"
        exit 1
    fi

    #launch gsnap - samtools and read count:
    for read1 in $(ls 02_trimmed/*R1.paired.fastq.gz ) ; do
        ../00_scripts/03_gsnap_SE.sh $genome $read1 2>&1 |tee gsnap_"$read1"_log
    done 

fi

rm -rf Trimmomatic-0.39  
