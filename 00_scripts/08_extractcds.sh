#!/bin/bash
#miniscript to extract cds from fasta and convert output into protein file
#Author: QR
#Date: 23-01-23

if [ $# -ne 2  ]; then
    echo "USAGE: $0 gtf_file genome_file "
    echo "Expecting the gtf file from TSEBRA and the genome" 
    exit 1
else
    gtf=$1
    genome=$2
fi


#script to extract fasta sequence given the gtf coordinates
#source code: http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread

gtf=$1
genome=$2
output=${gtf%.gtf}.cds.fa

gffread -w $output -g $genome $gtf 

#then convert also the file to its cds:
input=$output
transeq -sequence $input -outseq ${input%.fa**}.prot
