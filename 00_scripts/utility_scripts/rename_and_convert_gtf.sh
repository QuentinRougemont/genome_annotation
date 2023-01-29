#!/bin/bash
#Date: 18-01-23
#Author:QR
#purpose script to  reshape the gft and extract the cds based on gtf and fasta file
if [ $# -ne 3 ]; then
    echo "USAGE: $0 gtffile fasta_file prefix"
    echo "Expecting the following values on the command line, in that order"
    echo "gtffile : name of the gtf file"
    echo "fastafile: name of the fasta "
    echo "prefix: prefix for the gene in the new gtf"
    exit 1
else 
    gtf=$1 #gtf file
    fasta=$2 #folder where final sfs will appear
    prefix=$3
    echo "gtf file  is $gtf"
    echo "fasta file is $fasta"
    echo "prefix for the genes will be $prefix"
    echo "\n"
fi

#first fix the input name with rename_gtf.py from TSEBRA
#capture variable:
input=$gtf
genome=$fasta

output=${input%.gtf}.renamed.gtf

# -- first fix the input name with rename_gtf.py from TSEBRA
rename_gtf.py --gtf $input --prefix Mintermedium --translation_tab translation.tab --out ${output}

input=$output 
bedfile=${input%.gtf}.bed 
gtffasta=${input%.gtf}.cds.fasta

# -- convert the gtf2gff using augustus --- #
gtf2gff.pl <$input --out="${input%.gtf}".gff3 --gff3

#----- step1  convert gft2bed ---- #
echo input gtf for gtf2bed is $input
gtf2bed  <$input  > $bedfile 


#----- step2  extract fasta sequence matching bed coordinates -----#
#keep only the exons
awk '$8=="CDS" {print $0}' $bedfile > tmp
bedtools getfasta -fi $genome -bed tmp -fo  $gtffasta #split -s 

#clean-up
rm tmp

