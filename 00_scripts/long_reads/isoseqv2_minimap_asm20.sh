#!/bin/bash

#IsoSeq pipeline 
#run minimap - cupcake - stringtie and gmst from gene mark to find gene 
#use this script if IsoSeq read differs from the reference genome

##  ------------------------ general parameters --------------------------------  ##
while [ $# -gt 0 ] ; do
  case $1 in
    -i | --input) input="$2" ;echo "the species fastq is $input" >&2;;
    -g | --genome) genome="$2" ; echo "reference genome is $genome" >&2;;
    -h | --help) echo -e "Option required:
    -i/--input \t the species Isoseq fastq to process
    " >&2;exit 1;;
    esac
    shift
done

if [ -z "$input" ] || [ -z "$genome" ] ; then
	echo >&2 "Fatal error: input file require (isoseq fastq)\n
	see manual with -h or --help"
exit 2
fi

set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

#genome=03_genome/your_genome.fasta 

echo "input is $input" 
gunzip $input
isoseq=${input%.gz}
name=$(basename ${input%.fastq.gz}) 

LR=07_isoseq_mapped/

mkdir $LR 2>/dev/null

echo "running minimaps on $isoseq "


NCPUS=20

minimap2 -ax splice -asm20 -t $NCPUS -uf --secondary=no $genome $isoseq > $LR/$name.sam
sort -k 3,3 -k 4,4n $LR/$name.sam > $LR/$name.sorted.sam

#collapse redundant IsoForms with CupCake:
eval "$(conda shell.bash hook)"
conda activate anaCogent

echo "will collapse isoform"
collapse_isoforms_by_sam.py \
        --input $isoseq \
	--fq  \
       	-s $LR/$name.sorted.sam  \
	-c 0.99 -i 0.95\
	 --dun-merge-5-shorter \
	 -o $LR/cupcake."$name" # --cpus 24\

echo "sucessffuly collapsed redundand isoforms"

#You will find the collapsed transcripts in cupcake.collapsed.gff .
#Run GeneMarkS-T to predict protein-coding regions in the transcripts:
stringtie2fa.py -g "$genome" \
		-f $LR/cupcake."$name".collapsed.gff \
		-o $LR/cupcake."$name".fa

echo "strintie2fa successfully done"

echo "running gmst now"
gmst.pl --strand direct $LR/cupcake."$name".fa.mrna --output $LR/gmst."$name".out --format GFF

echo "gmst done"
echo -e "\n Use the GeneMarkS-T coordinates and the long-read transcripts to create a gene set in GTF format. \n"
gmst2globalCoords.py -t $LR/cupcake."$name".collpased.gff -p $LR/gmst."$name".out -o $LR/gmst."$name".global.gtf -g $genome

#gffread and transeq to extract the final prot
mkdir 08_transcripts
gffread -g $genome -w 08_transcripts/$name.cds 07_isoseq_mapped/gmst.$name.global.gtf
#eventually run transeq here
