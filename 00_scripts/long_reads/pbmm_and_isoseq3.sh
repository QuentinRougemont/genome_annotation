#!/bin/bash
#SBATCH --job-name=isoseq
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#script to run pbmm + isoseq3 + braker(genemark/Augustus)
#AUTHOR: QR
#Date: 10-03-2023


conda activate isoseq3

## -- external variables -- ##
input=$1 #Isoseq data #in fastq.gz format!! 
name=$(basename ${input%.fastq.gz} )
genome=$2 #reference genome 


## -- prepare architecture -- ##
aln=07_isoseq_mapped
mkdir $aln 2>/dev/null 
mkdir 08_collapsed 2>/dev/null

echo -e "input is $input \n"

echo "running pbmm2"
pbmm2 align \
       --preset ISOSEQ \
	--sort $input \
	$genome \
	$aln/$name.bam

bam="$name".bam

echo "collapse read from $bam "
isoseq3 collapse "$aln"/"$bam" 08_collapsed/"$name".gff

conda activate anaCogent
echo "#Run GeneMarkS-T to predict protein-coding regions in the transcripts:"
~/software/Augustus-3.5.0/scripts/stringtie2fa.py -g "$genome" -f 08_collapsed/"$name".gff -o 08_collpased/"$name".fa

echo "strintie2fa successfully done"

echo "running gmst now"
#export PATH="/scratch/qrougemont/software/BRAKER/scripts/:$PATH" #note : already in path
gmst.pl --strand direct 08_collapsed/"$name".fa.mrna --output 08_collapsed/gmst."$name".out --format GFF

echo "gmst done"
echo -e"\n Use the GeneMarkS-T coordinates and the long-read transcripts to create a gene set in GTF format. \n"
gmst2globalCoords.py \
	-t 08_collapsed/"$name".gff \
	-p 08_collapsed/gmst."$name".out \
	-o 08_collapsed/gmst."$name".global.gtf \
	-g $genome

