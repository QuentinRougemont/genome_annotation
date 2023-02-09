#!/bin/bash
#script to run braker
#AUTHOR: QR
#Date: 01-01-2023

#external data
if [ $# -ne 1  ]; then
    echo "USAGE: $0 reference_genome database name rm_unknwon(yes/no)"
    echo -e "Expecting the reference_genome as input\n"
    exit 1
else
    genome=$1
    echo -e "reference genome is $genome \n"
    echo " "
fi

base=$(basename $genome )

TIME=$(date +%Y-%m-%d_%Hh%Mm%Ss)
alnBAM="04_mapped/Aligned.out.ss.bam"
NCPUS=8
relatProt="/path/to/closely_related/protein.fa" 

## --------- step 1 : BRAKER WITH RNA SEQ  ---------  ##
#Run braker with RNAseq 
#round 1:
wd=06_braker/rnaseq
mkdir -p $wd
braker.pl --species="$base" --fungus --genome="$genome" --cores="$NCPUS"  --softmasking --bam="$alnBAM" --workingdir=$wd 


##  --------- step 2 : BRAKER WITH REFERENCE DATABASE USING THREE ROUNDS --------- ## 

#prepare architecture:
FOLDER1=06_braker/round1_braker_on_refprot #_$TIME
FOLDER2=06_braker/round2_braker_on_refprot #_$TIME
FOLDER3=06_braker/round3_braker_on_refprot #_$TIME 
FOLDER4=06_braker/round4_braker_on_refprot #_$TIME
FOLDER5=06_braker/round5_braker_on_refprot #_$TIME

mkdir -p $FOLDER1 $FOLDER2 $FOLDER3 $FOLDER4 $FOLDER5 2>/dev/null

# ---- round 1 : ----- #
wd=${FOLDER1}
braker.pl --species=$base"_refprot" --fungus --genome="$genome" --cores="$NCPUS" --prot_seq=$relatProt --workingdir=$wd --softmasking  

# ---- round 2 : ----- #
wd=${FOLDER2}
braker.pl  --fungus --genome="$genome" --cores="$NCPUS" --workingdir=$wd --hints=${FOLDER1}/hintsfile.gff --softmasking 

# ---- round 3 : ----- #
wd=${FOLDER3}
braker.pl  --fungus --genome="$genome" --cores="$NCPUS" --workingdir=$wd --hints=${FOLDER2}/hintsfile.gff --softmasking  

wd=${FOLDER4}
braker.pl  --fungus --genome="$genome" --cores="$NCPUS" --workingdir=$wd --hints=${FOLDER3}/hintsfile.gff --softmasking  

wd=${FOLDER5}
braker.pl  --fungus --genome="$genome" --cores="$NCPUS" --workingdir=$wd --hints=${FOLDER4}/hintsfile.gff --softmasking  


