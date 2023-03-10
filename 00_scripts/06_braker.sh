#!/bin/bash
#PURPOSE: script to run braker
#AUTHOR: QR
#Date updated: 10-03-2023
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT


#  ---- external data required arguments --------------- #
if [ $# -ne 3 ] ; then
        echo "USAGE: $0 reference_genome species RNAseq(YES/NO)" 
        echo -e "Expecting the following parameters:\n
                1 - the reference genome\n
                2 - a species name\n
                3 - a string YES/NO for wether RNAseq should be considered or not\n\n"
        exit 1
else
    genome=$1 #the reference genome here! "genome.wholemask.fa"
    species=$2
    RNAseq=$3 #YES/NO
    NCPUS=$4
fi

base=$(basename $genome )

if [[ -z "$NCPUS" ]]
then
    NCPUS=8
fi

TIME=$(date +%Y-%m-%d_%Hh%Mm%Ss)
alnBAM="04_mapped/RNAseq.bam"
relatProt="/path/to/closely_related/protein.fa" 

if [[ -d 06_braker ]] 
then
	echo "WARNING directory 06_braker already exists! check its content first"
	exit 1
fi
## --------- step 1 : BRAKER WITH RNA SEQ  ---------  ##
#Run braker with RNAseq 
#round 1:
wd=06_braker/rnaseq
mkdir -p $wd

if [[ $RNAseq = "YES" ]]
then
    braker.pl --species="$species" --fungus --genome="$genome" --cores="$NCPUS"  --softmasking --bam="$alnBAM" --workingdir=$wd 
fi 


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
braker.pl --species=$species"_refprot" --fungus --genome="$genome" --cores="$NCPUS" --prot_seq=$relatProt --workingdir=$wd --softmasking  

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


