#!/bin/bash
#PURPOSE: script to run braker
#AUTHOR: QR
#Date updated: 10-03-2023
set -e
source ../config/config

eval "$(conda shell.bash hook)"
conda activate braker_env

echo "relatedProt is $RelatedProt"

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

#  ---- external data required arguments --------------- #
if [ $# -ne 4 ] ; then
	echo "USAGE: $0 reference_genome species RNAseq(YES/NO) fungus(YES/NO)" 
        echo -e "Expecting the following parameters:\n
                1 - the reference genome\n
                2 - a species name\n
                3 - a string YES/NO for wether RNAseq should be considered or not\n\n
		4 - a string YES/NO stating wether data are from fungus or not \n\n" 
        exit 1
else
    genome=$1 #the reference genome here! "genome.wholemask.fa"
    species=$2
    RNAseq=$3 #YES/NO
    fungus=$4
    NCPUS=$5
fi

base=$(basename $genome )

if [[ -z "$NCPUS" ]]
then
    NCPUS=8
fi

TIME=$(date +%Y-%m-%d_%Hh%Mm%Ss)

#previously we used a concatenated and filtered bam.
#however, providing a list of all bam gives similart results with a significant gain in time:
alnBAM=$(echo 04_mapped/*sorted.bam |sed 's/ /,/g' )

#simpl --species="$species"_round1ifying this: 
#this hard-coded variable should be repl --species="$species"_round1aced in the config file":
relatProt="$RelatedProt" 

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
    if [[ $fungus = "YES" ]]
    then
            braker.pl --species="$species"_round1 --species="$species" --fungus --genome="$genome" --threads="$NCPUS"  --softmasking --bam="$alnBAM" --workingdir=$wd 
    else
            braker.pl --species="$species"_round1 --species="$species" --genome="$genome" --threads="$NCPUS"  --softmasking --bam="$alnBAM" --workingdir=$wd 
    fi
fi 


##  --------- step 2 : BRAKER WITH REFERENCE DATABASE USING THREE ROUNDS --------- ## 

#prepare architecture:
FOLDER1=06_braker/round1_braker_on_refprot #_$TIME
FOLDER2=06_braker/round2_braker_on_refprot #_$TIME
FOLDER3=06_braker/round3_braker_on_refprot #_$TIME 
FOLDER4=06_braker/round4_braker_on_refprot #_$TIME
FOLDER5=06_braker/round5_braker_on_refprot #_$TIME

mkdir -p $FOLDER1 $FOLDER2 $FOLDER3 $FOLDER4 $FOLDER5 2>/dev/null

echo "----------- round 1 ------------" 
echo AUGUSTUS_SCRIPTS_PATH is $AUGUSTUS_SCRIPTS_PATH 
echo AUGUSTUS_BINS_PATH is $AUGUSTUS_BIN_PATH 
echo AUGUSTUS_CONFIG_PATH is $AUGUSTUS_CONFIG_PATH 


wd=${FOLDER1}

if [[ $fungus = "YES" ]]
then
    braker.pl --species="$species"_round1  --genome="$genome" --threads="$NCPUS"  --softmasking --prot_seq=$relatProt --workingdir=$wd --fungus  
else
    braker.pl --species="$species"_round1  --genome="$genome" --threads="$NCPUS"  --softmasking --prot_seq=$relatProt --workingdir=$wd 
fi

echo "----------- round 2 ------------" 
wd=${FOLDER2}
if [[ $fungus = "YES" ]]
then
    braker.pl --species="$species"_round2 --genome="$genome" --threads="$NCPUS"  --softmasking --prot_seq=$relatProt --workingdir=$wd --fungus #--hints=${FOLDER1}/hintsfile.gff 
else
    braker.pl --species="$species"_round2 --genome="$genome" --threads="$NCPUS"  --softmasking --prot_seq=$relatProt --workingdir=$wd #--hints=${FOLDER1}/hintsfile.gff 
fi

echo "----------- round 3 ------------" 
wd=${FOLDER3}
if [[ $fungus = "YES" ]]
then
    braker.pl --species="$species"_round3 --genome="$genome" --threads="$NCPUS"  --softmasking --prot_seq=$relatProt --workingdir=$wd --fungus #--hints=${FOLDER2}/hintsfile.gff 
else
    braker.pl --species="$species"_round3 --genome="$genome" --threads="$NCPUS"  --softmasking --prot_seq=$relatProt --workingdir=$wd #--hints=${FOLDER2}/hintsfile.gff 
fi

echo "----------- round 4 ------------" 
wd=${FOLDER4}
if [[ $fungus = "YES" ]]
then
    braker.pl --species="$species"_round4 --genome="$genome" --threads="$NCPUS"  --softmasking --prot_seq=$relatProt --workingdir=$wd --fungus #--hints=${FOLDER3}/hintsfile.gff 
else
    braker.pl --species="$species"_round4 --genome="$genome" --threads="$NCPUS"  --softmasking --prot_seq=$relatProt --workingdir=$wd #--hints=${FOLDER3}/hintsfile.gff 
fi

echo "----------- round 5 ------------" 
wd=${FOLDER5}
if [[ $fungus = "YES" ]]
then
    braker.pl --species="$species"_round5 --genome="$genome" --threads="$NCPUS"  --softmasking --prot_seq=$relatProt --workingdir=$wd --fungus #--hints=${FOLDER4}/hintsfile.gff 
else
    braker.pl --species="$species"_round5 --genome="$genome" --threads="$NCPUS"  --softmasking --prot_seq=$relatProt --workingdir=$wd #--hints=${FOLDER4}/hintsfile.gff 
fi
echo "----------- finished------------" 

