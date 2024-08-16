#!/bin/bash

source ../config/colors
#Date: March-2023
#Author: QR
#Purpose: script to blast the CDS against A1 and A2 of M-lagerheimii-1253-A1/A2 

input=$1 #fasta file containing CDS nucleotid
name=$(basename "$input" )

if [ -z "${input}" ] 
then
    echo "FATAL ERROR: no input file provided"
    exit
fi

databaseA1=$2
diamond makedb --in "${databaseA1}" -d  "${databaseA1}" 

diamond blastx -d "$databaseA1" \
        -q "$input" --ultra-sensitive\
        --outfmt 6 qseqid sseqid pident length qstrand mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq \
        --threads 20 \
        -o matches."$name".on_"$(basename "$databaseA1" )".tsv 

#databaseA2=/scratch/qrougemont/diamond_blast_resultsA1_A2_HD_PR/HD_PR_MvSv_A2.prot.fa
databaseA2=$3 
diamond makedb --in "${databaseA2}" -d  "${databaseA2}" 

diamond blastx -d "$databaseA2" \
        -q "$input" --ultra-sensitive\
        --outfmt 6 qseqid sseqid pident length qstrand mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq \
        --threads 20 \
        -o matches."$name".on_"$(basename "$databaseA2" )".tsv 

echo -e "\n${BLU}---- diamond finished----\n"

#----------------- optional:  ---------------------------
#--------------- uncomment this to search HD/PR ---------
#test if HD and PR are present finally:
#RED='\033[0;31m'
#NC='\033[0m' # No Color
#echo -e "${RED}----searching HD----"
#cat matches.*.tsv |awk '$12<1e-8 && $2 ~/HD/ {print $1"\t"$2"\t"$3"\t"$12}' |LC_ALL=C sort -g -k4 |tee > HD.searches
#echo -e "HD1 matches are :"
#grep "HD1" HD.searches 
#echo -e "HD2 matches are :"
#grep "HD2" HD.searches 
#echo -e "--------------------\n\n\n"
#
#echo -e "${RED}----searching PR----"
#cat matches.*.tsv |awk '$12<1e-8 && $2 ~/PR/ {print $1"\t"$2"\t"$3"\t"$12}' |LC_ALL=C sort -g -k4 |tee > PR.searches
#echo -e "PR match is :"
#grep "PR" PR.searches 
echo -e "--------------------\n${NC}"


