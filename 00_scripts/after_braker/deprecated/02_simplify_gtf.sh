#!/bin/bash
#AUTHOR: QR
#DATE: 11-10-23
#Purpose: clean gtf to keeep only wanted cds (i.e. the longest transcript present in the protein folder"
species=$1 #name of the species 

#-------------------------------
gtf="$species".ok.gtf
protpath=03_cleaned_prot
#----------------------------------

#declare full path to input:
prot="$protpath"/"$species"_prot.fa

echo -e "there is $(wc -l $gtf |awk '{print $1}') lines in $gtf" 
#----------------------------------------------------
# subset our gtf to keep only the cds in the cds files!
grep ">" $prot |sed 's/>//g' |sed 's/_1$//g'  > wanted.cds.tmp 
sed 's/.t[0-9]$//g' wanted.cds.tmp > wanted.gene.tmp

#ok so here I tried many awk solution to have both the wanted transcript and the gene using NF == FNR, etc, nothing work so I used a
#two step process creating a p1 file for the whole stuf and p2 for the gene, which is not really what I wanted.
grep -Ff wanted.cds.tmp "$gtf" > p1 #"$name".sub.gtf
grep -f wanted.gene.tmp  <(awk '$3=="gene" ' $gtf)  > p3

#ideally I want to sort on the gene, then transcript, then CDS, then exon as well
cat p1 p3 |LC_ALL=C sort -k1,1 -k4,4n -k5,5n > "$species".longest_transcript.gtf

echo -e there is $(wc -l "$species".longest_transcript.gtf |awk '{print $1}' ) lines in "$species".longest_transcript.gtf 

