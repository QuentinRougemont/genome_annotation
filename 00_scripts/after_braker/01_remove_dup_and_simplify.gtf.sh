#!/bin/bash

#to run me on all species do:
#for file in $(find . -name "*ok.gtf" ) ; do base=$(basename $file ) ; path=${file%/*} ; echo -e "cd $path\n../../00_script/02_remove_redundant_genemark.sh $base ; \ncd ../../" >> clean.gtf.tmp  ; done
#bash clean.gtf.tmp 

gtf=$1
species=${gtf%.ok.gtf}

echo -e "species is $species"

mkdir 01_species_cds 2>/dev/null
mkdir 02_species_prot 2>/dev/null

#extract the cds and protein from it:
output="$species"_cds.fa
echo output cds is $output
gffread -w 01_species_cds/$output -g 03_genome/$species.fa $gtf

#here capture error if genome and gtf do not match and exit with error 1 

#here if any error occurs capture it!
#then convert also the file to its cds:
echo "translate CDS into amino acid "
transeq -sequence 01_species_cds/$output -outseq 02_species_prot/$species.prot

eval "$(conda shell.bash hook)"
conda activate /home/genouest/cnrs_umr5175/qrougemont/mes_envs/samtools1.18/
samtools faidx 02_species_prot/$species.prot


echo "-----------------------------------------------------------------"
echo "extract longest transcript" 
echo "-----------------------------------------------------------------"

cd 02_species_prot
#assumption : all transcript finishes by ".t1, .t2, .t3 so the dot (.) iis the delimiter
awk -F "." '{print $1"\t"$0}' $species.prot.fai |\
        awk '$3>max[$1]{max[$1]=$3; row[$1]=$2} END{for (i in row) print row[i]}' > longest.transcript.tmp

#linearize file so that the next command will work:
awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' "$species".prot > "$species".prot.lin.fasta

grep -A1 -Ff longest.transcript.tmp "$species".prot.lin.fasta > "$species".longest_transcript.fa
rm longest.transcript.tmp
rm $species.prot.fai
cd ../

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#clean based on CDS now:

echo "-----------------------------------------------------------------"
echo "remove redundant Gene in the CDS and Protein fasta files"
echo "-----------------------------------------------------------------"

#the purpose of the line of code is to remove redundant GeneMark annotation, that is GeneMark gene with same start and end position as those from Augustus:
sort -k1,1 -k4,4n  $gtf |cut -f 1-5,9 |awk '$3~"CDS" {print $1"_"$4"_"$5"\t"$1"_"$4"\t"$1"_"$5"\t"$0}' > tmp1
awk 'NR == FNR {count[$1]++; next} count[$1]>1' tmp1 tmp1 |grep "GeneMark" |cut -f 9 |cut -d ";" -f 1 |uniq |sed -e 's/transcript_id //' -e 's/"//g' -e 's/ //g' > toremove
awk 'NR == FNR {count[$2]++; next} count[$2]>1' tmp1 tmp1 |grep "GeneMark" |cut -f 9 |cut -d ";" -f 1 |uniq |sed -e 's/transcript_id //' -e 's/"//g' -e 's/ //g' >> toremove
awk 'NR == FNR {count[$3]++; next} count[$3]>1' tmp1 tmp1 |grep "GeneMark" |cut -f 9 |cut -d ";" -f 1 |uniq |sed -e 's/transcript_id //' -e 's/"//g' -e 's/ //g' >> toremove


rem=$(wc -l toremove |awk '{print $1}')
echo -e "fin a total of $rem redondant transcript\n------------------------------------"
mkdir 03_cleaned_prot 2>/dev/null
sed 's/ CDS=.*//g' 02_species_prot/"$species".longest_transcript.fa|sed '/^--$/d' |awk -vRS=">" -vOFS="\t" ' {$1=">"$1; print $0}' |grep -vFf toremove |sed 's/\t/\n/g' |sed 1d  > 03_cleaned_prot/"$species"_prot.fa
#rm tmp1 toremove

cut -f 9 "$gtf" |sed 's/^gene_id//g' |sed 's/transcript_id//g' |cut -d ";" -f 1 |sed 's/.t[1-9]//g' |sort |uniq -c |awk '$1==1 {print}'  > fragmented_gene.txt 

loss=$(wc -l fragmented_gene.txt |awk '{print $1}' )
echo " there is $loss fragmented gene"
echo "finished"


#species=$1 #name of the species 

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


gffread -w "$species".spliced_cds.fa -g 03_genome/$species.fa "$species".longest_transcript.gtf
echo "translate CDS into amino acid "
transeq -sequence "$species".spliced_cds.fa -outseq "$species"_prot.fa

