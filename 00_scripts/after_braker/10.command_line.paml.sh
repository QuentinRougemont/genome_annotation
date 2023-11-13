#!/bin/bash

#to run me: ./10.command_line.paml.sh species1.cds.fa species2.cds.fa  2>&1 |tee log
###############################################################

#To Do: here insert a help file with command lines ####################
#To Do: capture error and exit code at putative error places

outgroup=$1 #something ending in "spliced_cds.fa" generated from the previous scripts 
ingroup=$2 

###########################################################

in_base=$(basename ${ingroup%.spliced_cds.fa})
echo ingroup base name is $in_base
ou_base=$(basename ${outgroup%.spliced_cds.fa})
echo outgroup base name is $ou_base

sed -i 's/ CDS=.*$//g' $outgroup
sed -i 's/ CDS=.*$//g' $ingroup

#from orthofinder:
awk 'NF==4 && (( $2 ~/contig_8/ || $2 ~/contig_11/ ))  && $3 !~/lag/ ' orthofinder/Results_*/Orthogroups/Orthogroups.txt > single.copy.orthologs

dos2unix single.copy.orthologs

rm HD_and_PR* 2>/dev/null


cut -d " " -f3 single.copy.orthologs > HD_and_PR.$ou_base.txt
cut -d " " -f4 single.copy.orthologs  > HD_and_PR.$in_base.txt


##---  linearise the cds file ---#
cat "${outgroup}" | \
	awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}'  > "$ou_base".linearised.cds
cat "$ingroup" | \
	awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}'  > "$in_base".linearised.cds


##---- recover the wanted sequence HD and PR for in and outgroup #
rm sorted* 2>/dev/null
while read -r pattern ; 
do 
    grep -w -A1 "$pattern" "$ou_base".linearised.cds  >> sorted."$ou_base".wanted_cds.fa ;
done < HD_and_PR.$ou_base.txt 
#
while read -r pattern ; 
do 
    grep -w -A1 "$pattern" "$in_base".linearised.cds >> sorted."$in_base".wanted_cds.fa ; 
done < HD_and_PR.$in_base.txt 
#

### ---- then run paml and dnds.sh  ----- #
#paml does not like long fasta name so we shorten them for the sake of computation 
#grep ">" sorted."$ou_base".wanted_cds.fa >> id.A1
#grep ">" sorted."$in_base".wanted_cds.fa >> id.A2
#sed -i 's/Myosoton-A1_/MyoA1/g' sorted* 
#sed -i 's/Myosoton-A2_/MyoA2/g' sorted* 
#sed -i 's/_//g' sorted* 
#sed -i 's/-//g' sorted* 

#------ part specific to PAML --------------------- #
fasta1=sorted."$ou_base".wanted_cds.fa  
fasta2=sorted."$in_base".wanted_cds.fa

f1=$(basename $fasta1)
f2=$(basename $fasta2)
newf1=${f1%.fa**}.nostopcodon.fasta
newf2=${f2%.fa**}.nostopcodon.fasta
#here add test to linearize the fasta if necessary
#awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' 

##1 ----- remove stop codon -------
awk '{if($1 !~ /^>/) {if( substr($0, length($0)-2, length($0)) ~ /(TGA|TAG|TAA)/){ print substr($0, 1, length($0)-3)} else  {print $0 }}else{print $0}}' $fasta1 > $newf1 
awk '{if($1 !~ /^>/) {if( substr($0, length($0)-2, length($0)) ~ /(TGA|TAG|TAA)/){ print substr($0, 1, length($0)-3)} else  {print $0 }}else{print $0}}' $fasta2 > $newf2 

##2 ------ split and cat pairwise sequence -------
#split 
rm -rf sequence_files 2>/dev/null
rm wanted_sequence 2>/dev/null
rm ID1 ID2 2>/dev/null

mkdir sequence_files

grep ">"  $newf1 > ID1
grep ">"  $newf2 > ID2
paste ID1 ID2 |sed 's/>//g' > wanted_sequence

while IFS=$'\t' read -r -a line 
do 
	mkdir sequence_files/tmp.${line[0]}.vs.${line[1]}

	grep -A1 ${line[0]}  $newf1 > sequence_files/tmp.${line[0]}.vs.${line[1]}/sequence.fasta
	grep -A1 ${line[1]}  $newf2 >> sequence_files/tmp.${line[0]}.vs.${line[1]}/sequence.fasta

	translatorx_vLocal.pl -i sequence_files/tmp.${line[0]}.vs.${line[1]}/sequence.fasta -o sequence_files/tmp.${line[0]}.vs.${line[1]}/results

	cp 00.scripts/yn00_template.ctl sequence_files/tmp.${line[0]}.vs.${line[1]}/

        cd sequence_files/tmp.${line[0]}.vs.${line[1]}/

	path=$(pwd)
	echo $path
	sed -i "s|PATH|$path|g" 	 yn00_template.ctl #"

	yn00 yn00_template.ctl #

        cd ../../
	awk '/\+\-/ && !/(dS|SE)/ {split(FILENAME, a, "."); print $(NF-2), $(NF), $(NF-5), $(NF-3),"'${line[0]}'","'${line[1]}'"}'  sequence_files/tmp.${line[0]}.vs.${line[1]}/out_yn00_orthogp >  sequence_files/tmp.${line[0]}.vs.${line[1]}/resultat_Yang_Nielsen_2000_method.orthogp.txt 


done < wanted_sequence
cat sequence_files/tmp.*/resultat_Yang_Nielsen_2000_method.orthogp.txt >> results_YN.txt
