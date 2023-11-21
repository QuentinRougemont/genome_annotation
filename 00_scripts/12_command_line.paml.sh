#!/bin/bash

#to run me: ./10.command_line.paml.sh species1.cds.fa species2.cds.fa  2>&1 |tee log
###############################################################

#To Do: here insert a help file with command lines ####################
#To Do: capture error and exit code at putative error places

species1=$1 #name of species1 
species2=$2 #name of species2 
scaffold=$3 #list of scaffold on the ancestral species 
###########################################################

#cds file:
cdsfile1=$species1/08_best_run/$species1.spliced_cds.fa
cdsfile2=$species2/08_best_run/$species2.spliced_cds.fa

#species1=$(basename ${species2%.spliced_cds.fa})
#echo species2 base name is $species1
#species2=$(basename ${species1%.spliced_cds.fa})
#echo species1 base name is $species2

#remove the CDS length info that is introduced by gffread:
sed -i 's/ CDS=.*$//g' $cdsfile1
sed -i 's/ CDS=.*$//g' $cdsfile2


#-- get single copy orthologs from orthofinder ---
#-- criteria: we want 1:1:1 orthologs between the ancestralspecies:species1:species2
#awk 'NF==4 && (( $2 ~/contig_8/ || $2 ~/contig_11/ ))  && $3 !~/lag/ ' orthofinder/Results_*/Orthogroups/Orthogroups.txt > single.copy.orthologs
#may be also add a $2 /ancestral_species/ in awk ?

#awk -v s1="$species1" -v s2="$species2" 'NF==4 && ( $3 ~ s2  && $4 s1 || NF==4 && $3 ~ s1  && $4 s2) ' orthofinder/Results_*/Orthogroups/Orthogroups.txt > single.copy.orthologs

#strong assumption: species were provided alphabetically. 
#in orthofinder the speices appears in each column alphabetically 
#awk -v s1="$species1" -v s2="$species2" '{if (NF==4 &&  $3 ~ s1  && $4 s2 ){print $0}else if(NF==4 && $3 ~ s2  && $4 s1 ){print $1,$2,$4,$3}} ' genespace/orthofinder/Results_*/Orthogroups/Orthogroups.txt | \
#	grep -Ff $scaffold - > paml/single.copy.orthologs
grep -Ff genespace/orthofinder/Results_*/Orthogroups/Orthogroups_SingleCopyOrthologues.txt genespace/orthofinder/Results_*/Orthogroups/Orthogroups.txt |grep -Ff $scaffold - > paml/single.copy.orthologs


#correct the output that seemed to have bugs sometimes:
dos2unix paml/single.copy.orthologs

rm HD_and_PR* 2>/dev/null

cut -d " " -f3 paml/single.copy.orthologs > paml/HD_and_PR.$species1.txt
cut -d " " -f4 paml/single.copy.orthologs > paml/HD_and_PR.$species2.txt


##---  linearise the cds file ---#
cat "${cdsfile2}" | \
	awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}'  > paml/"$species2".linearised.cds
cat "$cdsfile1" | \
	awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}'  > paml/"$species1".linearised.cds


##---- recover the wanted sequence HD and PR for in and species1 #
cd paml

#all is run from paml folder now

rm sorted* 2>/dev/null
while read -r pattern ; 
do 
    grep -w -A1 "$pattern" "$species2".linearised.cds  >> sorted."$species2".wanted_cds.fa ;
done < HD_and_PR.$species2.txt 
#
while read -r pattern ; 
do 
    grep -w -A1 "$pattern" "$species1".linearised.cds >> sorted."$species1".wanted_cds.fa ; 
done < HD_and_PR.$species1.txt 
#

### ---- then run paml and dnds.sh  ----- #
#paml does not like long fasta name so we shorten them for the sake of computation 
#to do here

#------ part specific to PAML --------------------- #
fasta1=sorted."$species1".wanted_cds.fa  
fasta2=sorted."$species2".wanted_cds.fa

f1=$(basename $fasta1)
f2=$(basename $fasta2)
newf1=${f1%.fa**}.nostopcodon.fasta
newf2=${f2%.fa**}.nostopcodon.fasta

##1 ----- remove stop codon -------
awk '{if($1 !~ /^>/) {if( substr($0, length($0)-2, length($0)) ~ /(TGA|TAG|TAA)/){ print substr($0, 1, length($0)-3)} else  {print $0 }}else{print $0}}' ${fasta1} > ${newf1} 
awk '{if($1 !~ /^>/) {if( substr($0, length($0)-2, length($0)) ~ /(TGA|TAG|TAA)/){ print substr($0, 1, length($0)-3)} else  {print $0 }}else{print $0}}' ${fasta2} > ${newf2}

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

	translatorx_vLocal.pl -i sequence_files/tmp.${line[0]}.vs.${line[1]}/sequence.fasta -o sequence_files/tmp.${line[0]}.vs.${line[1]}/results 2>&1 |tee log.translator

	cp ../../00_scripts/yn00_template.ctl sequence_files/tmp.${line[0]}.vs.${line[1]}/

        cd sequence_files/tmp.${line[0]}.vs.${line[1]}/

	path=$(pwd)
	echo $path
	sed -i "s|PATH|$path|g" 	 yn00_template.ctl #"

	yn00 yn00_template.ctl  

        cd ../../
	awk '/\+\-/ && !/(dS|SE)/ {split(FILENAME, a, "."); print $(NF-2), $(NF), $(NF-5), $(NF-3),"'${line[0]}'","'${line[1]}'"}'  sequence_files/tmp.${line[0]}.vs.${line[1]}/out_yn00_orthogp >  sequence_files/tmp.${line[0]}.vs.${line[1]}/resultat_Yang_Nielsen_2000_method.orthogp.txt 


done < wanted_sequence 2>&1 |tee log.paml 
cat sequence_files/tmp.*/resultat_Yang_Nielsen_2000_method.orthogp.txt >> results_YN.txt
