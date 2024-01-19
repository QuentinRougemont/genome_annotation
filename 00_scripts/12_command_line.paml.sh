#!/bin/bash

#purpose: compute ds based on paml using yn00 model
#to run me: ./12.command_line.paml.sh haplo1.cds.fa haplo2.cds.fa scaffold ancestral_genome 2>&1 |tee log

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo -e "master script to: \n run translatorX and paml "
   echo " "
   echo "Usage: $0 [-h1|-h2|-s|-a|-h|]"
   echo "options:"
   echo " -h|--help: Print this Help."
   echo " -h1|--haplo1: the name of the first  focal haplotype\t "
   echo " -h2|--haplo2: the name of the second focal haplotype\t "
   echo " -a |--ancestral_genome: the name of the ancestral haplo to infer orthology and plot gene order"
   echo " -s |--scaffold: the name of scaffold of interest"
   echo " "
   echo "dependancies: paml (yn00), translatorX "
}

############################################################
# Process the input options.                               #
############################################################
while [ $# -gt 0 ] ; do
  case $1 in
    -h1 | --haplo1) haplo1="$2" ; echo -e "haplotype 1 Name is ***${haplo1}*** \n" >&2;;
    -h2 | --haplo2) haplo2="$2" ; echo -e "haplotype 2 Name is ***${haplo2}*** \n" >&2;;
    -a  | --ancestral_genome) ancestral_genome="$2" ; echo -e "ancestral haplo  Name is ***${ancestral_genome}*** \n" >&2;;
    -s  | --scaffold) scaffold="$2"  ; echo -e "target scaffold is ${scaffold} \n" >&2;;
    -h  | --help ) Help ; exit 2 ;;
   esac
   shift
done 

if [ -z "${haplo1}" ] && [ -z "${haplo2}" ]  && [ -z "${scaffold}" ]    ; then
	Help
	exit 2
fi

#test if ancestral genome is provided or not:

#------------------------------ step 1 prepare input files  -------------------------------------#
#haplo1=$1 	#name of haplo1 
#haplo2=$2 	#name of haplo2 
#scaffold=$3 	#list of scaffold on the ancestral haplo 
#ancestral_genome=$4 #ancestral genome 

#cds file:
cdsfile1=$haplo1/08_best_run/$haplo1.spliced_cds.fa
cdsfile2=$haplo2/08_best_run/$haplo2.spliced_cds.fa


#remove the CDS length info that is introduced by gffread:
sed -i 's/ CDS=.*$//g' $cdsfile1
sed -i 's/ CDS=.*$//g' $cdsfile2

#-- get single copy orthologs from orthofinder ---
#-- criteria: we want 1:1:1 orthologs between the ancestralhaplo:haplo1:haplo2
#single copy path:

scopy=$(echo "genespace/orthofinder/Results_*/Orthogroups/Orthogroups_SingleCopyOrthologues.txt" ) 

if [ -n "$ancestral_sp" ] ; then
    echo "using ancestral genome"
    ancestral_vs_hap1=$(echo "genespace/orthofinder/Results_*/Orthologues/Orthologues_"$ancestral_genome"/"$ancestral_genome"__v__"$haplo1".tsv ")
    ancestral_vs_hap2=$(echo "genespace/orthofinder/Results_*/Orthologues/Orthologues_"$ancestral_genome"/"$ancestral_genome"__v__"$haplo2".tsv ")

    paste <(grep -Ff "$(echo $scopy )" "$(echo $ancestral_vs_hap1 )" )  <(grep -Ff "$(echo $scopy )" "$(echo $ancestral_vs_hap2 )" )  |\
        grep -Ff $scaffold - |\
        awk '{ if ($1 == $4) { print $1"\t"$2"\t"$3"\t"$6; } else { print $0"\tdifference exitst -- error"; } }' > paml/single.copy.orthologs 
        
        sed -i -e "s/\r//g"  paml/single.copy.orthologs

        cut  -f3 paml/single.copy.orthologs > paml/sco.$haplo1.txt
        cut  -f4 paml/single.copy.orthologs > paml/sco.$haplo2.txt

else
    echo "no ancestral genome"
    hap1_vs_hap2=$(echo "genespace/orthofinder/Results_*/Orthologues/Orthologues_"$haplo1"/"$haplo1"__v__"$haplo2".tsv ")

    paste <(grep -Ff "$(echo $scopy )" "$(echo $hap1_vs_hap2)" ) |\
	    grep -f <(cut -f 2 $scaffold ) - > paml/single.copy.orthologs  
       
        sed -i -e "s/\r//g"  paml/single.copy.orthologs

        cut  -f2 paml/single.copy.orthologs > paml/sco.$haplo1.txt
        cut  -f3 paml/single.copy.orthologs > paml/sco.$haplo2.txt

fi

##---  linearise the cds file ---#
cat "${cdsfile2}" | \
	awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}'  > paml/"$haplo2".linearised.cds
cat "$cdsfile1" | \
	awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}'  > paml/"$haplo1".linearised.cds


##---- recover the wanted sequences in the CDS file #
cd paml

#all is run from paml folder now

rm sorted* 2>/dev/null
while read -r pattern ; 
do 
    grep -w -A1 "$pattern" "$haplo2".linearised.cds  >> sorted."$haplo2".wanted_cds.fa ;
done < sco.$haplo2.txt 
#
while read -r pattern ; 
do 
    grep -w -A1 "$pattern" "$haplo1".linearised.cds >> sorted."$haplo1".wanted_cds.fa ; 
done < sco.$haplo1.txt 
#

### ---- then run paml and dnds.sh  ----- #
#paml does not like long fasta name so we shorten them for the sake of computation 
#to do here

#------ part specific to PAML --------------------- #
fasta1=sorted."$haplo1".wanted_cds.fa  
fasta2=sorted."$haplo2".wanted_cds.fa

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

	cp ../config/yn00_template.ctl sequence_files/tmp.${line[0]}.vs.${line[1]}/

        cd sequence_files/tmp.${line[0]}.vs.${line[1]}/

	path=$(pwd)
	echo $path
	sed -i "s|PATH|$path|g" 	 yn00_template.ctl #"

	yn00 yn00_template.ctl  

        cd ../../
	awk '/\+\-/ && !/(dS|SE)/ {split(FILENAME, a, "."); print $(NF-2), $(NF), $(NF-5), $(NF-3),"'${line[0]}'","'${line[1]}'"}'  sequence_files/tmp.${line[0]}.vs.${line[1]}/out_yn00_orthogp >  sequence_files/tmp.${line[0]}.vs.${line[1]}/resultat_Yang_Nielsen_2000_method.orthogp.txt 


done < wanted_sequence 2>&1 |tee log.paml 
cat sequence_files/tmp.*/resultat_Yang_Nielsen_2000_method.orthogp.txt >> results_YN.txt
