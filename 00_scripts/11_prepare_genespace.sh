#!/bin/bash

#Purpose:
#master script to prepare bed files, laucnch GeneSpace, run paml and launch downstream Rscript 
#will check existence of all dependencies
#Date: 2023
#Author: QR

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo -e "master script to: \n 1 - create bed files, \n 2 - launch GeneSpace, \n 3 - run paml and launch downstream Rscripts (Rideogram, plot paml, etc)"
   echo " "
   echo "Usage: $0 [-s1|-s2|-f|-a|-h|]"
   echo "options:"
   echo " -h|--help: Print this Help."
   echo " -s1|--species1: the name of the first  focal species\t "
   echo " -s2|--species2: the name of the second focal species\t "
   echo " -a |--ancestral_sp: the name of the ancestral species to infer orthology and plot gene order"
   echo " -f|--folderpath: the path to the global folder containing species1 and species 2"
   echo " "
   echo "dependancies: orthofinder, mcscanx, GeneSpace, paml (yn00), Rideogram, translatorX "
}


###########################################################
## to do: add more support to handle the ancestral species: it could be either only the genome + gtf or directly a bed + protein file for instance
## test if MCScanX is install
## test if orthoFinder is installed
## test if paml/yn00 is installed
##########################################################


############################################################
# Process the input options.                               #
############################################################
while [ $# -gt 0 ] ; do
  case $1 in
    -s1 | --species1) species1="$2" ; echo -e "Species 1 Name is ***${species1}*** \n" >&2;;
    -s2 | --species2) species2="$2" ; echo -e "Species 1 Name is ***${species1}*** \n" >&2;;
    -a  | --ancestral_sp) ancestral_sp="$2" ; echo -e "ancestral species  Name is ***${ancestral_sp}*** \n" >&2;;
    -f  | --folderpath  ) folderpath="$2"  ; echo -e "global folder is  ${folderpath} \n" >&2;;
    -h  | --help ) Help ; exit 2 ;;
   esac
   shift
done 

if [ -z "${species1}" ] || [ -z "${species2}" ] || [ -z "${folderpath}" ]  || [ -z "${ancestral_sp}" ]    ; then
	Help
	exit 2
fi

#------------------------------ step 1 prepare bed file for each species -------------------------------------#
#really simple:

cd $folderpath

#remove any existing folder:
rm genespace peptide paml plots -rf 2>/dev/null
#to do: insert alert message for the user.
mkdir -p genespace/bed genespace/peptide paml plots

# create bed
awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$10}' $species1/08_best_run/$species1.longest_transcript.gtf |sed 's/"//g' > genespace/bed/$species1.bed
awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$10}' $species2/08_best_run/$species2.longest_transcript.gtf |sed 's/"//g' > genespace/bed/$species2.bed

# simplify the protein file to match the bed (i.e. remove the _1 inserted by transeq and the CDS length info):
sed 's/_1 CDS=.*$//g' $species1/08_best_run/"$species1"_prot.fa > genespace/peptide/$species1.fa
sed 's/_1 CDS=.*$//g' $species2/08_best_run/"$species2"_prot.fa > genespace/peptide/$species2.fa

#verify that IDs in bed and fasta file are matching - else exit  
grep ">" genespace/peptide/$species1.fa |sed 's/>//g' > tmp1
grep ">" genespace/peptide/$species2.fa |sed 's/>//g' > tmp2

check1=$(grep -Ff tmp1 genespace/bed/$species1.bed |wc -l )
check2=$(grep -Ff tmp2 genespace/bed/$species2.bed |wc -l )

echo -e "check2 size is $check2"
echo -e "check1 size is $check1"

bedsize1=$(wc -l genespace/bed/$species1.bed |awk '{print $1}' )
bedsize2=$(wc -l genespace/bed/$species2.bed |awk '{print $1}' )

echo -e "bedisze2  size is $bedsize2"
echo -e "bedisze1  size is $bedsize1"

#check that all is matching:
if [ "$bedsize1" = "$check1" ]
then
	echo "input1 is ok" 
	rm tmp1
else
	echo "input1 is not ok"
	echo "check your data"
	exit 2
fi

if [ "$bedsize2" = "$check2" ]
then
	echo "input2 is ok" 
	rm tmp2
else
	echo "input2 is not ok"
	echo "check your data"
	exit 2
fi

# -- handling ancestral species ------
# -- this part assumes that a nice bed and peptide file are existant for the ancestral species
# -- here we used a genome annotated with the same pipeline relying on braker 

cd genespace/bed/
ln -s ../../../"$ancestral_sp"/"$ancestral_sp".bed . 
cd ../peptide
ln -s ../../../"$ancestral_sp"/"$ancestral_sp".prot.fa "$ancestral_sp".fa

cd ../../

#------------------------------ step 2 run GeneSpace ---------------------------------------------------------#
cd genespace 
Rscript ../../00_scripts/Rscripts/01.run_geneSpace.R
Rscript ../../00_scripts/Rscripts/02.plot_geneSpace.R
#plot genespace subspace of target chromosomes: 
cd ../
#------------------------------ step 3 run paml  -------------------------------------------------------------#

scaffold=scaffold.txt #hardcoded variable to be passed as an argument for later
cp ../$scaffold .
echo $(pwd)
echo species1 is "$species1"
echo species2 is "$species2"

../00_scripts/12_command_line.paml.sh "$species1" "$species2" $scaffold 

pamlsize=$(wc -l paml/results_YN.txt |awk '{print $1}' ) 
scpo=$(wc -l paml/single.copy.orthologs |awk '{print $1}' )

echo -e "there is $pamlsize results for PAML \n"
echo -e "there is $scpo single copy orthologs \n" 


#-- step4 -- plot paml results 
#test if previous step was successfull else plot or exit with high levels of pain
echo $ancestral_sp 

Rscript ../00_scripts/Rscripts/03.plot_paml.R $ancestral_sp $species1 $species2


#-- step5 -- plot ideogram 
#test if previous step was successfull else plot or exit with high levels of pain
samtools faidx $species1/03_genome/"$species1".fa
samtools faidx $species2/03_genome/"$species2".fa

Rscript ../00_scripts/Rscripts/04.ideogram.R $species1 $species2
 
#-- step6 -- run minimap between the genomes 
#run minimap on the genome 
#assumption : each genome MUST BE located in folder 03-genome
