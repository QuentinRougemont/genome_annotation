#!/bin/bash

#DATE: 19-10-23
#Author: QR
#Purpose: run miniprot and extract PR

#
############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo -e "script to run miniprot aligning a genome against a reference protein file.\n This allow to find important protein present in the reference genome but not predicted in the previous step\n in details this will allow to: \n 1 - allign the genome against prot \n 2 - search for the missing protein \n 3 - add it to the gtf file "
   echo " "
   echo "Usage: $0 [-g|p|r|f|]"
   echo "options:"
   echo "-h|--help: Print this Help."
   echo "-g|--genomes: the name of the reference genome for the target species missing the genes of interest"
   echo "-r|--pr    : the name of the gene of interest that is missing (e.g. name of PR gene)"
   echo "-p|--prot  : the name of the protein file from which the gene of interest (e.g. PR) will be extracted"
   echo "-f|--gff   : the gff file in which the gene of interest is missing"
   echo " "
   echo "dependancies: TSEBRA, samtools, gffread, transeq "
}


############################################################
# general parameters arguments                             #
############################################################

while [ $# -gt 0 ] ; do
  case $1 in
    -g | --genome) refgenome="$2" ;echo "the genome file of the target speices  is: $genome" >&2;;
    -p | --prot)   protfile="$2" ;echo "the protein file from the outgroup species containing PR with be $protfile" >&2;;
    -r | --pr )    PRinPROT="$2" ; echo "the PR gene in the protein file will be $PRinPROT" >&2;;
    -f | --gff )   gff_file="$2" ; echo "gff_file missing the PR will be $gff_file" >&2;; 
    -h | --help) Help ; exit 2;;
    esac
    shift
done

if [ -z "$refgenome" ] || [ -z "$protfile" ] || [ -z "$PRinPROT" ] || [ -z "$gff_file" ]; then
	Help	
	exit 2
fi

protbase=$(basename $protfile )
refbase=$(basename $refgenome )
base=$(basename $gff_file )

#---------------------step 0: ---------------------------------------------------

#you MUST have run orthofinder to find the closest relative species 
#this is the species (or set of species) we will use to input the missing PR gene

#-------------------- step1: run miniprot  --------------------------------------

miniprot --gff -I $refgenome $protfile > whole."$protbase"_on_"$refbase".gff


#-------------------- step2: extract PR  ----------------------------------------
mkdir 09_gff_final

#some note for external users:
#sed 1d => remove the PAF info
#sed 1p to duplicate the mRNA colomn and transform into "transcript"
#sed -e .... is just removing unwanted stuff, we could have reach the same in awk I guess
#then awk to the magic renaming :) 
#Note we will add a "fake" 'transcript type ".t1" into the column $9 to closely follow AUGUSTUS prediction from braker  
grep -A1 "$PRinPROT" whole."$protbase"_on_"$refbase".gff |\
	sed 1d |grep -v "#PAF" |\
	sed '1p' |\
	sed -e '1s/mRNA/gene/' -e '2s/mRNA/transcript/' -e 's/ID=//' -e 's/Parent=//'  -e 's/;//' -e 's/Rank=.*//'  |\
	awk '{ 
	if ($3=="gene") 
		{ for(i=1; i<=NF-1; i++) printf $i"\t"; print "gene_id \"" $1 "_" $9 "\"" } 
	else if ($3 == "transcript" ) 
		{ for(i=1; i<=NF-1; i++) printf $i"\t"; print "transcript_id \""$1 "_" $9 ".t1\"" } 
	else { for(i=1; i<=NF-1; i++) printf $i"\t"; print "transcript_id \""$1 "_" $9 ".t1\"; gene_id \"" $1 "_" $9 "\"" } }'  > PR.${refbase%.fa}.gff


#here test if the results was sucessul else exit with error 1

if [ ! -s "PR.${refbase%.fa}.gff" ]; then
    echo "Error: File is empty"
    echo "minimap did not align the PR on your genome - ckeck parameters"
    exit 1
else
	echo "base is $base"
	cat $gff_file PR.${refbase%.fa}.gff > 09_gff_final/"$base"
	
	#-------------------------------- step 3 ---------------------------------------"
	#extract spliced CDS:
	gffread -g $refgenome -w 09_gff_final/${base%.gtf}_cds.fa 09_gff_final/$base

	#extract Proteins!!
	gffread -g $refgenome -y 09_gff_final/${base%.gtf}_prot.fa 09_gff_final/$base
	
fi
