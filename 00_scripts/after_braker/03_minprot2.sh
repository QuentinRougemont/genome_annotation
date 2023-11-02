#!/bin/bash
#DATE: 19-10-23
#Author: QR
#Purpose: run miniprot and extract PR
#optionally: add PR to the existing filtered gtf and extract the corresponding protein so we get a cleaned dataset! 


#to do:
#make a parser and test if all argument have been provided, otherwise exti
##  ------------------------ general parameters --------------------------------  ##
while [ $# -gt 0 ] ; do
  case $1 in
    -g | --genome) refgenome="$2" ;echo "the genome file of the target speices  is: $genome" >&2;;
    -p | --prot)   protfile="$2" ;echo "the protein file from the outgroup species containing PR with be $protfile" >&2;;
    -r | --pr )    PRinPROT="$2" ; echo "the PR gene in the protein file will be $PRinPROT" >&2;;
    -f | --gff )   gff_file="$2" ; echo "gff_file missing the PR will be $gff_file" >&2;; 
    -h | --help) echo -e "Option required:
    -g/--genome \t the reference genome of the species without PR
    -p/--prot \t the name of the protein file from which PR will be extracted
    -r/--pr \t   the name of PR
    -f/--gff \t  the gff file of the species lacking PR
    #optional: 
    #-s/--species \t the species name (used for database building and basename in several steps)
    #-r/--ref \t a string stating wether RNAseq data should be used (YES/NO) -- default: NO
    " >&2;exit 1;;
    esac
    shift
done

#protfile=$1  #Msco-A1.fa #proteinfile
#refgenome=$2 #Msco-A2.fa
#PRinPROT=$3  #Mscorzo-A1_tig00000006_g6796.t1
#gff_file=$4  #the gff file on which we want to add the PR! 
if [ -z "$refgenome" ] || [ -z "$protfile" ] || [ -z "$PRinPROT" ] || [ -z "$gff_file" ]; then
	echo >&2 "Fatal error: Ref genome (-g), proteinfile (-p) PRgene (-r) and gff file (-f) not defined\n
	see manual with -h or --help"
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
	#make a backup of the original alignment to PR:
	grep -A1 "$PRinPROT" whole."$protbase"_on_"$refbase".gff > PR.miniprot.bkp
	cat $gff_file PR.${refbase%.fa}.gff > 09_gff_final/"$base"
	
	#-------------------------------- step 3 ---------------------------------------"
	#note: here we wante to run the ancient code to remove duplicate genes etc ----!
	# 		TO DO
	#
	#

	#--------------------------------step 4 ---------------------------------------#
	#extract spliced CDS:
	gffread -g $refgenome -w 09_gff_final/${base%.gtf}_cds.fa 09_gff_final/$base

	#extract Proteins!!
	gffread -g $refgenome -y 09_gff_final/${base%.gtf}_prot.fa 09_gff_final/$base

	eval "$(conda shell.bash hook)"
	conda activate busco_env #cluster
	#conda activate env_busco #on me
	#script to run busco
	cd 09_gff_final
	input_fa=${base%.gtf}_prot.fa 

	#for braker:
	busco -c10 -o busco_"$input_fa" -i "$input_fa" -l basidiomycota_odb10 -m protein #
	#count the number of protein:
  	cd ../	

fi

exit
#note: here we wante to run the ancient code to remove duplicate genes etc ----!

#extract spliced CDS:
gffread -g $refgenome -w ${base%.gtf}_cds.fa new.$base

#extract Proteins!!
gffread -g $refgenome -y ${base%.gtf}_prot.fa new.$base
