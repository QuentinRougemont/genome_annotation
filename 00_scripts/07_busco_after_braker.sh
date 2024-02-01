#!/bin/bash

#purpose: run busco after braker in order to assess quality of the results
lineage=$1 
RNAseq=$2  #YES/NO

#activate braker 
eval "$(conda shell.bash hook)"
conda activate braker_env

#--- test braker version-------- 
# this will allow capturing the good amino-acid file for busco test:
#in braker2 the file name is "augustus.hints.aa"
#in braker3 the file name is "braker.aa"

v=$(braker.pl --version|awk '{print $3}' )
if [[ $v > 3 ]] ; then input_fa="braker.aa" ; else input_fa="augustus.hints.aa" ; fi


#activate busco
eval "$(conda shell.bash hook)"
conda activate busco_env

if [ -z "$RNAseq" ] ; then 
  #RNAseq=NO
  echo -e "\tWARNING: No RNAseq info provided\nassuming no RNAseq data were used" 
  
   #run for the database:
   cd 06_braker/ 
   for i in round* ; do
	   echo "----- running busco on: $i ------" 
       cd $i
       busco -c8 -o busco_augustus -i $input_fa -l $lineage -m protein #--updata-data #to update the database if there's a warning
       cd ../
   done 
fi


#run for rnaseq if rnaseq is existent:
if [[ $RNAseq = "YES" ]]
then
   echo "running busco on RNAseq data"
   cd 06_braker/rnaseq
   busco -c8 -o busco_augustus -i $input_fa -l $lineage -m protein #--updata-data #to update the database if there's a warning
   cd ../../
fi

