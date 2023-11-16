#!/bin/bash

#purpose: run busco after braker in order to assess quality of the results
lineage=~/basidiomycota_odb10 #to be passed as an argument #

#source /local/env/envbusco-5.2.2.sh 
#activate busco:
eval "$(conda shell.bash hook)"
conda activate busco_env

RNASeq=$1

if [ -z "$RNAseq" ] ; then 
  RNAseq=NO
  echo -e "WARNING: No RNAseq info provided\n
           assuming no RNAseq data were used" 
fi


#the braker2 results is called "augustus.hints.aa"
input_fa=augustus.hints.aa 

#run for rnaseq if rnaseq is existent:
if [[ $RNAseq = "YES" ]]
then
   echo "running busco on RNAseq data"
   cd 06_braker/rnaseq
   busco -c8 -o busco_augustus -i $input_fa -l ~/basidiomycota_odb10 -m protein #--updata-data #to update the database if there's a warning
fi

#run for the database:
cd ../
for i in round* ; do
    cd $i
    busco -c8 -o busco_augustus -i $input_fa -l ~/basidiomycota_odb10 -m protein #--updata-data #to update the database if there's a warning
    cd ../
done 
