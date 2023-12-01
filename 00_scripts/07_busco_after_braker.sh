#!/bin/bash

#purpose: run busco after braker in order to assess quality of the results
lineage=$1 #~/basidiomycota_odb10 #to be passed as an argument #
RNASeq=$2  #YES/NO

#activate busco:
eval "$(conda shell.bash hook)"
conda activate busco_env


if [ -z "$RNAseq" ] ; then 
  #RNAseq=NO
  echo -e "WARNING: No RNAseq info provided\n
           assuming no RNAseq data were used" 
  
   #run for the database:
   cd 06_braker/ 
   input_fa=augustus.hints.aa 
   for i in round* ; do
       cd $i
       busco -c8 -o busco_augustus -i $input_fa -l $lineage -m protein #--updata-data #to update the database if there's a warning
       cd ../
   done 
fi


#the braker2 results is called "augustus.hints.aa"
input_fa=augustus.hints.aa 

#run for rnaseq if rnaseq is existent:
if [[ $RNAseq = "YES" ]]
then
   echo "running busco on RNAseq data"
   cd 06_braker/rnaseq
   busco -c8 -o busco_augustus -i $input_fa -l $lineage -m protein #--updata-data #to update the database if there's a warning

    #run for the database:
    cd ../
    for i in round* ; do
        cd $i
        busco -c8 -o busco_augustus -i $input_fa -l $lineage -m protein #--updata-data #to update the database if there's a warning
        cd ../
    done 
fi

