#!/bin/bash
set -e

current_command=$BASH_COMMAND
last_command=""

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

#purpose: run busco after braker in order to assess quality of the results
lineage=$1 
RNAseq=$2  #YES/NO

#------------- CONDA ACTIVATION  -------------- #
eval "$(conda shell.bash hook)"
conda activate superannot

#--- test braker version-------- 
# this will allow capturing the good amino-acid file for busco test:
#in braker2 the file name is "augustus.hints.aa"
#in braker3 the file name is "braker.aa"

v=$(braker.pl --version|awk '{print $3}' |cut -d "." -f 1 )
if [[ $v -ge 3 ]] ; 
then
    input_fa="braker.aa" ; 
else 
    input_fa="augustus.hints.aa" ; 
fi

echo braker input is $input_fa 


#activate busco
eval "$(conda shell.bash hook)"
conda activate busco571

run_busco=$(echo -e "busco -c8 -o busco_augustus -i "$input_fa" -l "$lineage" -m protein" -f)

if [ -z "$RNAseq" ] ; then 
  #RNAseq=NO
  echo -e "\tWARNING: No RNAseq info provided\nassuming no RNAseq data were used" 
  
   #run for the database:
   cd 06_braker/ 
   for i in round* ; do
        echo -e "----- running busco on: $i ------" 
        cd "$i"
        for file in busco_augustus/short_summary*txt
         do
              if [ ! -s "${file}" ]
              then 
                  $run_busco
              else
                  echo "busco already ok"
              fi
         done
         cd ../
   done
fi


#run for rnaseq if rnaseq is existent:
if [[ "$RNAseq" = "YES" ]]
then
    echo "running busco on RNAseq data"
    cd 06_braker/rnaseq
    for file in busco_augustus/short_summary*txt
    do
        if [ ! -s "${file}" ]
        then
            $run_busco
        else
            echo "busco already ok"
        fi
    cd ../
    done

   #run for the database:
   for i in round* ; do
        echo -e "----- running busco on: $i ------" 
        cd "$i"
        for file in busco_augustus/short_summary*txt
        do
            if [ ! -s "${file}" ]
            then
                $run_busco
            else
                echo "busco already ok"
            fi
        done
        cd ../
   done 

fi
