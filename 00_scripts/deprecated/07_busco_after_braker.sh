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

#activate busco
eval "$(conda shell.bash hook)"
conda activate busco570

run_busco=$(echo -e "busco -c8 -o busco_augustus -i "$input_fa" -l ""$lineage"" -m protein -f ")

if [ -z "$RNAseq" ] ; then 
  #RNAseq=NO
  echo -e "\tWARNING: No RNAseq info provided\nassuming no RNAseq data were used" 
  
   #run for the database:
   cd 06_braker/ 
   for i in round* ; do
        echo -e "----- running busco on: $i ------" 
        cd "$i"
        if [[ -d busco_augustus ]]
        then
            echo -e "WARNING directory busco_augustus already exists! check its content first
            \t Do you wish to remove it?\n
            \t if Yes this will rerun  busco computation\n
            \t if No this will skip busco computation \n"
            select yn in "Yes" "No"; do
                    case $yn in
                    Yes ) rm -rf; $run_busco ; break ;;
                    No  ) break ;;
                    esac
            done
        else
            bash $run_busco
        fi

        cd ../
   done 
fi


#run for rnaseq if rnaseq is existent:
if [[ "$RNAseq" = "YES" ]]
then
    echo "running busco on RNAseq data"
    cd 06_braker/rnaseq
    if [[ -d busco_augustus ]]
    then
        echo -e "WARNING directory busco_augustus already exists! check its content first
        \t Do you wish to remove it?\n
        \t if Yes this will rerun  busco computation\n
        \t if No this will skip busco computation \n"
        select yn in "Yes" "No"; do
                case $yn in
                Yes ) rm -rf; $run_busco ; break ;;
                No  ) echo "skipping busco" ; break ;;
                esac
        done
    else
        bash $run_busco
    fi
    cd ../

   #run for the database:
   for i in round* ; do
        echo -e "----- running busco on: $i ------" 
        cd "$i"
        if [[ -d busco_augustus ]]
        then
            echo -e "WARNING directory busco_augustus already exists! check its content first
            \t Do you wish to remove it?\n
            \t if Yes this will rerun  busco computation\n
            \t if No this will skip busco computation \n"
            select yn in "Yes" "No"; do
                    case $yn in
                    Yes ) rm -rf; $run_busco ; break ;;
                    No  ) echo "skipping busco" ;break ;;
                    esac
            done
        else
            bash $run_busco
        fi

        cd ../
   done 

fi
