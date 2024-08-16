#!/bin/bash
#miniscript to run tsebra
#Author: QR
#Date: 23-01-23

if [ $# -ne 2  ]; then
    echo "USAGE: $0 species_name best_database_run "
    echo "Expecting a basename for the species and the name of the best run 
    from braker using database from protein "
    exit 1
else
    name=$1
    best_run=$2
    echo "Name is : ${name}"
    echo " "
    echo "best run is ${best_run}"
    echo "*****"
fi

run1=06_braker/rnaseq
run2=06_braker/$best_run 

b1=$run1/braker.gtf
b2=$run2/braker.gtf


new_b1=$run1/braker1.gtf
new_b2=$run2/braker2.gtf

# ----- fix braker gtf ids ----- #
fix_gtf_ids.py --gtf "$b1" --out "$new_b1"
fix_gtf_ids.py --gtf "$b2" --out "$new_b2"

b1="$new_b1"
b2="$new_b2"

rm -rf 07-tsebra_results/ 2>/dev/null
mkdir 07-tsebra_results/
cp default.cfg 07-tsebra_results/

# ----- run TSEBRA ----- #
tsebra.py -g "${b1}","${b2}" \
          -c default.cfg \
          -e  "${run1}"/hintsfile.gff,"${run2}"/hintsfile.gff\
          -o 07-tsebra_results/"$name".combined.gtf

