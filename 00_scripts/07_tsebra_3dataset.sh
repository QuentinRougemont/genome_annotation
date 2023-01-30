#!/bin/bash
#miniscript to run tsebra
#Author: QR
#Date: 23-01-23

if [ $# -ne 3  ]; then
    echo "USAGE: $0 species_name best_database_run best_database2_run "
    echo "Expecting a basename for the species and the name of the best  run from braker using database from protein "
    exit 1
else
    name=$1
    best_run=$2
    best_run_db2=$3
    echo "Name is : ${name}"
    echo " "
    echo "best run is ${best_run}"
    echo ""
    echo "best run for second database is ${best_run_db2}"
    echo "***************"
fi

run1=06_braker/rnaseq
run2=$best_run #round2_braker_on_refprot_2023-01-29_21h27m06s round3_braker_on_refprot_2023-01-17_17h34m03s
run3=$best_run_db2

b1=$run1/braker.gtf
b2=$run2/braker.gtf
b3=$run3/braker.gtf

new_b1=$run1/braker1.gtf
new_b2=$run2/braker2.gtf
new_b3=$run3/braker3.gtf


# ----- fix braker gtf ids ----- #
fix_gtf_ids.py --gtf $b1 --out $new_b1
fix_gtf_ids.py --gtf $b2 --out $new_b2
fix_gtf_ids.py --gtf $b3 --out $new_b3

b1=$new_b1
b2=$new_b2
b3=$new_b3

rm -rf 07-tsebra_results_3dataset/ 2>/dev/null
mkdir  07-tsebra_results_3dataset/
cp default.cfg 07-tsebra_results_3dataset/

# ----- run TSEBRA ----- #
tsebra.py -g "${b1}","${b2}","${b3}" \
          -c default.cfg \
          -e  "${run1}"/hintsfile.gff,"${run2}"/hintsfile.gff,"${run3}"/hintsfile.gff\
          -o 07-tsebra_results_3dataset/$name.combined.gtf

