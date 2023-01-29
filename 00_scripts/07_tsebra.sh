#!/bin/bash

#miniscript to run tsebra

name=$1 #name of the species

run1=06_braker/rnaseq
run2=06_braker/round3_braker_on_refprot_2023-01-17_17h34m03s

b1=$run1/braker.gtf
b2=$run2/braker.gtf


new_b1=$run1/braker1.gtf
new_b2=$run2/braker2.gtf

# ----- fix braker gtf ids ----- #
fix_gtf_ids.py --gtf $b1 --out $new_b1
fix_gtf_ids.py --gtf $b2 --out $new_b2

b1=$new_b1
b2=$new_b2

mkdir 07-tsebra_results/
cp default.cfg 07-tsebra_results/

# ----- run TSEBRA ----- #
tsebra.py -g "${b1}","${b2}" \
	  -c default.cfg \
	  -e  "${run1}"/hintsfile.gff,"${run2}"/hintsfile.gff\
	  -o 07-tsebra_results/$name.combined.gtf 

exit

# ----- using AUGSUSTUS hints instead ----- #  
tsebra.py -g "${run1}"/augustus.hints.gtf,"${run2}"/augustus.hints.gtf \
	  -c default.cfg \
	  -e  "${run1}"/hintsfile.gff,"${run2}"/hintsfile.gff\
	  -o $name.combined.gtf 
