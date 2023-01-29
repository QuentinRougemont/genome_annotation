#!/bin/bash

#activate busco:
#conda activate busco_env
#script to run busco
input_fa=$1 #augustus.hints.aa #fasta acid-aminée/protein/genome file

#for braker:
busco -c2 -o busco_augustus -i $input_fa -l ~/basidiomycota_odb10 -m protein #--updata-data #to update the database if there's a warning

exit

#base commande:
busco -c2 -o output_name -i input_fasta_file -l ~/basidiomycota_odb10 -m protein #--updata-data #to update the database if there's a warning

#users guide: https://busco.ezlab.org/busco_userguide.html#lineage-datasets

#use -m geno for genome
#use -m protein for protein
