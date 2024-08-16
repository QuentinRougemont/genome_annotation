#!/bin/bash
#script to run busco

input_fa=$1 #ex: augustus.hints.aa 


busco -c2 -o busco_augustus -i "$input_fa" -l basidiomycota_odb10 -m protein #--update-data #to update the database if there's a warning

exit
#use -m geno for genome
#use -m protein for protein
