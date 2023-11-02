#!/bin/bash

genome=$1 #name of the genome
new_name=$2 #name for the scaffold ids

sed "s/^>/>"$new_name"_/g" $genome > "$new_name".fa
exit
awk -v var="$new_name"  '{gsub(">[0-9A-Za-z-]*",">" var ,$1); print  }' $genome  > new_genome.fa 
