#!/bin/bash

#script to split fasta
#warning it must be linearized!!

fasta=$1

awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' "$fasta" > tmp

mkdir fasta_files

i=1;
while read -r line ; do
  if [ ${line:0:1} == ">" ] ; then
    filename=$(echo "$line" | cut -d ":" -f1 | tr -d ">")
    touch fasta_files/"$filename".fasta 
    echo "$line" >> fasta_files/"${filename}".fasta
  else
    echo "$line" >> fasta_files/"${filename}".fasta
    ((i++))
  fi
done < tmp

#cleanup
rm tmp

