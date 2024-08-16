#!/bin/bash

fasta=$1
output=${fasta%.fa**}.lin.fasta

awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' "$fasta" > "$output"

