#!/bin/bash

#script to compute the length of a fasta
fasta=$1
awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' $fasta >> fasta_length.txt
