#!/bin/bash

#script to extract fasta sequence given the gtf coordinates
#source code: http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread

gtf=$1
genome=$2
output=${gtf%.gtf}.cds.fa

gffread -w $output -g $genome $gtf 

#then convert also the file to its cds:
input=$output
transeq -sequence $input -outseq ${input%.fa}.prot

exit
exemple:
gffread -w M.superbum.A1.combined.renamed.cds.fa -g ../../MvDp-1065-A1.contigs.fasta.masked.masked.masked.masked.simpl.fa  M.superbum.A1.combined.renamed.gtf 

