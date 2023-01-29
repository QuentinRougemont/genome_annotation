#!/bin/bash

#split gff by chromosomes
gff=$1
#gff=05-tsebra_results/M.superbumA1.combined.renamed.gff3 

mkdir gff_files 2>/dev/null
#cut -f 1 $gff |sort |uniq > list.chromosomes

for chr in $(cat list.chromosomes ) ; do
	grep $chr $gff > gff_files/"$chr".gff3
	echo -e "$chr" > gff_files/"$chr".scaffIDwithCDS
done

sed -i 's/.gff3//g' gff_files/*CDS

