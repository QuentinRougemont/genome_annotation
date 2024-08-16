#!/bin/bash

#split gff by chromosomes
gff=$1

mkdir gff_files 2>/dev/null
cut -f 1 "$gff" |sort |uniq > list.chromosomes

#for chr in $(cat list.chromosomes ) ; do
while read -r line ; do
	grep "$line" "$gff" > gff_files/"$line".gff3
	echo -e "$line" > gff_files/"$line".scaffIDwithCDS
done < list.chromosomes

sed -i 's/.gff3//g' gff_files/*CDS

