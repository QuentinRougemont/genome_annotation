#!/bin/bash

#script to extract the longest transcript and evaluate the final quality with busco
#Note that AGAT can do the same probably better, but this does not require any supplementary installation.


cd 06_braker
best_round=$(grep "C:" round*_braker_on_refprot/busco_augustus/short_summary.specific.basidiomycota_odb10.busco_augustus.txt |\
	sed -e 's/%//g' -e 's/\[/,/g' -e 's/]//g' -e 's/:/\t/g' -e  's/,/\t/g' |\
	LC_ALL=C sort -nr -k3 -k5 -n -k7 -k9 -k11  |sed -n 1p |cut -d "/" -f 1 ) 

echo best round is "$best_round"
 
cd $best_round


samtools faidx augustus.hints.aa
awk -F "." '{print $1"\t"$0}' augustus.hints.aa.fai |\
	awk '$3>max[$1]{max[$1]=$3; row[$1]=$2} END{for (i in row) print row[i]}' > longest.transcript


awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' augustus.hints.aa > augustus.hints.aa.lin.fasta 

grep -A1 -Ff longest.transcript augustus.hints.aa.lin.fasta > longest_transcript.fa

cp longest_transcript.fa ../../

#then run busco
eval "$(conda shell.bash hook)"
conda activate busco_env #activate conda
cd ../../
00_scripts/utility_scripts/busco.sh longest_transcript.fa 

