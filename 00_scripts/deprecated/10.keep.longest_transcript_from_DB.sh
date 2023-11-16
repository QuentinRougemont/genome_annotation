#!/bin/bash
#Date: 10-03-2023
#Author: QR
#Purpose:
#script to extract the longest transcript and evaluate the final quality with busco
#Note that AGAT can do the same probably better, but this does not require any supplementary installation.

if (( $# < 1 )) ; then
	echo "USAGE: $0 data_type  (optional :protein file name)" 
        echo -e "Expecting the following parameters:\n
	1 - the type of data: (i) either TSEBRA or (ii)  DB\n
        2 - if data type == TSBEBRA then a name for the protein file should be provided\n"
        exit 1
else
    data_type=$1 #the reference genome here! "genome.wholemask.fa"
    protein_file_name=$2
fi


if [[ -z $data_type ]] ; then
	echo >&2 "Fatal error: please provide the type of data \n
	this must be DB or TSEBRA \n
	if TSEBRA please also provid tha name of the protein file as reconstructed from script 08"
	exit 2
fi

if [[ $data_type == "DB" ]] ; then

	echo "assuming data derived from protein database only"
	echo "will extract longest transcript from braker output 'augustus_hints.aa' file "

	cd 06_braker
	#capturing best run of busco from the previous script 
	best_round=$(grep "C:" round*_braker_on_refprot/busco_augustus/short_summary.specific.basidiomycota_odb10.busco_augustus.txt |\
		sed -e 's/%//g' -e 's/\[/,/g' -e 's/]//g' -e 's/:/\t/g' -e  's/,/\t/g' |\
		LC_ALL=C sort -nr -k3 -k5 -n -k7 -k9 -k11  |sed -n 1p |cut -d "/" -f 1 ) 

	echo best round is "$best_round"
 
	cd $best_round
		
	samtools faidx augustus.hints.aa
	#extract longest transcript for the best run: 
	awk -F "." '{print $1"\t"$0}' augustus.hints.aa.fai |\
		awk '$3>max[$1]{max[$1]=$3; row[$1]=$2} END{for (i in row) print row[i]}' > longest.transcript

	
	#linearize file so that the next command will work:
	awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' augustus.hints.aa > augustus.hints.aa.lin.fasta 

	grep -A1 -Ff longest.transcript augustus.hints.aa.lin.fasta > longest_transcript.fa

	cp longest_transcript.fa ../../
	cd ../../

else

	echo "assuming data type are derived from TSEBRA (RNAseq + DB type) "
	echo "protein file is $protein_file_name"

	if [[ -z $protein_file_name ]] ; then
		echo "error no protein_file_name provided this can't work"
		echo "please provide a dataset name"	
		exit 1
	else

	cd 07-tsebra_results

	samtools faidx $protein_file_name
	awk -F "." '{print $1"\t"$0}' $protein_file_name.fai |\
	awk '$3>max[$1]{max[$1]=$3; row[$1]=$2} END{for (i in row) print row[i]}' > longest.transcript

	awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' $protein_file_name > "$protein_file_name".lin.fasta

	grep -A1 -Ff longest.transcript "$protein_file_name".lin.fasta  > longest_transcript.fa

	cp longest_transcript.fa ../
        cd ../

        fi
   
fi

#then run busco
eval "$(conda shell.bash hook)"
conda activate busco_env #activate conda
00_scripts/utility_scripts/busco.sh longest_transcript.fa 

