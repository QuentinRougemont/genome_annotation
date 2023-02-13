#!/bin/bash
#Date:13-02-23
#Author QR

genome=$1 #
species=$2
mask=YES
RNAseq=NO

# ------------------ run RepeatModeler and RepeatMasker ------------------  ##
./00_scripts/05_repeatmodeler.sh "$genome" "$species" "$mask" 2>&1 |tee log_rm

#here test if the genome should be mask or not, then test its existence
#if it does not exist then exit
#it it exist run braker with the appropriate parameter

# -------------------- run Braker  ---------------------------- #
./00_scripts/06_braker.sh 03_genome/genome.wholemask_no_unknown.fa $species $RNAseq 2>&1 |tee braker_log  #NO for no rnaseq  


# -------------------- run Busco  ---------------------------- #
eval "$(conda shell.bash hook)"
conda activate busco_env #activate conda

cd 06_braker
for i in round* ; do cd $i ;  ~/busco.sh augustus.hints.aa ; cd ../ ; done

if [[ $RNAseq = YES]] ; then
    cd rnaseq ;
    ~/busco.sh augustus.hints.aa 
    cd ../
fi

#choose the best round
best_round=$(grep "C:" round*/busco*/short*txt |\
	sed -e 's/%//g' -e 's/\[/,/g' -e 's/]//g' -e 's/:/\t/g' -e  's/,/\t/g' |\
	LC_ALL=C sort -nr -k3 -k5 -k7 -k9 -k11 |sed -n 1p |cut -f 1 )

cd ../

if [[ $NRAseq = "NO" ]] ; then
	echo "work is finished "
	echo best_round is $best_round
else
	echo "now combining RNAseq and Protein with TSEBRA"
	#to do
	#now combine tesbra best_round and RNAseq
	./00_scripts/07_tsebra.sh $species $best_round 
fi

## Fix TSEBRA output (source: https://github.com/Gaius-Augustus/BRAKER/issues/457 )
# Fix lack of gene_id and transcript_id tags in gtf file column 9
cat ./07_tsebra_results/${species}.combined.gtf | ./00_scripts/Fix_Augustus_gtf.pl > ./07_tsebra_results/${species}.combined.fixed.gtf


./00_scripts/08_extractcds.sh ./07_tsebra_results/${species}.combined.fixed.gtf 03_genome/genome.wholemask_no_unknown.fa 

exit 
#Further stuff that is done:
#1 - evaluate final quality with busco
#2 - run several of the agat module to check quality
#3 - run interproscan

