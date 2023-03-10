#!/bin/bash
#Date:13-02-23
#Author QR
# WARNING STILL BETA #! 

##  ------------------------ general parameters --------------------------------  ##
while [ $# -gt 0 ] ; do
  case $1 in
    -g | --genome) genome="$2" ;echo "the genome file  is: $genome" >&2;;
    -s | --species) species="$2" ;echo "the species name will be $species" >&2;;
    -m | --mask )   Mask="$2" ; echo "unknown TE will be removed after repeatemasking" >&2;;
    -r | --rna )    RNAseq="$2" ; echo "No RNAseq provided" >&2;; 
    -h | --help) echo -e "Option required:
    -g/--genome \t the reference genome file 
    -s/--species\t the species name (used for database building and basename in several steps)
    Optional:
    -m/--mask\t a string stating wether unknown TE should be removed for genome annotation (YES/NO) -- default: YES
    -r/--ref \t a string stating wether RNAseq data should be used (YES/NO) -- default: NO
    " >&2;exit 1;;
    esac
    shift
done

if [ -z "$genome" ] || [ -z "$species" ] ; then
	echo >&2 "Fatal error: Ref genome (-g), and species name (-s) not defined\n
	see manual with -h or --help"
exit 2
fi

#optional parameters:
if [ -z "$Mask" ] ; then
    Mask=YES
fi
if [ -z "$RNAseq" ] ; then 
  RNAseq=NO
fi

# ------------------ run RepeatModeler and RepeatMasker ------------------  ##
./00_scripts/05_repeatmodeler.sh "$genome" "$species" "$Mask" 2>&1 |tee log_rm

#here test if the genome should be mask or not, then test its existence
#if it does not exist then exit
#it it exist run braker with the appropriate parameter

# -------------------- run Braker  ---------------------------- #
./00_scripts/06_braker.sh 03_genome/"$species".genome.wholemask_no_unknown.fa $species $RNAseq 2>&1 |tee braker_log  #NO for no rnaseq  


# -------------------- run Busco  ---------------------------- #
eval "$(conda shell.bash hook)"
conda activate busco_env #activate conda

cd 06_braker
for i in round* ; do cd $i ;  ../00_scripts/utility_scripts/busco.sh augustus.hints.aa ; cd ../ ; done

if [[ $RNAseq = "YES" ]] ; then
    cd rnaseq ;
    ../../00_scripts/utility_scripts/busco.sh augustus.hints.aa 

    #also generate a report here using braker scripts:
    ../../00_scripts/utility_scripts/GenerateReport.py augustus.hints.gtf hintsfile.gff reportRNAseq.$species.pdf
    cd ../
fi
conda deactivate busco_env #activate conda

#choose the best round
best_round=$(grep "C:" round*/busco*/short*txt |\
	sed -e 's/%//g' -e 's/\[/,/g' -e 's/]//g' -e 's/:/\t/g' -e  's/,/\t/g' |\
	LC_ALL=C sort -nr -k3 -k5 -n -k7 -k9 -k11  |sed -n 1p |cut -d "/" -f 1 )

#note: you may wan to change the sorting to pick a better run

cd ../

if [[ $RNAseq = "NO" ]] ; then
	echo "work is finished "
	echo best_round is $best_round

	cd 06_braker/
	cd $best_round
	../../00_scripts/utility_scripts/GenerateReport.py augustus.hints.gtf hintsfile.gff report.proteinDB.$species.pdf
	exit 1
else
	echo "now combining RNAseq and Protein with TSEBRA"
	#to do
	#now combine tesbra best_round and RNAseq
	./00_scripts/07_tsebra.sh $species $best_round 
	
	#rename the gene id in the gtf:
	rename_gtf.py --gtf ./07-tsebra_results/${species}.combined.gtf --prefix $species --translation_tab 07-tsebra_results/translation_tab --out 07-tsebra_results/$species.renamed.gtf 


	## Fix TSEBRA output (source: https://github.com/Gaius-Augustus/BRAKER/issues/457 )
	# Fix lack of gene_id and transcript_id tags in gtf file column 9
	cat ./07-tsebra_results/${species}.renamed.gtf | ./00_scripts/utility_scripts/Fix_Augustus_gtf.pl > ./07-tsebra_results/${species}.renamed.fixed.gtf


	./00_scripts/08_extractcds.sh ./07-tsebra_results/${species}.renamed.fixed.gtf 03_genome/"$species".genome.wholemask_no_unknown.fa 
	
	eval "$(conda shell.bash hook)"
	conda activate busco_env #activate conda
	cd 07-tsebra_results
	#running busco
	 ../00_scripts/utility_scripts/busco.sh ${species}.renamed.fixed.cds.prot 
fi

exit 
#Further stuff that will be add:
#1 - evaluate final quality with busco
#2 - run several of the agat module to check quality
#3 - run interproscan

