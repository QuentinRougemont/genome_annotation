#!/bin/bash

#script to rename the scaffold in the gtf, rename the gtf, create a synchronised genome
#extract longest protein
#remove redundant GeneMark gene 

folder=$1   #the folder containing the species	
species=$2  #the new name of the species that we want to use

cd $folder

#----------------------------------------------------------------------#
#step2 get the best braker run
#capturing best run of busco from the previous script 
cd 06_braker
best_round=$(grep "C:" round*_braker_on_refprot/busco_augustus*/short_summary.specific.basidiomycota_odb10.busco_augustus*.txt |\
	sed -e 's/%//g' -e 's/\[/,/g' -e 's/]//g' -e 's/:/\t/g' -e  's/,/\t/g' |\
	LC_ALL=C sort -nr -k3 -k11n -k9n -k7n -k5  |tail -n 1 |cut -d "/" -f 1 ) 
	#LC_ALL=C sort -nr -k11 -k3 -n -k5 -k7 -k9  |tail -n 1 |cut -d "/" -f 1 ) 

echo best_round is $best_round

cd ../

run_id=${best_round%_braker_on_refprot} 
mkdir 08_best_run"$run_id"
cp 06_braker/$best_round/braker.gtf 08_best_run"$run_id"/$species.gtf

sed -i 's/_pilon//g' 08_best_run"$run_id"/$species.gtf

file=08_best_run$run_id/$species.gtf

rename_gtf.py --gtf ${file} --prefix ${species} --translation_tab translation_tab.$species --out 08_best_run"$run_id"/${species}.renamed.gtf 

### Fix TSEBRA output (source: https://github.com/Gaius-Augustus/BRAKER/issues/457 )
## Fix lack of gene_id and transcript_id tags in gtf file column 9
cat 08_best_run"$run_id"/${species}.renamed.gtf | Fix_Augustus_gtf.pl > 08_best_run"$run_id"/${species}.renamed.fixed.gtf


#-------------------------------------------------------------------------------#
#final reshaping with awk
gtf=08_best_run"$run_id"/${species}.renamed.fixed.gtf #$1 #input gtf file in general this should be "braker.gtf"
newname=$species #$2  #new name for the gtf and basename for the scaffold
newgtf=${newname}.ok.gtf

echo -e " new name is $newname \n"
echo "-----------------------"
#test if tig is present or directly replacr the start of the chromosome in chromosome 1 by the new pattern
awk '{gsub("_id \"[A-Za-z0-9-]*_","_id \""  $1 "_", $0); print }' $gtf > "$newgtf"

exit
awk -v val1=$newname  'BEGIN{FS=OFS="\t"}
        {$1=val1"_"$1; print $0}' $gtf |\
        awk '{gsub("_id \"[A-Za-z0-9-]*_","_id \""  $1 "_", $0); print }'  > "$newgtf"

