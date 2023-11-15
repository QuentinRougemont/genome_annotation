#!/bin/bash

#script to rename the scaffold in the gtf, rename the gtf, create a synchronised genome
#extract longest protein

#To do: insert help here 


#required variable here:
if [ $# -ne 2  ]; then
    echo "USAGE: $0 species_name RNAseq "
    echo "Expecting the species name and a YES/NO variable stating stating wether RNAseq was provided"
    exit 1
else
    species=$1
    RNAseq=$2
    echo "Species Name is : ${species}"
    echo "*******************"
    echo "annotation was performed with RNAseq ? ${RNAseq}"
    echo "*******************"
fi

#test if rnaseq ici 

#----------------------------------------------------------------------#
#1 -- find best run
echo -e "\nfinding best run \n" 
cd 06_braker
best_round=$(grep "C:" round*_braker_on_refprot/busco_augustus*/short_summary.specific.basidiomycota_odb10.busco_augustus*.txt |\
sed -e 's/%//g' -e 's/\[/,/g' -e 's/]//g' -e 's/:/\t/g' -e  's/,/\t/g' |\
LC_ALL=C sort -nr -k3 -k11n -k9n -k7n -k5  |tail -n 1 |cut -d "/" -f 1 ) 
#LC_ALL=C sort -nr -k11 -k3 -n -k5 -k7 -k9  |tail -n 1 |cut -d "/" -f 1 ) 

echo -e "best_round is $best_round\n----------------------------------------------"

cd ../

if [[ $RNAseq = "YES" ]]
then
	echo -e "\nrunning tsebra\n" 

	#2 -- run tsebra
	./00_scripts/08_tsebra.sh $species $best_round

	mkdir 08_best_run
	cd 08_best_run
	ln -s ../07-tsebra_results/$species.combined.gtf $species.gtf
	cd ../
	file=08_best_run/$species.gtf
	nb_genes=$(awk '$3=="gene" ' $file |wc -l)
	echo -e "\nthere is $nb_genes genes after running tsebra\n" 

else 
	#----------------------------------------------------------------------#
	#2 -- copy best run based on database in a final folder
	#run_id=${best_round%_braker_on_refprot} 
	#mkdir 08_best_run"$run_id"
	#cp 06_braker/$best_round/braker.gtf 08_best_run"$run_id"/$species.gtf
	#sed -i 's/_pilon//g' 08_best_run"$run_id"/$species.gtf
	#file=08_best_run$run_id/$species.gtf
	
	echo "running on protein only - not running tsebra" 

	mkdir 08_best_run
	cp 06_braker/$best_round/braker.gtf 08_best_run/$species.gtf
	file=08_best_run/$species.gtf

	nb_genes=$(awk '$3=="gene" ' $file |wc -l)
        echo -e "\nthere is $nb_genes in the best protein run\n"
fi

#------------------------------ step 3 ---------------------------------------------------------#
echo -e  "\n------------renaming and fixing braker output now------------\n" 
#then this is a part common to both RNAseq + Proteins or Proteins only:
#rename_gtf.py --gtf ${file} --prefix ${species} --translation_tab translation_tab.$species --out 08_best_run"$run_id"/${species}.renamed.gtf 
rename_gtf.py --gtf ${file} --prefix ${species} --translation_tab translation_tab.$species --out 08_best_run/${species}.renamed.gtf 

# Fix TSEBRA output (source: https://github.com/Gaius-Augustus/BRAKER/issues/457 )
# Fix lack of gene_id and transcript_id tags in gtf file column 9
#cat 08_best_run"$run_id"/${species}.renamed.gtf | Fix_Augustus_gtf.pl > 08_best_run"$run_id"/${species}.renamed.fixed.gtf
cat 08_best_run/${species}.renamed.gtf | Fix_Augustus_gtf.pl > 08_best_run/${species}.renamed.fixed.gtf


#------------------------------ step 4 --------------------------------------------------------#
#final reshaping with awk
#gtf=08_best_run"$run_id"/${species}.renamed.fixed.gtf #$1 #input gtf file in general this should be "braker.gtf"
gtf=08_best_run/${species}.renamed.fixed.gtf #$1 #input gtf file in general this should be "braker.gtf"

newgtf=${species}.ok.gtf

#test if tig is present or directly replacr the start of the chromosome in chromosome 1 by the new pattern
awk '{gsub("_id \"[A-Za-z0-9-]*_","_id \""  $1 "_", $0); print }' $gtf > 08_best_run/"$newgtf"

#exit
#awk -v val1=$newname  'BEGIN{FS=OFS="\t"}
#        {$1=val1"_"$1; print $0}' $gtf |\
#        awk '{gsub("_id \"[A-Za-z0-9-]*_","_id \""  $1 "_", $0); print }'  > "$newgtf"
#to run me on all species do:
#for file in $(find . -name "*ok.gtf" ) ; do base=$(basename $file ) ; path=${file%/*} ; echo -e "cd $path\n../../00_script/02_remove_redundant_genemark.sh $base ; \ncd ../../" >> clean.gtf.tmp  ; done
#bash clean.gtf.tmp 

#------------------------------ step 5 extracting prot and cds from new files  --------------------------------------------------------#
echo -e "\n#------------------------------ step 5 extracting prot and cds from new files  --------------------------------------------------------#"

mkdir 08_best_run/01_species_cds 2>/dev/null
mkdir 08_best_run/02_species_prot 2>/dev/null

gtf=$newgtf
gtffull=08_best_run/$gtf

#extract the cds and protein from it:
output="$species"_cds.fa
echo output cds is $output
gffread -w 08_best_run/01_species_cds/$output -g 03_genome/"${species}".fa $gtffull

#here capture error if genome and gtf do not match and exit with error 1 
#here if any error occurs capture it!

#then convert also the file to its cds:
echo "translate CDS into amino acid "
transeq -sequence 08_best_run/01_species_cds/$output -outseq 08_best_run/02_species_prot/$species.prot

samtools faidx 08_best_run/02_species_prot/$species.prot

echo "-----------------------------------------------------------------"
echo "extract longest transcript" 
echo "-----------------------------------------------------------------"

cd 08_best_run/02_species_prot
#assumption : all transcript finishes by ".t1, .t2, .t3 so the dot (.) iis the delimiter
awk -F "." '{print $1"\t"$0}' $species.prot.fai |\
        awk '$3>max[$1]{max[$1]=$3; row[$1]=$2} END{for (i in row) print row[i]}' > longest.transcript.tmp

#linearize file so that the next command will work:
awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' "$species".prot > "$species".prot.lin.fasta

grep -A1 -Ff longest.transcript.tmp "$species".prot.lin.fasta > "$species".longest_transcript.fa
rm longest.transcript.tmp
rm $species.prot.fai
cd ../../


echo -e "\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ step 6 : cleaning the GTF based on non-overlapping transcript ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n#"
echo "-----------------------------------------------------------------"
echo "remove redundant Gene in the CDS and Protein fasta files"
echo "-----------------------------------------------------------------"

#the purpose of the line of code is to remove redundant GeneMark annotation, that is GeneMark gene with same start and end position as those from Augustus:
cd 08_best_run

#sort -k1,1 -k4,4n  $gtf |cut -f 1-5,9 |awk '$3~"CDS" {print $1"_"$4"_"$5"\t"$1"_"$4"\t"$1"_"$5"\t"$0}' > tmp1
#awk 'NR == FNR {count[$1]++; next} count[$1]>1' tmp1 tmp1 |grep "GeneMark" |cut -f 9 |cut -d ";" -f 1 |uniq |sed -e 's/transcript_id //' -e 's/"//g' -e 's/ //g' > toremove
#wc -l toremove
#exact_overlap=$(wc -l toremove |awk '{print $1}' )
#echo -e "\nthere is $exact_overlap GeneMark CDS with the exact same chromosome and end position"
#
#
#rem=$(wc -l toremove |awk '{print $1}')
#echo -e "find a total of $rem redondant transcript\n------------------------------------"
#
#mkdir 03_cleaned_prot 2>/dev/null
#sed 's/ CDS=.*//g' 02_species_prot/"$species".longest_transcript.fa|\
#	sed '/^--$/d' |\
#	awk -vRS=">" -vOFS="\t" ' {$1=">"$1; print $0}' |\
#	grep -vFf toremove |\
#	sed 's/\t/\n/g' |\
#	sed 1d  > 03_cleaned_prot/"$species"_prot.fa
#
##rm tmp1 toremove

cut -f 9 "$gtf" |sed 's/^gene_id//g' |sed 's/transcript_id//g' |cut -d ";" -f 1 |sed 's/.t[1-9]//g' |sort |uniq -c |awk '$1==1 {print}'  > fragmented_gene.txt 

loss=$(wc -l fragmented_gene.txt |awk '{print $1}' )
echo " there is $loss fragmented gene"
echo "finished"

#species=$1 #name of the species 

#-------------------------------
#protpath=03_cleaned_prot
protpath=02_species_prot
#----------------------------------

#declare full path to input:
#prot="$protpath"/"$species"_prot.fa
prot="$protpath"/"$species".longest_transcript.fa

echo -e "there is $(wc -l $gtf |awk '{print $1}') lines in $gtf" 
#----------------------------------------------------
# subset our gtf to keep only the cds in the cds files!
#grep ">" $prot |sed 's/>//g' |sed 's/_1$//g'  > wanted.cds.tmp 
grep ">" $prot |sed 's/>//g' |sed 's/_1 CDS=.*//g'  > wanted.cds.tmp 
sed 's/.t[0-9]$//g' wanted.cds.tmp > wanted.gene.tmp

#ok so here I tried many awk solution to have both the wanted transcript and the gene using NF == FNR, etc, nothing work so I used a
#two step process creating a p1 file for the whole stuf and p2 for the gene, which is not really what I wanted.
pwd

grep -Ff wanted.cds.tmp "$gtf" > p1 #"$name".sub.gtf
grep -f wanted.gene.tmp  <(awk '$3=="gene" ' $gtf)  > p3

#ideally I want to sort on the gene, then transcript, then CDS, then exon as well
cat p1 p3 |LC_ALL=C sort -k1,1 -k4,4n -k5,5n > "$species".longest_transcript.gtf

echo -e there is $(wc -l "$species".longest_transcript.gtf |awk '{print $1}' ) lines in "$species".longest_transcript.gtf 

#declare new gtf for new work:
gtf="$species".longest_transcript.gtf

#--------------------- now removing fully overlapping CDS #
sort -k1,1 -k4,4n  $gtf |cut -f 1-5,9 |awk '$3~"CDS" {print $1"_"$4"_"$5"\t"$1"_"$4"\t"$1"_"$5"\t"$0}' > tmp1
awk 'NR == FNR {count[$2]++; next} count[$2]>1' tmp1 tmp1 |cut -f 9 |cut -d ";" -f 1 |uniq |sed -e 's/transcript_id //' -e 's/"//g' -e 's/ //g' > toremove
exact_overlap=$(wc -l toremove |awk '{print $1}' )
echo -e "\nthere is $exact_overlap CDS with the exact same chromosome and start  "
awk 'NR == FNR {count[$3]++; next} count[$3]>1' tmp1 tmp1 |cut -f 9 |cut -d ";" -f 1 |uniq |sed -e 's/transcript_id //' -e 's/"//g' -e 's/ //g' >> toremove
exact_overlap=$(wc -l toremove |awk '{print $1}' )
echo -e "\nthere is $exact_overlap CDS with the exact same chromosome and end  "


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ step 7 : re-extracting the proteins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo -e "\n#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ step 7 : re-extracting the proteins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n"

gffread -w "$species".spliced_cds.fa -g ../03_genome/$species.fa "$species".longest_transcript.gtf
echo "translate CDS into amino acid "
source /local/env/envemboss-6.6.0.sh
transeq -sequence "$species".spliced_cds.fa -outseq "$species"_prot.fa

#now run busco to validate. the score should be close to the initial score 
