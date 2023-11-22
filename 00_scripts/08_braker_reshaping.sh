#!/bin/bash

#Purpose:
#script to rename the scaffold in the gtf, rename the gtf, create a synchronised genome
#extract longest protein
#Date: 2023
#Author: QR

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo -e "script to post-process braker gtf:\n 1 - run tsebra \n 2 - rename genes id \n 3 - find longest transcript \n 4 - create gtf containing non overlapping transcript \n 5 - extract protein and CDS fasta"
   echo " "
   echo "Usage: $0 [-s|r|]"
   echo "options:"
   echo "-h|--help: Print this Help."
   echo "-s|--haplo: the name of the focal haplo\t will be used to rename the genes"
   echo "-r|--RNAseq: a float YES/NO stating whether RNAseq was used to annotate the genome "
   echo " "
   echo "dependancies: TSEBRA, samtools, gffread, transeq busco "
}


## --------------- TEST THAT ALL DEPENDENCIES ARE HERE --------------------------------- ##
## for the section above; insert test at each command to verify that installations were successfull 
## else exit the code

command='gffread'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will try a manual installation through git"
    git clone https://github.com/gpertea/gffread
    cd gffread
    make release
    #if command was successfull then add to path:
    if [ $? -eq 0 ]; then
        echo gffread installation worked successfully
        path=$(pwd)
        echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
        source ~/.bashrc  
        cd ../

    else
       echo installation failed\nmake sur to have git, make and gclib
       exit 1
    fi
fi

command='transeq'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will try a manual installation through conda/mamba"
    #use mamba here
    mambe create -p braker_env emboss=6.6.0
    #if command was successfull then add to path:
    if [ $? -eq 0 ]; then
        echo emboss and transeq installation successfull

    else
       echo installation failed\nmake sur to have git, make and gclib
       exit 1
    fi
fi

#samtools: 
command='samtools'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will try a manual installation through conda/mamba"
    #use mamba here
    eval "$(conda shell.bash hook)"
    mamba activate braker_env 
    mamba install -c bioconda samtools samtools=1.18
    #if command was successfull then add to path:
    if [ $? -eq 0 ]; then
        echo samtools installation successfull

    else
       echo installation failed\nmake sur to have git, make and gclib
       exit 1
    fi
fi

command='busco'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will try a manual installation through conda/mamba"
    #use mamba here
    eval "$(conda shell.bash hook)"
    mamba activate braker_env 
    mamba install -c bioconda busco=5.5.0

    #if command was successfull then add to path:
    if [ $? -eq 0 ]; then
        echo busco installation successfull

    else
       echo installation failed\nmake sur to have git, make and gclib
       exit 1
    fi
fi


#tsebra:
#git clone https://github.com/Gaius-Augustus/TSEBRA + export path to bin
command='tsebra'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will try a manual installation through conda/mamba"
    git clone https://github.com/Gaius-Augustus/TSEBRA 
    #+ export path to bin
    cd TSEBRA/bin
    if [ $? -eq 0 ]; then
        echo tsebra installation worked successfully
        path=$(pwd)
        echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
        source ~/.bashrc  
        cd ../

    else
       echo installation failed\nmake sur to have git,
       exit 1
    fi
fi

############################################################
# Process the input options.                               #
############################################################
while [ $# -gt 0 ] ; do
  case $1 in
    -s | --haplo ) haplo="$2" ; echo -e "haplotype Name is ***${haplo}*** \n" >&2;;
    -r | --rnaseq  ) RNAseq="$2"  ; echo -e "annotation was performed with RNAseq ? ${RNAseq} \n" >&2;;
    -h | --help ) Help ; exit 2 ;;
   esac
   shift
done 


if [ -z "${haplo}" ] || [ -z "${RNAseq}" ]  ; then
	Help
	exit 2 
fi

#------------------------------ step 1 find best run ---------------------------------------------------------#
echo -e  "\n-----------------------------------------------------------------"
echo -e "\nfinding best run \n" 
echo -e  "-----------------------------------------------------------------\n"

cd $haplo/06_braker

#there is a bit of variance between braker run train on a protein database,
#therefore we perform a few replicate run and choose the 'best run' where best is defined as the run providing 
#the highest score in terms of busco completness 
#we sort on k3 (highest completeness), k11 (lowest amount of missing), k9 (lowest amount of fragmented genes)  
best_round=$(grep "C:" round*_braker_on_refprot/busco_augustus*/short_summary.specific.basidiomycota_odb10.busco_augustus*.txt |\
sed -e 's/%//g' -e 's/\[/,/g' -e 's/]//g' -e 's/:/\t/g' -e  's/,/\t/g' |\
LC_ALL=C sort -k3nr -k11n -k9n -k7n -k5  |head -n 1 |cut -d "/" -f 1 ) 
#alternative sorting: 
#LC_ALL=C sort  -k11 -k3 -n -k5 -k7 -k9  |tail -n 1 |cut -d "/" -f 1 ) #minimising missing first 


echo -e "best_round is $best_round\n----------------------------------------------"

cd ../

#------------------------------ step 2 finding runs --------------------------------------------#

if [[ $RNAseq = "YES" ]]
then
	echo -e "\nrunning tsebra\n" 

	#2 -- run tsebra
	#test if default.cfg can be found:
	tsebraconf=config/default.cfg
	#[ -f $tsebraconf ] && { echo " $tsebraconf exist";  cp "$tsebraconf" . } || {echo "no config file for tsebra"; exit } 
	if [ -f $tsebraconf ] ; then 
		echo " $tsebraconf exist";  cp "$tsebraconf" . 
	else 	
		echo "FATAL ERROR: no config file for tsebra"
		exit
	fi

	#then run tsebra:
	./00_scripts/09_tsebra.sh $haplo $best_round

	mkdir 08_best_run
	cd 08_best_run
	ln -s ../07-tsebra_results/$haplo.combined.gtf $haplo.gtf
	cd ../
	file=08_best_run/$haplo.gtf
	nb_genes=$(awk '$3=="gene" ' $file |wc -l)
	echo -e "\nthere is $nb_genes genes after running tsebra\n" 

else 
	#----------------------------------------------------------------------#
	#2 -- copy best run based on database in a final folder

	#DEPRCETATED :
	#run_id=${best_round%_braker_on_refprot} 
	#mkdir 08_best_run"$run_id"
	#cp 06_braker/$best_round/braker.gtf 08_best_run"$run_id"/$haplo.gtf
	#sed -i 's/_pilon//g' 08_best_run"$run_id"/$haplo.gtf
	#file=08_best_run$run_id/$haplo.gtf
	
	echo "running on protein only - not running tsebra" 

	mkdir 08_best_run
	cp 06_braker/$best_round/braker.gtf 08_best_run/$haplo.gtf
	file=08_best_run/$haplo.gtf

	nb_genes=$(awk '$3=="gene" ' $file |wc -l)
        echo -e "\nthere is $nb_genes in the best protein run\n"
fi

#------------------------------ step 3 ---------------------------------------------------------#
echo -e  "\n-----------------------------------------------------------------"
echo -e  "\n------------renaming and fixing braker output now------------\n" 
echo -e  "-----------------------------------------------------------------\n"

#then this is a part common to both RNAseq + Proteins or Proteins only:
rename_gtf.py --gtf ${file} --prefix ${haplo} --translation_tab translation_tab.$haplo --out 08_best_run/${haplo}.renamed.gtf 

# Fix TSEBRA output (source: https://github.com/Gaius-Augustus/BRAKER/issues/457 )
# Fix lack of gene_id and transcript_id tags in gtf file column 9
cat 08_best_run/${haplo}.renamed.gtf | ../00_scripts/Fix_Augustus_gtf.pl > 08_best_run/${haplo}.renamed.fixed.gtf


#------------------------------ step 4 --------------------------------------------------------#
#rename the gene/cds/transcript/exons/ ID in the gtf so they contain the chromosome and haplo name (easier to parse)
#1 - declare the gtf:
gtf=08_best_run/${haplo}.renamed.fixed.gtf #current gtf
newgtf=${haplo}.ok.gtf #future gtf

#2 - renaming with a simple command: 
awk '{gsub("_id \"[A-Za-z0-9-]*_","_id \""  $1 "_", $0); print }' $gtf > 08_best_run/"$newgtf"

#awk -v val1=$newname  'BEGIN{FS=OFS="\t"}
#        {$1=val1"_"$1; print $0}' $gtf |\
#        awk '{gsub("_id \"[A-Za-z0-9-]*_","_id \""  $1 "_", $0); print }'  > "$newgtf"

#------------------- step 5 extracting prot and cds from new files  ---------------------------#
echo -e  "\n-----------------------------------------------------------------"
echo "extract protein and cds from the renamed gtf" 
echo -e  "-----------------------------------------------------------------\n"

mkdir 08_best_run/01_haplo_cds 2>/dev/null
mkdir 08_best_run/02_haplo_prot 2>/dev/null

gtf=$newgtf
gtffull=08_best_run/$gtf

#extract the cds and protein from it:
output="$haplo"_cds.fa
echo output cds is $output
gffread -w 08_best_run/01_haplo_cds/$output -g 03_genome/"${haplo}".fa $gtffull

#here capture error if genome and gtf do not match and exit with error 1 
#here if any error occurs capture it!

#then convert also the file to its cds:
echo "translate CDS into amino acid "
transeq -sequence 08_best_run/01_haplo_cds/$output -outseq 08_best_run/02_haplo_prot/$haplo.prot

samtools faidx 08_best_run/02_haplo_prot/$haplo.prot
#note: will replace this by awk code to remove samtools dependency here

echo -e "\n-----------------------------------------------------------------"
echo "extract longest transcript" 
echo -e  "-----------------------------------------------------------------\n"

cd 08_best_run/02_haplo_prot
#assumption : all transcript finishes by ".t1, .t2, .t3 so the dot (.) iis the delimiter

#note to do: change code to relax the assumption of a single dot delimiter
awk -F "." '{print $1"\t"$0}' $haplo.prot.fai |\
        awk '$3>max[$1]{max[$1]=$3; row[$1]=$2} END{for (i in row) print row[i]}' > longest.transcript.tmp

#linearize file so that the next command will work:
awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' "$haplo".prot > "$haplo".prot.lin.fasta

grep -A1 -Ff longest.transcript.tmp "$haplo".prot.lin.fasta > "$haplo".longest_transcript.fa
rm longest.transcript.tmp
rm $haplo.prot.fai
cd ../../


#~~~~~~~~~~~~~~~~~~~~~~ step 6 : cleaning the GTF based on non-overlapping transcript ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo -e "\n-----------------------------------------------------------------"
echo "remove redundant Gene in the CDS and Protein fasta files"
echo "-----------------------------------------------------------------\n"

cd 08_best_run

#look for putative fragented gene (this should be zero):
cut -f 9 "$gtf" |\
	sed 's/^gene_id//g' |\
	sed 's/transcript_id//g' |\
	cut -d ";" -f 1 |\
	sed 's/.t[1-9]//g' |\
	sort |uniq -c |\
	awk '$1==1 {print}'  > fragmented_gene.txt 

loss=$(wc -l fragmented_gene.txt |awk '{print $1}' )
echo " there is $loss fragmented gene"


#  declare full path to input:
protpath=02_haplo_prot
prot="$protpath"/"$haplo".longest_transcript.fa
echo -e "there is $(wc -l $gtf |awk '{print $1}') lines in $gtf" 

#----------------------------------------------------
# subset our gtf to keep only the cds in the cds files!
#grep ">" $prot |sed 's/>//g' |sed 's/_1$//g'  > wanted.cds.tmp 
grep ">" $prot |sed 's/>//g' |sed 's/_1 CDS=.*//g'  > wanted.cds.tmp 
sed 's/.t[0-9]$//g' wanted.cds.tmp > wanted.gene.tmp

#we keep the cds:
grep -Ff wanted.cds.tmp "$gtf" > p1 
#now the genes:
grep -f wanted.gene.tmp  <(awk '$3=="gene" ' $gtf)  > p3

#ideally I want to sort on the gene, then transcript, then CDS, then exon as well
#concatenate gene and cds/etc:
cat p1 p3 |LC_ALL=C sort -k1,1 -k4,4n -k5,5n > "$haplo".longest_transcript.gtf

echo -e there is $(wc -l "$haplo".longest_transcript.gtf |awk '{print $1}' ) lines in "$haplo".longest_transcript.gtf 

#declare new gtf for new work:
gtf="$haplo".longest_transcript.gtf

# now removing fully overlapping CDS 
#create an ID of the type chromosome_start and chromosome_end to filter :
sort -k1,1 -k4,4n  $gtf |cut -f 1-5,9 |awk '$3~"CDS" {print $1"_"$4"_"$5"\t"$1"_"$4"\t"$1"_"$5"\t"$0}' > tmp1

#we remove any CDS that has identical chromosome and start 
awk 'NR == FNR {count[$2]++; next} count[$2]>1' tmp1 tmp1 |cut -f 9 |cut -d ";" -f 1 |uniq |sed -e 's/transcript_id //' -e 's/"//g' -e 's/ //g' > toremove
exact_overlap=$(wc -l toremove |awk '{print $1}' )
echo -e "\nthere is $exact_overlap CDS with the exact same chromosome and start  "

#we remove any CDS that has identical chromosome and end position:
awk 'NR == FNR {count[$3]++; next} count[$3]>1' tmp1 tmp1 |cut -f 9 |cut -d ";" -f 1 |uniq |sed -e 's/transcript_id //' -e 's/"//g' -e 's/ //g' >> toremove
exact_overlap=$(wc -l toremove |awk '{print $1}' )
echo -e "\nthere is $exact_overlap CDS with the exact same chromosome and end  "


#~~~~~~~~~~~~~~~~~~~~~~ step 7 : re-extracting the proteins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
echo -e "\n-----------------------------------------------------------------"
echo "re-extracting protein and CDS from the final non-redundan gtf" 
echo "-----------------------------------------------------------------\n"


gffread -w "$haplo".spliced_cds.fa -g ../03_genome/$haplo.fa "$haplo".longest_transcript.gtf
echo "translate CDS into amino acid "
transeq -sequence "$haplo".spliced_cds.fa -outseq "$haplo"_prot.fa

#now run busco to validate. the score should be close to the initial score 
#insert call to busco here
