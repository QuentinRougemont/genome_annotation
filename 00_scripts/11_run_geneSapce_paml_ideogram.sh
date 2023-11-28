#!/bin/bash

#Purpose:
#master script to prepare bed files, laucnch GeneSpace, run paml and launch downstream Rscript 
#will check existence of all dependencies
#Date: 2023
#Author: QR

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo -e "master script to: \n 1 - create bed files, \n 2 - launch GeneSpace, \n 3 - run paml and launch downstream Rscripts (Rideogram, plot paml, etc)"
   echo " "
   echo "Usage: $0 [-s1|-s2|-f|-a|-h|]"
   echo "options:"
   echo " -h|--help: Print this Help."
   echo " -s1|--haplo1: the name of the first  focal haplotype\t "
   echo " -s2|--haplo2: the name of the second focal haplotype\t "
   echo " -a |--ancestral_sp: the name of the ancestral haplo to infer orthology and plot gene order"
   echo " -f|--folderpath: the path to the global folder containing haplo1 and haplo 2"
   echo " "
   echo "dependancies: orthofinder, mcscanx, GeneSpace, paml (yn00), Rideogram, translatorX minimap2"
}


###########################################################
## to do: add more support to handle the ancestral haplo: it could be either only the genome + gtf or directly a bed + protein file for instance



## --------------- TEST THAT ALL DEPENDENCIES ARE HERE --------------------------------- ##
## test if MCScanX is install

## for the section above; insert test at each command to verify that installations were successfull 
## else exit the code

command='MCScanX'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will try a manual installation through git"
    git clone https://github.com/wyp1125/MCScanX
    cd MSCanX ; make -j 8
    #if command was successfull then add to path:
    path=$(pwd)
    echo -e "\n#Path to MScanX\nexport PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../
    #exit 1
fi

## test if orthoFinder is installed
command='orthofinder'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will try a manual installation through wget"
    wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.5/OrthoFinder.tar.gz
    tar zxf OrthoFinder.tar.gz
    cd OrthoFinder
    #if command was successfull then add to path:
    path=$(pwd)
    echo -e "\n#Path to $command\nexport PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    
    #also add diamond, fastme and mcl which are present within orthofinder:
    cd bin/
    path=$(pwd)
    echo -e "\n#Path to diamond\n export PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../../
    #exit 1
fi


# test if paml/yn00 is installed
command='yn00'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will try a download and manual installation through wget"
    wget https://github.com/abacus-gene/paml/releases/download/4.10.7/paml-4.10.7-linux-X86_64.tgz
    tar zxvf paml-4.10.7-linux-X86_64.tgz
    cd paml-4.10.7/bin
    path=$(pwd)
    echo -e "\n#Path to $command\n export PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../../
fi

# test if translatorx is installed:
command='translatorx_vLocal.pl'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will try a manual installation through wget"
    wget http://161.111.160.230/cgi-bin/translatorx_vLocal.pl
    path=$(pwd)
    echo -e "\n#Path to $command\n export PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
fi

#for it we must all install muscle :
#which muscle 
command='muscle'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will try a manual installation through wget"
    #muscle5.1.linux_intel64.230/cgi-bin/translatorx_vLocal.pl
    mkdir muscle ; cd muscle
    wget wget https://github.com/rcedgar/muscle/releases/download/5.1.0/muscle5.1.linux_intel64
    ln -s muscle5.1.linux_intel64 muscle
    chmod +x muscle
    path=$(pwd)
    echo -e "\n#Path to $command\n export PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../
fi

command='R'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "I am becoming too lazy, please install R"
    exit 1
fi

command='minimap2'
if ! command -v $command &> /dev/null
then
    echo "$command could not be found"
    echo "will try a download and manual installation through wget"
    curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf -
    cd ./minimap2-2.26_x64-linux 
    path=$(pwd)
    echo -e "\n#Path to $command\n export PATH=\$PATH:$path" >> ~/.bashrc 
    source ~/.bashrc  
    cd ../
fi

echo "all dependencies succeffuly checked"


############################################################
# Process the input options.                               #
############################################################
while [ $# -gt 0 ] ; do
  case $1 in
    -s1 | --haplo1) haplo1="$2" ; echo -e "haplotype 1 Name is ***${haplo1}*** \n" >&2;;
    -s2 | --haplo2) haplo2="$2" ; echo -e "haplotype 2 Name is ***${haplo2}*** \n" >&2;;
    -a  | --ancestral_sp) ancestral_sp="$2" ; echo -e "ancestral haplo  Name is ***${ancestral_sp}*** \n" >&2;;
    -f  | --folderpath  ) folderpath="$2"  ; echo -e "global folder is  ${folderpath} \n" >&2;;
    -h  | --help ) Help ; exit 2 ;;
   esac
   shift
done 

if [ -z "${haplo1}" ] || [ -z "${haplo2}" ] || [ -z "${folderpath}" ]  || [ -z "${ancestral_sp}" ]    ; then
	Help
	exit 2
fi

#------------------------------ step 1 prepare bed file for each haplo -------------------------------------#
#really simple:

cd $folderpath

#remove any existing folder:
rm genespace peptide paml plots -rf 2>/dev/null
#to do: insert alert message for the user.
mkdir -p genespace/bed genespace/peptide paml plots

# create bed
awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$10}' $haplo1/08_best_run/$haplo1.longest_transcript.gtf |sed 's/"//g' > genespace/bed/$haplo1.bed
awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$10}' $haplo2/08_best_run/$haplo2.longest_transcript.gtf |sed 's/"//g' > genespace/bed/$haplo2.bed

# simplify the protein file to match the bed (i.e. remove the _1 inserted by transeq and the CDS length info):
sed 's/_1 CDS=.*$//g' $haplo1/08_best_run/"$haplo1"_prot.fa > genespace/peptide/$haplo1.fa
sed 's/_1 CDS=.*$//g' $haplo2/08_best_run/"$haplo2"_prot.fa > genespace/peptide/$haplo2.fa

#verify that IDs in bed and fasta file are matching - else exit  
grep ">" genespace/peptide/$haplo1.fa |sed 's/>//g' > tmp1
grep ">" genespace/peptide/$haplo2.fa |sed 's/>//g' > tmp2

check1=$(grep -Ff tmp1 genespace/bed/$haplo1.bed |wc -l )
check2=$(grep -Ff tmp2 genespace/bed/$haplo2.bed |wc -l )

echo -e "check2 size is $check2"
echo -e "check1 size is $check1"

bedsize1=$(wc -l genespace/bed/$haplo1.bed |awk '{print $1}' )
bedsize2=$(wc -l genespace/bed/$haplo2.bed |awk '{print $1}' )

echo -e "bedisze2  size is $bedsize2"
echo -e "bedisze1  size is $bedsize1"

#check that all is matching:
if [ "$bedsize1" = "$check1" ]
then
	echo "input1 is ok" 
	rm tmp1
else
	echo "input1 is not ok"
	echo "check your data"
	exit 2
fi

if [ "$bedsize2" = "$check2" ]
then
	echo "input2 is ok" 
	rm tmp2
else
	echo "input2 is not ok"
	echo "check your data"
	exit 2
fi

# -- handling ancestral haplo ------
# -- this part assumes that a bed and peptide file are existant for the ancestral haplo
# -- here we used a genome annotated with the same pipeline relying on braker 

cd genespace/bed/
ln -s ../../../"$ancestral_sp"/"$ancestral_sp".bed . 
cd ../peptide
ln -s ../../../"$ancestral_sp"/"$ancestral_sp".prot.fa "$ancestral_sp".fa

cd ../../

#------------------------------ step 2 run GeneSpace ---------------------------------------------------------#
cd genespace 

Rscript ../../00_scripts/Rscripts/01.run_geneSpace.R

if [ $? -eq 0 ]; then
    echo genespace worked successfully
else
    echo genespace failed
    exit 1
fi

#plot genespace subspace of target chromosomes: 
Rscript ../../00_scripts/Rscripts/02.plot_geneSpace.R

cd ../
#------------------------------ step 3 run paml  -------------------------------------------------------------#

scaffold=scaffold.txt #hardcoded variable to be passed as an argument for later
cp ../$scaffold .
echo $(pwd)
echo haplo1 is "$haplo1"
echo haplo2 is "$haplo2"

../00_scripts/12_command_line.paml.sh -h1 "$haplo1" -h2 "$haplo2" -s "$scaffold" -a "$ancestral_sp"

#here insert a test to veryfy that previous code was successful and else exit
if [ $? -eq 0 ]; then
    echo paml worked successfully
else
    echo paml failed
    exit 1
fi


pamlsize=$(wc -l paml/results_YN.txt |awk '{print $1}' ) 
scpo=$(wc -l paml/single.copy.orthologs |awk '{print $1}' )

echo -e "there is $pamlsize results for PAML \n"
echo -e "there is $scpo single copy orthologs \n" 


#----------------------------------- step4 -- plot paml results  -----------------------------------------#
#test if previous step was successfull else plot or exit with high levels of pain
echo $ancestral_sp 

Rscript ../00_scripts/Rscripts/03.plot_paml.R $ancestral_sp $haplo1 $haplo2


# -- step5 -- plot ideogram 
#test if previous step was successfull else plot or exit with high levels of pain
samtools faidx $haplo1/03_genome/"$haplo1".fa
samtools faidx $haplo2/03_genome/"$haplo2".fa

Rscript ../00_scripts/Rscripts/04.ideogram.R $haplo1 $haplo2 #add links!

#
## --------------------------------Make Synteny table -----------------------------------------------
is_anc='TRUE'
if [ -n ${ancestral_sp} ] ; then
	is_anc='TRUE'
else

	is_anc='FALSE'
fi
	
path_orthofinder='genespace/orthofinder/Results_*/'
path_bed='genespace/bed/'

python3 ../00_scripts/utility_scripts/02.Make_synteny_table.py ${haplo1} ${haplo2} ${path_orthofinder} ${path_bed} ${is_anc} ${ancestral_sp}



# ---------------------------------- step6 -- create circos plot ----------------------------------------#
#circos plot here:
Rscript ../00_scripts/Rscripts/05_plot_circos.R $haplo1 $ancestral_sp $scaffolds $genes_plot
Rscript ../00_scripts/Rscripts/05_plot_circos.R $haplo2 $ancestral_sp $scaffolds $genes_plot


#-- step7 -- run minimap between the genomes 
#run minimap on the genome 
#assumption : each genome MUST BE located in folder 03-genome
minimap2 -cx asm5 $haplo1/03_genome/"$haplo1".fa $haplo2/03_genome/"$haplo2".fa > aln."$haplo1"_"$haplo2".paf 

if [ -n ${ancestral_sp} ] ; then
    minimap2 -cx asm5 $ancestral_sp/$ancestral_sp.fa $haplo2/03_genome/"$haplo2".fa > aln."$ancestral_sp"_"$haplo2".paf 
    minimap2 -cx asm5 $ancestral_sp/$ancestral_sp.fa $haplo1/03_genome/"$haplo1".fa > aln."$ancestral_sp"_"$haplo1".paf 
    #preparing scaffold to highlight in dotplot:
    awk '{gsub("_","\t",$0) ; print $2"_"$3"_"$4"\t"$6"_"$7}' paml/single.copy.orthologs|sort |uniq -c|awk '$1>10 ' > scaff.anc.haplo1.txt
    awk '{gsub("_","\t",$0) ; print $2"_"$3"_"$4"\t"$9"_"$10}' paml/single.copy.orthologs|sort |uniq -c|awk '$1>10 ' > scaff.anc.haplo2.txt
    awk '{gsub("_","\t",$0) ; print $6"_"$7"\t"$9"_"$10}' paml/single.copy.orthologs|sort |uniq -c|awk '$1>10 ' > scaff.haplo1.haplo2.txt 

    Rscript ../00_scripts/Rscripts/dotplot_paf.R  aln."$haplo1"_"$haplo2".paf 
    Rscript ../00_scripts/Rscripts/dotplot_paf.R  aln."$ancestral_sp"_"$haplo1".paf 
    Rscript ../00_scripts/Rscripts/dotplot_paf.R  aln."$ancestral_sp"_"$haplo2".paf 

    Rscript ../00_scripts/Rscripts/synteny_plot.R aln."$ancestral_sp"_"$haplo1".paf scaff.anc.haplo1.txt 
    Rscript ../00_scripts/Rscripts/synteny_plot.R aln."$ancestral_sp"_"$haplo2".paf scaff.anc.haplo2.txt 
    Rscript ../00_scripts/Rscripts/synteny_plot.R aln."$haplo1"_"$haplo2".paf scaff.haplo1.haplo2.txt 

else 
    awk '{gsub("_","\t",$0) ; print $2"_"$3"\t"$5"_"$6}' paml/single.copy.orthologs|sort |uniq -c|awk '$1>10 ' > scaff.haplo1.haplo2.txt
    
    #then run pafr to generate a whole genome dotplot and eventually dotplot for some target scaffold:
    Rscript ../00_scripts/Rscripts/dotplot_paf.R  aln."$haplo1"_"$haplo2".paf 
    Rscript ../00_scripts/Rscripts/synteny_plot.R aln."$haplo1"_"$haplo2".paf scaff.haplo1.haplo2.txt 

fi


#we can also run Rideogram here
