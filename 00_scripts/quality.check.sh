
#Purpose:
#making a few checks for quality (in addition to the classical busco)
#Date: 2023
#Author: QR
source config/config

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo -e "script to make a few quality checks: make blast against uniprot - run interproscan - make some other check (length of protein, etc)"
   echo " "
   echo "Usage: $0 [-s|]"
   echo "options:"
   echo "-h|--help: Print this Help."
   echo "-s|--haplo: the name of the focal haplo\t"
   echo " "
   echo "dependancies: diamond interproscan "
}


############################################################
# Process the input options.                               #
############################################################
while [ $# -gt 0 ] ; do
  case $1 in
    -s | --haplo ) haplo="$2" ; echo -e "haplotype Name is ***${haplo}*** \n" >&2;;
    -h | --help ) Help ; exit 2 ;;
   esac
   shift
done 

#cds input
input="$haplo".spliced_cds.fa

#--------------------- get uniprot ------------------#
mkdir uniprot
cd uniprot
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

#index:
diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot.fasta

cd ../

#moove to directory:
cd $haplo/08_best_run/

#---------------- run blast/diamond against the cds ----#
name=$(basename $input )

output=diamond_blastx
mkdir $output 2>/dev/null

diamond blastx -d uniprot_sprot.fasta \
        -q $input --ultra-sensitive\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq \
        --threads 20 \
        -o $output/matches."$name".tsv


#------ Inter pro scan --------------#
if [[ $interpro = "YES" ]]
then
    #get interpro-scan
    #interproscan: 
    cd ../../
    
    mkdir my_interproscan
    cd my_interproscan
    wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.66-98.0/interproscan-5.66-98.0-64-bit.tar.gz
    wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.66-98.0/interproscan-5.66-98.0-64-bit.tar.gz.md5
    
    #md5sum -c interproscan-5.66-98.0-64-bit.tar.gz.md5
    
    tar -pxvzf interproscan-5.66-98.0-*-bit.tar.gz
    python3 setup.py -f interproscan.properties
    
    path=$(pwd)
    echo -e "\n#Path to interproscan\nexport PATH=\$PATH:$path" >> ~/.bashrc
    source ~/.bashrc
    
    cd ../../$haplo/08_best_run
    interproscan.sh -i "$haplo"_prot.clean.fa -goterms 2>&1 |tee interpro.$i.log
    
fi

#make other checks:
#plot cds length/intron/exons conts, etc...
