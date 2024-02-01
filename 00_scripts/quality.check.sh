
#Purpose:
#making a few checks for quality (in addition to the classical busco)
#Date: 2024
#Author: QR
source ../../config/config

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


#---------------- run blast/diamond against the cds ----#
name=$(basename $input )

output=diamond_blastx
mkdir $output 2>/dev/null

diamond blastx -d ../../uniprot/uniprot_sprot.fasta \
        -q $input --ultra-sensitive\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq \
        --threads 20 \
        -o $output/matches."$name".tsv


#------ Inter pro scan --------------#
if [[ $interpro = "YES" ]]
then
    interproscan.sh -i "$haplo"_prot.clean.fa -goterms 2>&1 |tee interpro.log
    
fi

#make other checks:
#plot cds length/intron/exons conts, etc...
