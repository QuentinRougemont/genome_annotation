#downlaod uniprot and index

#this will be used for later quality checks
rm -r uniprot 2>/dev/null #remove any existing folder to download latest release

mkdir uniprot
cd uniprot

echo -e "\n\n------downloading uniprot--------\n"
wget -q https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
if [ $? -eq 0 ]; then
    echo -e "\n---- download sucessfull\n make diamond db for later checks ------\n"
    gunzip uniprot_sprot.fasta.gz

    #index:
    diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot.fasta
else
   echo -e "download of uniprot failed\n
       check your internet connexion"   
   exit 
fi


