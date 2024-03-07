#downlaod uniprot and index

#this will be used for later quality checks
if [ -f uniprot/uniprot_sprot.fasta ] ;
then
    echo uniprot already present 

    #check diamond
    if [ -f uniprot/uniprot_sprot.fasta.dmd ] ;
    then
        echo uniprot already indexed 
    else#index:
        diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot.fasta
    fi

else
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
fi

