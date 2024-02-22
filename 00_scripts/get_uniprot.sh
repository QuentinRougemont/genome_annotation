#downlaod uniprot and index

#this will be used for later quality checks
rm -r uniprot 2>/dev/null #remove any existing folder to download latest release

mkdir uniprot
cd uniprot
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

#index:
diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot.fasta
