#downlaod uniprot and index

#this will be used for later quality checks
mkdir uniprot
cd uniprot
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

#index:
diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot.fasta

