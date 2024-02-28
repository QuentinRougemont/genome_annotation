
* If **annotation** is required then the following data are recquired:

	* a list of protein from the same or closely related species in fasta format

	* if possible: a custom database of TE for the TE annotation steps

	* a lineage name for busco evaluation of the genome annotation quality.  

The list of lineage can be obtained [here](https://busco-data.ezlab.org/v5/data/lineages/) or using :

```shell
busco --list-dataset
``` 

in the config/config file set the **busco_lineage** by providing the name of the species that is closer to your study organism 

* A list of scaffold from the targeted region of interest (sex chromosome supergenes). Should be as follows:


| Genome name        | chromosome      |     Order  | 
|:-------------------|:-----------------|-----------| 
| __haplo1__         | chrX_scaff01    	| N         | 
| __haplo1__         | chrX_scaff02     | R         | 


Header is not needed. 

N = Normal, R = Reversed, meaning that the scaffold orientation should be reversed.



* optional:  1 ancestral genome in **fasta** format with its annotation in **gtf** format.  

**!keep a single transcript per gene!**  

the file is as follows:

```shell
cat config/config
# config file
#--- COMPULSORY MINIMAL LEVEL OF INFORMATION REQUIRED -----
genome1=""     #full path to current genome1 assembly (fasta for haplotype1 - compressed or not)
genome2=""     #full path to current genome2 assembly (fasta for haplotype2 - compressed or not)

#----- optional --------
current_name1=" "   #current basename of scaffolds in the genome (ignore if names are already short like "chr1, chr2, etc")
haplotype1=""       #name1 [name of haplotype1 - will be used to rename the genome and the contigs inside the genome]
current_name2=" "   #current basename of scaffolds in the genome2 (ignore if name are already like "chr1, chr2, etc")
haplotype2=""       #name2 [name of haplotype2 - will be used to rename the genome and the contigs inside the genome]

ancestral_genome="" #full path to the ancestral genome used as proxy of ancestral gene order (fasta, compressed or not)
ancestral_gff=""    #full path to gff for the ancestral species (uncompressed)

#--- annotate or not #
annotate=""  #a string (YES/NO)? if annotation = YES then annotation of the genomes will be performed
             #else gtf and fastafiles are expected and only the paml/ds/genespace etc will be performed

RelatedProt="/path/to/related_protein.fa" #a full path to a set of external protein data (fasta format) for braker

# ---- orthoDB ---- #
orthoDBspecies=" " #one of : "Metazoa" "Vertebrata" "Viridiplantae" "Arthropoda" "Eukaryota" "Fungi" "Alveolata" "Stramenopiles"
#if a species is given then it will be used for annotation (in addition to eventual RNAseq and orthProteins) 

#if annotate = NO then gtf should be provided: 
gtf1=""
gtf2=""

#RNASeq data ?
RNAseq=""   #YES/NO a string stating wether rnaseq data is available or not
RNAseqlist=" " #/full/path/to/rnaseq.list.txt" listing the rnaseq data available

#TE INFO:
TEdatabase="" #a full path to a custom database of species/genus species TE
#NCBI species for de-novo TE:
ncbi_species=""

#BUSCO SPECIFIC ARGUMENTS:
#busco_lineage name:
busco_lineage=""
```

set all variables accordingly, provide full path when needed

list the RNAseq reads, if available, in a file called "rnaseq.list.txt"
this file is as follows:

```shell
cat rnaseq.list_SE.txt
/path/to/data/rnaseq/rnaseq1.R1.fq.gz
/path/to/data/rnaseq/rnaseq2.R1.fq.gz
/path/to/data/rnaseq/rnaseq3.R1.fq.gz
/path/to/data/rnaseq/rnaseq4.R1.fq.gz
```

with PAIRED END:
```sh
cat rnaseq.list_PE.txt
/path/to/data/rnaseq/rnaseq1.R1.fq.gz	/path/to/data/rnaseq/rnaseq1.R2.fq.gz
/path/to/data/rnaseq/rnaseq2.R1.fq.gz   /path/to/data/rnaseq/rnaseq2.R2.fq.gz
/path/to/data/rnaseq/rnaseq3.R1.fq.gz   /path/to/data/rnaseq/rnaseq3.R2.fq.gz
/path/to/data/rnaseq/rnaseq4.R1.fq.gz   /path/to/data/rnaseq/rnaseq4.R2.fq.gz
```

can be created simply by a command like so:   
`readlink -f your_rna_folder > rnaseq.list_SE.txt ` 

with paired-end awk can do the trick of puting reads on two columns: 
`readlink -f your_rna_folder |awk 'ORS=NR%2?FS:RS ' > rnaseq.list.PE.txt `  


