# genome_annotation - GeneSpace - dotplot construction -  Ds computation between haplotypes - changepoint analysis - whole genome alignements 


# Purpose:
##  sets of scripts to : 
 1 - perform de-novo TE detection, 

 2 - genome annotation using braker2  

 3 - run GeneSpace  

 4 - plot Ds along the genome using ancestral genomes  

 5 - identify rearrangements  

 6 - plot ideograms   
 
 7 - perform changepoint analyses

in details these scripts will: 
* trim RNAseq reads
* map them to the reference genome
* mask the genome using repeatModeller + external evidence
* run braker on RNAseq + external protein data
* combine braker results with TSEBRA
* reshape the output

* Run GeneSpace (with orthoFinder and MCScanX)
* Compute Ds between single copy orthologs of interest assuming a 1:1:1 correspondance
* Plot Ds along the genome and genes order
* Perform changepoint analysis
* Plot Ideograms 

This code has been tested with linux. 


**Braker3** does not produced results of as good quality as braker for now, but could be used due to the simplicity of implementation through **singularity** 

# Dependencies: 

**braker2** and all of its dependencies available [here](https://github.com/Gaius-Augustus/BRAKER)

basic requirements: **git**, **gcc** , **R**, **make**, **cmake**, **wget** or **curl** , **java** 

conda or **mamba** can be used to install several dependencies, especially without root access 

### First thing is to get mamba or conda. 

I recommand mamba for linux: 

```
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge-pypy3-Linux-x86_64.sh
bash Miniforge-pypy3-Linux-x86_64.sh
#note: see here for other verions of mamba: https://github.com/conda-forge/miniforge
```

### minimal dependencies:  

see [here](https://github.com/QuentinRougemont/genome_annotation/blob/main/.infos/install.md) to have a list of all dependencies or run directly : 

```shell 
git clone https://github.com/QuentinRougemont/genome_annotation
cd genome annotation
./requirements.sh #will work with mamba. alternatively replace by conda (sed -i 's/mamba/conda/g' requirements.sh) 
``` 
 
# minimal input data needed: 

* 2 genomes to annotate :

	* with OR without RNAseq

* OR **2 genomes** already annotated with their corresponding **gtf** files

these must correspond to each haplotype you'd like to compare

* If annotation is required then the following data are recquired:

	* a lineage name for busco evaluation of the genome annotation quality.  

The list of lineage can be obtained [here](https://busco-data.ezlab.org/v5/data/lineages/) or using :

```shell
busco --list-dataset
``` 
in the config/config file set the busco_lineage name that is closer to your study organism

	* a list of protein from the same or closely related species in fasta format

	* if possible a custom database of TE


* optional:  1 ancestral genome with its annotation in gff format.  

	     keep a single transcript per gene.  

### note on input naming: 

we recommend to use short name for each haplotype and avoid any special characters appart from underscore.  

# Quick start:

***before running the pipeline make sure to have ALL dependencies installed***

***make sure to have the correct input data***

### step1 - edit the config file and set the correct path 

the config file is here: ./confi/config

the file is as follows:

```shell
cat config/config
# config file
#--- COMPULSORY MINIMAL LEVEL OF INFORMATION REQUIRED -----
genome1=""     #full path to current genome1 assembly (haplotype1 - compressed or not)
haplotype1=""  #name1 [name of haplotype1 - will be used to rename the genome and the contigs inside the genome]

#----- optional --------
ancestral_genome=""
genome2=""     #full path to current genome2 assembly (haplotype2)
haplotype2=""  #name2 [name of haplotype2 - will be used to rename the genome and the contigs inside the genome]

#--- annotate or not #
annotate=""  #a string (YES/NO)? if annotation = YES then annotation of the genomes will be performed
             #else gtf and fastafiles are expected and only the paml/ds/genespace etc will be performed

RelatedProt="/path/to/related_protein.fa" #a full path to a set of external protein data (fasta format) for braker


#if annotate = NO then gtf should be provided: 
gtf1=""
gtf2=""

#RNASeq data ?
RNAseq=""   #YES/NO a string stating wether rnaseq data is available or not
RNAseqlist=" "

#TE INFO:
TEdatabase="" #a full path to a custom database of species/genus species TE
#NCBI species for de-novo TE:
ncbi_species=""

#BUSCO SPECIFIC ARGUMENTS:
#busco_lineage name:
busco_lineage=""
```

set all variables accordingly, provide full path when needed

list the RNAseq reads if available in a file called "rnaseq.list.txt"
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


## step2 - run the master script 

once the config file is ready with your path and dataset correctly simply run 

```shell
./master.sh 2>&1 |tee log
```

this script should handle automatically the different use cases (annotation only, annotation and plot, etc):

for more details run: 

```shell
./master.sh --help or ./master -h
```





# detailed step by step guide: ------------------------------------------------------------------#
###  Steps 

your genomes should be in the `haplo1/03_genome` and `haplo2/03_genome` folders for each haplotype!  

## 0 - rename the contigs/scaffold/chromosome ID in your reference genome. 

We recommand to use very short name with simple separator such as `_`. 
Long name tend to cause bug with paml.
In the rest of the pipeline, we will insert the contigs/scaffold ID in the gene ID, so they are easier to track.


## 1 - trimmomatic: trim the reads  

you can count the number of retained reads using the scripts : `./00_scripts/utility_scripts/count_read_fastq.sh`

for instance: 
```sh
cd 02_trimmed ../00_scripts/utility_scripts/count_read_fastq.sh *fq >> read_count.txt
```


## 2 - create database for gsnap

```sh
cd haplo1
./00_scripts/02_gmap.sh 03_genome/your_genome.fa.gz
cd ../haplo2
./00_scripts/02_gmap.sh 03_genome/your_genome.fa.gz
cd ../
 ```

## 3 - alignment with gsnap:

for a given genome located in the folder `03_genome` and a set of input in `02_trimmed` ;  
simply loop over files:

```sh
cd haplo1
for i in 02_trimmed/*R1.paired.fastq.gz ; 
do 
./00_scripts/03_gsnap.sh 03_genome/your.genome.fa.gz $i ; 
done

#mapping against haplo2:
cd ../haplo2
for i in 02_trimmed/*R1.paired.fastq.gz ; 
do 
./00_scripts/03_gsnap.sh 03_genome/your.genome.fa.gz $i ; 
done
cd ../

```

#### DEPRECATED STEP (ALREADY DONE IN GSNAP )
## 4 - count the number of well mapped reads
#
#use the script:
#```
#./00_scripts/04_count_mapped_read.sh 
#```
#
#and compare it to the number of trimmed reads to evaluate the quality of the data


## 5 - TE discovery and masking

## Note: It is possible to run step 05-06-07-08 with ```run_step_05_06_07_08.sh``` but read everything before !

We will identify denovo repeat (from our genome) using *RepeatModeler* and mask the genome using known TE Library using *Repeatmasker* 

use this script:

```
./00_scripts/05_repeatmodeler.sh 2>&1 |tee RM.log
```

### WARNING:  
edit the script to provide path to your own database of repetetive regions. 
Here I use 3 custom libraries + online data you may have more or less of these, so comment or delete unessecary rounds of RepeatMasker 


## 6 - Runnning braker

## data: 
* RNAseq for the target species

* Protein database from several closely related species


### WARNING /!\ 

make sure to have all the dependencies installed as indicated on [braker](https://github.com/Gaius-Augustus/BRAKER#installation)

#### when all is ok:

Run: 
```./00_scripts/06_braker.sh 2>&1 |tee braker.log``` 

This will run Braker separately for RNAseq and the protein database.   

I use 5 runs for the proteiinDB and choose the one with best Busco score 


## 7 - Evaluate quality with busco - combine different run with TSEBRA - reshape braker output 

### /!\ WARNING /!\

read TSEBRA manual before running the script. 
set the parameter of tsebra accordingly

then run:
```sh
./00_scripts/08_braker_reshaping.sh -s species_name -r YES/NO

#with -s the haplotype_name 
# -r a YES/NO string stating whether RNEseq was used (YES) or NOT (NO)
 
```

## 8 - Write a report -- quality assesment and extraction of CDS

* annotate further with [interproscan](https://interproscan-docs.readthedocs.io/en/latest/index.html):  

```
interproscan.sh -i input.prot.fasta -goterms -cpu 16 2>&1 |tee interpro.log
```

Note: I had to install libdw1 (without root). 


## 9 - Run GeneSpace - compute and plot Ds - plot ideogram - plot arrangement 

* input needed: haplo1/haplo2 + the ancestral genome 

## 10 - Plot circos

to do 


## 11 - perform changepoint analyses

to o

## 12 - minimap alignements and plots of target region




# DEPRECATED PART FOR NOW: 

# ------   Under Construction ------------ ##
## Running with long-reads PacBio IsoSeq 


Here's some exploratory stuff combining long-reads from PacBio + RNAseq from short + Protein data
I've followed the protocol from braker with some modifications: https://github.com/Gaius-Augustus/BRAKER/blob/master/docs/long_reads/long_read_protocol.md



## ----- dependencies ---- 
### new braker :  
change braker to long read mode by cloning another version in a separate directory

```
git clone https://github.com/Gaius-Augustus/BRAKER
cd BRAKER/
git checkout long_reads
```

export it to your $PATH

### new tsebra:

```
git clone https://github.com/Gaius-Augustus/TSEBRA


cd TSEBRA/
git checkout long_reads
```

export it to your $PATH 


### GeneMark ST: 

download GeneMark ST http://exon.gatech.edu/GeneMark/license_download.cgi 

export ```gmst.pl``` 

### Minimap:

see [here](https://github.com/lh3/minimap2)
```
git clone https://github.com/lh3/minimap2
cd minimap2 && make
```


### Cupcake :

```
mamba create -n anaCogent 
mamba activate anaCogent
mamba install -c anaconda biopython

source activate anaCogent

git clone https://github.com/Magdoll/cDNA_Cupcake.git
cd cDNA_Cupcake
python setup.py build
python setup.py install



## then run the script for IsoSeq: 
```

run scripts:

```
00_scripts/long_reads/isoseqv2_minimap_target_species.sh
```

For species different from the reference genome I've added the `asm20` parameter in minimap. Need to check if this affect the mapping of such reads

```
00_scripts/long_reads/isoseqv2_minimap_asm20.sh
```

the resulting protein from the outgroup were combined to a database of external protein:

```cat 07_isoseq/*mRNA protein.fa >> all.proteins.fa```

Note: Running pbmm2 + Isoseq3 seems to provide fairly similar results



## then run braker on the protein database from external evidence + long read outgroups

```
genome=$1
protein="all.proteins.fa"
wd=08_braker_long_read/protein
mkdir -p $wd
braker.pl --genome=$genome --softmasking --cpu="$NCPUS" --epmode --prot_seq="$protein"  --workingdir="$wd" 2> $wdir/braker2.log
```

then use RNAseq (see script ```00_scripts/06_braker.sh``` )


## Combine dataset with TSEBRA


```sh
name=$1 #specie name
LR=07_isoseq_mapped/
wd=08_braker_long_read/protein

tsebra.py -g 06_braker/rnaseq/augustus.hints.gtf, $wd/augustus.hints.gtf -e 06_braker/rnaseq/hintsfile.gff,$wd/hintsfile.gff -l $LR/gmst."$name".gtf -c long_reads.cfg -o tsebra."$name".gtf
```

## then perform all the quality check as described above for classical data:

1 - busco   

2 - proportion of correct annotation (missing start/stop codon, length of transcript, length of intron, number of exon per genes, etc)  

3 - blast all protein against each other  

4 - blast against uniprot  

5 - annotate with Inter-pro  


## ACKNOWLEDGMENTS:

