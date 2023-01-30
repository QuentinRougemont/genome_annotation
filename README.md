# genome_annotation
genome annotation (RNAseq and TE) pipeline (draft)

# Purpose:
##  sets of scripts to perform genome annotation using braker2

in details these scripts will: 
* trim RNAseq reads
* map them to the reference genome
* mask the genome using repeatModeller + external evidence
* run braker on RNAseq + external protein data
* combine braker results with TSEBRA
* reshape the output

for now it is designed for microbotryum as some custom fungus database are necessary but it should work on any organism with RNAseq data



## Dependencies: 

**braker2** and all of its dependencies available [here](https://github.com/Gaius-Augustus/BRAKER)

**TSEBRA** available [here](https://github.com/Gaius-Augustus/TSEBRA)

**trimmomatic** software available [here](http://www.usadellab.org/cms/?page=trimmomatic)

**GSNAP** for read mappig available [here](http://research-pub.gene.com/gmap/)

[repeatmodeler](https://www.repeatmasker.org/RepeatModeler/) and [repeatmasker](https://www.repeatmasker.org/)

**ffread** to extract CDS from fasta [see](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)

**transeq** to convert fasta into protein [click_here](https://www.bioinformatics.nl/cgi-bin/emboss/help/transeq)

**optional: BUSCO** for quality assesment


# Steps 

your genome should be in the `03_genome` folder!  


## 1 - trimmomatic: trim the reads  

you can count the number of retained reads using the scripts : `./00_scripts/utility_scripts/count_read_fastq.sh`

for instance: 
```sh
cd 02_trimmed ../00_scripts/utility_scripts/count_read_fastq.sh *fq >> read_count.txt
```


## 2 - create database for gsnap

```sh
./00_scripts/01_gmap.sh 03_genome/your_genome.fa
 ```

Note: only works with uncompressed genome  
 

## 3 - alignment with gsnap:

for a given genome located in the folder `03_genome` and a set of input in `02_trimmed` ;  
simply loop over files:

```sh
for i in 02_trimmed/*R1.paired.fastq.gz ; 
do 
./00_scripts/02_gsnap.sh 03_genome/M_inter_1389.PBcR.20160424.quiver.finished.fasta $i ; 
done
```

## 4 - count the number of well mapped reads

use the script:
```
./00_scripts/03_count_mapped_read.sh 
```

and compare it to the number of trimmed reads to evaluate the quality of the data


## 5 - TE discovery and masking

We will identify denovo repeat (from our genome) using *RepeatModeler* and mask the genome using known TE Library using *Repeatmasker* 

use this script:

```
./00_scripts/05_repeatmodeler.sh
```

### WARNING:  
edit the script to provide path to your own database of repetetive regions. 
Here I use 3 custom libraries + online data you may have more or less of these, so comment or delete unessecary rounds of RepeatMasker 


## 6 - Runnning braker

## data: 
* RNAseq for the target species

* Protein database from several closely related species


The script ```./00_scripts/06_braker.sh``` will run Braker separately for RNAseq and the protein database. 
I use 5 runs for the proteiinDB and choose the one with best Busco score 



## 7 -  Combining different run with TSEBRA

### /!\ WARNING /!\

read TSEBRA manual before running the script. 
set the parameter of tsebra accordingly

then run:
```sh
./00_scripts/07_tsebra.sh species_name best_database_run 
``

## 8 - Write a report -- quality assesment and extraction of CDS

run busco

run braker script to obtain a report

annotate further with interproscan
