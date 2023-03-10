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


#### minimal braker dependencies:

**Protint** 

wget https://github.com/gatech-genemark/ProtHint/releases/download/v2.6.0/ProtHint-2.6.0.tar.gz 

then add to ~/.bashrc

**Diamond**

#wget https://github.com/bbuchfink/diamond/releases/download/v2.1.1/diamond-linux64.tar.gz

wget wget https://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz

then add to ~/.bashrc

**cdbfasta**:

```
git clone https://github.com/gpertea/cdbfasta.git`
cd cdbfasta
make all 
```


**Augustus**

see details at https://github.com/Gaius-Augustus/Augustus
with a bit of luck this could work:

```
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus
make augustus
```

For augustus I also had to install bamtools: 
bamtools : https://github.com/pezmaster31/bamtools
some change may be necessary to install bamtools properly....

** boost ** was also necessary on some cluster ....


the rest was OK but you may need to carrefully read the Readme of Augustus to make it work


### add augsutus script/config/bin to the .bashrc ###

```
export AUGUSTUS_CONFIG_PATH=/home/path/to/augustus/config
export AUGUSTUS_BIN_PATH=/home/path/to/augustus/bin/
export AUGUSTUS_SCRIPTS_PATH=/home/path/to/augustus_scripts
```


**genemark** 

wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_pxuuc/gmes_linux_64.tar.gz

the gm_key is necessary
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_pxuuc/gm_key_64.gz

it must be copied to your home
```sh
gunzip gm_key_64.gz 
mv gm_key_64 ~/.gm_key
```

Note: you must register [online](http://exon.gatech.edu/GeneMark/license_download.cgi) for genemark


#Other dependencies: 

**TSEBRA** available [here](https://github.com/Gaius-Augustus/TSEBRA)

**trimmomatic** software available [here](http://www.usadellab.org/cms/?page=trimmomatic)

**GSNAP** for read mappig available [here](http://research-pub.gene.com/gmap/)

[repeatmodeler](https://www.repeatmasker.org/RepeatModeler/) and [repeatmasker](https://www.repeatmasker.org/)

**ffread** to extract CDS from fasta [see](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)

**transeq**  from EMBOSS to convert fasta into protein [click_here](https://www.bioinformatics.nl/cgi-bin/emboss/help/transeq)

**BUSCO** for quality assesment (https://busco.ezlab.org/)


**optional:**


**[agat](https://agat.readthedocs.io/en/latest/index.html)** for statistics, quality assesment, errors

**interproscan** for [annotation](https://interproscan-docs.readthedocs.io/en/latest/index.html):  



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

install all perl-dependencies! 
like so:

```
conda install -c anaconda perl biopython
conda install -c bioconda perl-app-cpanminus perl-hash-merge perl-parallel-forkmanager \
    perl-scalar-util-numeric perl-yaml perl-class-data-inheritable \
    perl-exception-class perl-test-pod perl-file-which  perl-mce \
    perl-threaded perl-list-util perl-math-utils cdbtools \
    perl-list-moreutils

#additional install:
conda install -c bioconda perl-file-homedir
```

install [Augustus](https://github.com/Gaius-Augustus/Augustus) and [Bamtools](https://github.com/pezmaster31/bamtools). 

Note: Without root privilege I had to find some tricks for Bamtools and Augustus. 
In case of bugs with Augustus see details [here](https://github.com/Gaius-Augustus/Augustus/blob/master/docs/INSTALL.md)  

export augustus ```config/bin/script```  path (add them to your ~/.bashrc):  

```
export AUGUSTUS_CONFIG_PATH=/home/path/to/augustus/config
export AUGUSTUS_BIN_PATH=/home/path/to/augustus/bin/
export AUGUSTUS_SCRIPTS_PATH=/home/path/to/augustus_scripts
```




#### when all is ok:

Run: 
```./00_scripts/06_braker.sh 2>&1 |tee braker.log``` 

This will run Braker separately for RNAseq and the protein database.   


I use 5 runs for the proteiinDB and choose the one with best Busco score 



## 7 -  Combining different run with TSEBRA

### /!\ WARNING /!\

read TSEBRA manual before running the script. 
set the parameter of tsebra accordingly

then run:
```sh
./00_scripts/07_tsebra.sh species_name best_database_run 
```

## 8 - Write a report -- quality assesment and extraction of CDS

* When using TSEBRA we no longer have a consensus file for amino-acid (`augustus.hints.aa` nor the coding seq `augustus.hints.codingseq`),   
  so we extract them again from the fasta using gffread and convert them with transeq.
  
See exemple script ```00_scripts/08_extractcds.sh```  to do this

* run busco

* run braker script to obtain a report

* annotate further with [interproscan](https://interproscan-docs.readthedocs.io/en/latest/index.html):  

```
interproscan.sh -i input.prot.fasta -goterms -cpu 16 2>&1 |tee interpro.log
```

Note: I had to install libdw1 (without root). 


# ------   Under Construction ------------ ##
## Running with long-reads PacBio IsoSeq 


Here's some exploratory stuff combining long-reads from PacBio + RNAseq from short + Protein data
I've followed the protocol from braker with some modifications: https://github.com/Gaius-Augustus/BRAKER/blob/master/docs/long_reads/long_read_protocol.md



## ----- dependencies ---- 
### new braker : ??
change braker to long read mode by cloning another version in a separate directory

```
git clone https://github.com/Gaius-Augustus/BRAKER
cd BRAKER/
git checkout long_reads
```

export it to your $PATH

###??new tsebra:

```
git clone https://github.com/Gaius-Augustus/TSEBRA


cd TSEBRA/
git checkout long_reads
```

export it to your $PATH 


###??GeneMark ST: 

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

For species different from the reference genome I've added the `asm20`??parameter in minimap. Need to check if this affect the mapping of such reads

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

then use RNAseq (see script ```00_scripts/06_braker.sh```??)


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

