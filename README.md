# genome annotation - Synteny - Ds computation between haplotypes - changepoint analysis - whole genome alignements 
====================================================================================

   * [Purpose](#purpose)
   * [Dependencies](#dependencies)
   * [Input data](#input-data)
   * [Example input data](#example-input-data)
   * [Quick start](#quick-start)
   * [Example plot](#example-plot)
   * [Detailed steps](#detailed-steps)


# Purpose:
##  sets of scripts to : 
 1 - perform de-novo TE detection, 

 2 - genome annotation using braker2  

 3 - identify synteny blocks and rearragnements (GeneSpace and Minimap)  

 4 - plot Ds along the genome using ancestral genomes  

 5 - perform changepoint analyses


<img src="https://github.com/QuentinRougemont/genome_annotation/blob/main/pictures/Fig1.png" width = "490" heigth = "490">


**Braker3** does not produced results of as good quality as braker for now, but could be used due to the simplicity of implementation through **singularity** 

# Dependencies: 

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

 
# Input data: 


* minimal input:
	* 2 genomes assemblies (**fasta files**) 

	* additional settings (RNAseq data, etc) can be declared in the **config file** 

more details on the config file are [here](https://github.com/QuentinRougemont/genome_annotation/blob/main/.infos/input_data_info) 

### note on input naming: 

we recommend to use short name for each haplotype and avoid any special characters appart from underscore.  


# Example input data:  

   
**Compulsory** 

	* genome1: `example_data/genome1.fa.gz`
	* genome2: `example_data/genome2.fa.gz`
    * config file: `example_data/example.config`


**Optional data** 

	* ancestral genome: `ancestral_genome/Mlag129A1.fa.gz`
	* ancestral gff: `ancestral_genome/Mlag129A1.gtf.gz`
	* RNAseq data:  `example_data/rnaseq.fq.gz`


**TE data (compulsory)** 

	* custom database: `example_data/TE.fa.gz`


**Proteins for annotation (optional)**

	* example_data/prot.fa 
	  note: if no protein data are available we will use orthoDB11  (downloaded automatically)



# Quick start:

***before running the pipeline make sure to have ALL dependencies installed***

***make sure to have the correct input data***

### step1 - edit the config file and set the correct path 

the config file is here: `./config/config`




## step2 - run the master script 

once the config file is ready with your path and dataset correctly simply run 

```shell
./master.sh 2>&1 |tee log
```

this script should handle automatically the different **use cases**

Use case will be infered from the config file automatically. These can be:
 * annotation only
 * annotation, synteny, arrangement and Ds inference
 * synteny, arrangement and Ds inference  
 

for more details run: 

```shell
./master.sh --help or ./master -h
```



# Example Results 


   * 1 - Genome Annotations 

insert results here 


   * 2 - Genome Annotation quality assesment based on busco: 

insert busco plot here 


   * 3a -  minimap based divergence along the mating type chromosomes :

![Fig2.png](https://github.com/QuentinRougemont/genome_annotation/blob/main/pictures/Fig2.png)

   * 3b - minimap based whole genome alignment : 
	
![Fig3.png](https://github.com/QuentinRougemont/genome_annotation/blob/main/pictures/Fig3.png)



   * 4 - Synteny plot from GeneSpace

![Fig4.png](https://github.com/QuentinRougemont/genome_annotation/blob/main/pictures/Fig4.png)


   * 5 - Ds plot : 

![Fig5.png](https://github.com/QuentinRougemont/genome_annotation/blob/main/pictures/Fig5.png)



   * 6 - changePoint inference : 

insert some plot here 



# detailed steps  

### Step by step guide: 

# --------------------------------------------------------------------------


# list of operations and tools


| Operation                        |  Tools                         |   data type  | 
|:---------------------------------|:-------------------------------|-----------| 
| __read trimming__                |  Trimmomatic                   | RNAseq         | 
| __read mapping__                 |  gmap/gsnap                    | RNAseq          | 
| __sorting read__                 |  samtools                      | RNAseq        |
| __mapping quality assement__     |  samtools + R                  | RNAseq        |
| __TE detection and softmasking__ |  RepeatModeler + RepeadMasker  | genome assembly |
| __genome annotation__            |  Braker + tsebra               | genome assembly + protein + database |
| __quality assessment__           |  Busco + Blast + Inter Pro     | genome prediction |
| __Synteny and HOG__              |  GeneSpace (including OrthoFinder/MCScan) | gene prediction and proteins |
| __cds alignement__               |  muscle + translatorX          | gene prediction (single copy orthologs) | 
| __Ds computation__               |  paml                          | CDS alignment |
| __Ds plot/CIRCOS plot__          |  R                             | Ds and genome information |
| __whole genome alignemnt__       |  minimap2                      | genome assemblies |
| __gene microsynteny__            |  R                      | single copy orthologs |
| __changepoint analysis__         |  R                      | Ds values and gene order |



This code has been tested with linux. 



Normally, you should only run the script ```./master.sh```

below we provided a description of what will be done at each steps.

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

to do


## 12 - minimap alignements and plots of target region  

to do








