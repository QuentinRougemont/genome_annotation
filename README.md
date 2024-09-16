# genome annotation - Synteny - Ds computation - changepoint analysis - whole genome alignments 
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
I - Perform genome annotation using braker3

II - Identify synteny blocks and rearragnements (GeneSpace and Minimap2)

III - Plot Ds along the genome using ancestral genomes

IV - Perform changepoint analysis to obbjectively identify evolutionary strata

<img src="https://github.com/QuentinRougemont/genome_annotation/blob/main/pictures/Fig1.png" width = "490" heigth = "490">


# Full automated installation 


### Mamba

If you want to avoid potential conflicting versions or do not have root access on your device, you can use **conda** or **mamba** to install dependencies.

We recommend mamba for linux:

```
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge-pypy3-Linux-x86_64.sh
bash Miniforge-pypy3-Linux-x86_64.sh
#note: see here for other verions of mamba: https://github.com/conda-forge/miniforge
```

### Cloning the git

To get the workflow on your device:

```bash
git clone https://github.com/QuentinRougemont/genome_annotation
```

then run the following:
```sh
mamba env create annotation.yml  

#for busco:
mamba env create busco_env.yml

#for repeatmodeller:
mamba env repeatmodeler.yml

#and for non-conda dependencies:
bash ./dependencies.sh
```

### Manual installation of each software separately:

if you prefer, you can install all dependencies one by one, this takes more time but may help solving bugs. See [here](https://github.com/QuentinRougemont/genome_annotation/blob/main/.infos/install.md) to have a list of all dependencies.

# Launching the whole workflow

All options and input must be set in the config file: `config/config`

## General parameters set in config

Your **input data** may be

\- one genome containing both sex/mating type chromosomes

\- or two haplotypes containing each one of the sex/mating type chromosomes.

\+ RNAseq data for each genome (optional) + a custom TE database (compulsory for annotation) 


==Warning names XXX==
we recommend to use short name for each of your genome assemblies name and avoid any special characters appart from underscore.

**For instance:**

species-1.fasta will not be valid in GeneSpace. => Use **Species1.fasta** instead.

For the chromosome/contig/scaffold ids we recommand a standard naming including the Species name within it without any thing else than alhpanumeric character.



| option in config | description |
| --- | --- |
| *genome1* | Full path to the assembly of the genome or haplotype you wish to analyse. |
| *haplotype1* | Name of the genome1 |
| \[*genome2*\] | Full path to the assembly of the second haplotype you with to analyse. Only for the case where you have two haplotypes, each containing one of the sex/mating type chromosomes. |
| \[*haplotype2*\] | Compulsory if *genome2* is given. Name of the second haplotype. |

You can launch the whole workflow by typing:

`bash master.sh -o 1`

This will perform steps I to IV as follows:


# HOW TO USE:

After cloning the pipeline, please, work from within it to preserve the architecture  

We recommend that you clone the pipeline ***for each of your new project*** and work within it.   

Keep all project separated otherwise it will be difficult to recover your results.   


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
 * synteny, arrangement and Ds inference (if your already have annotations for your genomes) 
 

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


# list of operations and tools


| __Operation__                     |  __Tools__                         |  __data type__  | 
|:---------------------------------:|:------------------------------:|:-----------:| 
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


### Step by step guide: 

# --------------------------------------------------------------------------

All steps will be performed automatically by the ```master.sh``` script once the config file is set appropriately


## 0 - rename the contigs/scaffold/chromosome ID in your reference genome. 

We provide a script: 
```00_scripts/00_rename_fasta.py```

that will rename automatically the contigs_id if provided the current name and a new name of the following form:

``` [species]_[haplotype]```  

to create something like:

``` [species]_[haplotype]_[contigID] ```  

We recommand short name with simple separator such as ```_```      

Note: Long name tend to cause bug with paml.

In the rest of the pipeline, we will insert the contigs/scaffold ID in the gene ID, so they are easier to track.


## With RNA seq data: 

the script: 
```00_scripts/launch_rnaseq.sh```

will be launch automatically  

the following Steps will be performed:  


## 1 - trim the reads: trimmomatic 

the script will launch trimmomatic and recognized wether your data are Sinle-End or Paired-End
at the end the following script is launch to count the number of retained reads:

```sh
/00_scripts/utility_scripts/count_read_fastq.sh
```

## 2 - create database for gsnap

the following steps are performed automatically for each of your genome:

```sh
cd haplo1
./00_scripts/02_gmap.sh 03_genome/your_genome.fa.gz
cd ../haplo2
./00_scripts/02_gmap.sh 03_genome/your_genome.fa.gz
cd ../
 ```

## 3 - alignment with gsnap:

for a given genome located in the folder `03_genome` and a set of input in `02_trimmed` ;  

a script will launch **gsnap** automatically:  

* If PE:  ```00_scripts/03_gsnap_PE.sh```  
 
* If SE:  ```00_scripts/03_gsnap_SE.sh```  


## 4 - mapping quality assessment 

for each RNAseq data we will compute the sequencing depth and MAPQ along the genome and this will be plotted automatically.
This is implemented directly in the script ```03_gsnap_[P/S]E.sh```  

Insert example plot here


## 5 - TE discovery and masking


## TO UPDATE

This part will run the script ```launch_step05_to_08.sh```  

It will perform the following :  

   
* step 05:  

    * run **repeatmodeler**: denovo repeat identification from your genome   

    * run **repeatmasker** : mask the genome using known TE Libraries 



**data needed :** 

    * your genome  

    * NCBI database and species name on NCBI for TE 

    *Ideally: a custom (home made) TE library 
 
* step 06: running braker

* step 07: assessing quality 

* step 08: filtering and rehaping the data



## 6 - Runnning braker

**data needed :**
 
* Your genome 

* A busco lineage name 

At least one of the three following data must be provided : 

* 1 - RNAseq for the target species

and/or:

* 2 - An orthoDB species name to download external protein data

and/or: 

* 3 - A custom database of proteins from closely related species


### WARNING /!\ 

**all details should be provided in the config file** 

make sure to have all the dependencies installed as indicated on [braker](https://github.com/Gaius-Augustus/BRAKER#installation)

#### when all is ok:


The master script will run automatically:

```./00_scripts/06_braker.sh 2>&1 |tee braker.log``` 

This will run Braker separately for RNAseq and the protein database.   

I use 5 runs for the proteiinDB and choose the one with best Busco score  
 

## 7 - Evaluate quality  

4 tools can be used for quality assesment :

1 - Busco 

2 - Braker report on the raw hintsfile

3 - Blast against Uniprot (optional)

4 - UniProt (optional - more time consuming)

**1 - Busco**  
 
This will run the script 
```
00_scripts/07_busco_after_braker.sh
```

to assess the quality of each braker annotation.  


Again you need to provide the busco species name in the config file   



**2 - Braker Report :**  



On all hintsfile this will generate very usefull report including:  
 
	* Number of  Gene (Total, Single-exons and Multi-exons genes) 
 
	* Number of introns per genes 

	* support for the genes (count and %) 
 
	* complete genes (count and %)  
 
	* as well as various histograms usefull for error checking


**3 - Blast against Uniprot :** 

this will be performed automatically -  

comment the line 295 of the script in ```00_scripts/08_braker_reshaping.sh``` if you don't want to  


**4 - InterProsScan annotations :**  


to obtain InterProScan annotation set ```interpro="YES"``` in the config/config file  


This will run : 

```
interproscan.sh -i input.prot.fasta -goterms -cpu 16 2>&1 |tee interpro.log
```



Note that if you disabled blast against uniprot this will not be performed either 

Interproscan and blast against uniprot will in reality be run at the very end of the renaming process when single transcript have been selected 

 

## 8 - reshape the data : combine run with TSEBRA - rename genes 


Read TSEBRA manual before running the script. 

set the parameter of tsebra according to the levels of stringeancy that you want .
 
we provide a config file in ```config/default.cfg``` 

the levels of stringeancy have been reduced to keep genes of interest that can be highly degenerated


then our pipeline will automatically run: 

```sh
./00_scripts/08_braker_reshaping.sh -s species_name -r YES/NO
```

internally the gene will be renamed to insert the scaffold name in it. 

the longest transcript will further be kept as this is important for single copy orthologs identifications 
  


## 9 - Run GeneSpace - compute and plot Ds - plot ideogram - plot arrangement 


* **input needed:**  gene annotation for haplo1/haplo2 + the ancestral genome 

This will: 
* 1 - Identify single copy orthologs (OrthoFinder run from GeneSpace)  

* 2 - construct dotplot and riparian plot of whole genome (GeneSpace) 

* 3 - perform a riparian plot focussed on the sex chromosome (GeneSpace) 


* 4 - Align all CDS from the hap1 vs hap2 sex chromosomes (TranslatorX + muscle)

* 5 - Compute Ds & Dn from PAML

* 6 - Plot Ds values along the focal species (either ancestral genome or haplotype1) and gene rank order to display rearrangements 


insert some geneSpace plot here 


## 10 - Plot circos

to do:

* fix the code 

* insert circos here 

## 11 - perform changepoint analyses

to do

insert graph here and how to interpret  

explain the output

 
 

## 12 - minimap alignements and plots of target region


This is the easiest step implemented in : 

```
00_scripts/11_run_geneSapce_paml_ideogram.sh
```  
 

It will: 

* run **minimap**  between the two haplotype 

* run **minimap** between the two haplotype and ancestral genome if you have one  

* run pafR to construct whole genome dotplot and synteny plot on the focal scaffold (X/Y etc).


insert some graph here





