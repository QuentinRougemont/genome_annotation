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

After cloning the pipeline, please, work from within it to preserve the architecture  

We recommend that you clone the pipeline ***for each of your new project*** and work within it.   

Keep all project separated otherwise it will be difficult to recover your results.   


All options and input must be set in the config file: `config/config`

an example is provided [here](https://github.com/QuentinRougemont/genome_annotation/blob/main/example_data/example.config)

## General parameters set in config

Your **input data** may be

\- one genome containing both sex/mating type chromosomes

\- or two haplotypes containing each one of the sex/mating type chromosomes.

\+ RNAseq data for each genome (optional) + a custom TE database (compulsory for annotation) 


==** Warning names** ==
we recommend to use short name for each of your genome assemblies name and avoid any special characters appart from underscore.

**For instance:**

species-1.fasta will not be valid in GeneSpace. => Use **Species1.fasta** instead.

For the chromosome/contig/scaffold ids we recommand a standard naming including the Species name within it without any thing else than alhpanumeric character.

# Example input data:  
   
**Compulsory** 

see example data folder: 

	* genome1: `example_data/genome1.fa.gz`
	* genome2: `example_data/genome2.fa.gz`
    * config file: `example_data/example.config`

**Optional data** 

	* ancestral genome: `example_data/Mlag129A1.fa.gz`
	* ancestral gff: `example_data/Mlag129A1.gtf.gz`
	* RNAseq data:  `example_data/rnaseq.fq.gz`


**TE data (compulsory)** 

	* custom database: `example_data/TE.fa.gz`


**Proteins for annotation (optional)**

	* example_data/prot.fa 
	  note: if no protein data are available we will use orthoDB11  (downloaded automatically)



| option in config | description |
| --- | --- |
| *genome1* | Full path to the assembly of the genome or haplotype you wish to analyse. |
| *haplotype1* | Name of the genome1 |
| \[*genome2*\] | Full path to the assembly of the second haplotype you with to analyse. Only for the case where you have two haplotypes, each containing one of the sex/mating type chromosomes. |
| \[*haplotype2*\] | Compulsory if *genome2* is given. Name of the second haplotype. |

for more details run: 

```shell
./master.sh --help or ./master -h 
```
This will list the different option


You can launch the whole workflow by typing:

```shell
./master.sh -o 1 2>&1 |tee log
```

This will perform steps I to IV as follows:



# STEP I - RNA-seq alignment - TE masking - Gene prediction - Quality assessment

## Input of step I

### Parameters set in config

|     |     |
| --- | --- |
| **option in config** | **description** |
| *RelatedProt* | Full path to a fasta file containing protein sequences to be used for gene prediction. |
| \[*RNAseqlist*\] | Compulsory with option *a*. Full path to a .txt file containing the list of RNA-seq data files. |
| \[*bamlist1*\] | Compulsory with option *b*. Full path to a .txt file containing the list of bam files for *genome1* (alignment of RNA-seq data onto the fasta of *genome 1*). |
| \[*bamlist2*\] | Compulsory with option *b* and if *genome2* is given. Full path to a .txt file containing the list of bam files for genome2 (alignment of RNA-seq data onto the fasta of *genome 2*). |
| \[*orthoDBspecies*\] | "Metazoa" "Vertebrata" "Viridiplantae" "Arthropoda" "Eukaryota" "Fungi" "Alveolata" or "Stramenopiles". Will use a database from **orthoDB** for gene prediction. |
| *fungus* | "YES" or "NO" (default), whether your species is a fungus. |
| *TEdatabase* | Full path to a database of TE for your species/genus, used in TE prediction, in fasta format. |
| *ncbi_species* | Name of the ncbi species, used in TE prediction. ==list available [here](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi)== |
| \[*gtf1*\] | Compulsory with option *d*. Full path to a .gtf file for an existing gene prediction on *genome1*. |
| \[*gtf2*\] | Compulsory with option *d* and if *genome2* is given. Full path to a .gtf file for an  existing gene prediction on *genome2*. |
| *busco_lineage* | Lineage used for **busco** analysis. You can access the list of available lineages by typing `busco --list-dataset`. |
| *interpro* | YES or NO (default), whether **interproscan** quality check of the genome annotation should be performed. Warning: this can take up to several days for a big dataset. |

Note: If option *d* is set, only *gtf1* and *gtf2* (if you use 2 haplotypes) are needed.

==Give examples of files for RNAseqlist and bamlist HERE==

### Options

**(*a*) - Align RNA & Annotate:**  
Will perform alignment of provided RNA-seq data and use it as additional information for genome annotation.

**(*b*) - Annotate, use BAM of RNA**  
Will use provided BAM of already aligned RNA-seq data and use it as additional information for genome annotation.

**(*c*) - Annotate, no RNA**  
Will perform genome annotation without using RNA information.

**(*d*) - Skip**  
If you have already annotated your genome, will use provided gtf for the following steps, effectively skipping genome annotation.

## Operations of step I

With option *d*, all the following operations are skipped.

### 1\. Alignment of RNA-seq data (only with option *a*)

Corresponding script: `00_scripts/launch_rnaseq.sh`

- Reads trimming using **trimmomatic**

The script will detect whether the data is Single-End or Paired-End and launch trimmomatic, then count the number of retained reads.

- Creation of database for **gsnap** using **gmap**
- Alignment using **gsnap**

Corresponding scripts: `00_scripts/03_gsnap_PE.sh` for PE ; `00_scripts/03_gsnap_SE.sh` for SE

- Mapping quality assessment

Sequencing depth and MAPQ along the genome will be computed and plotted. The resulting plots can be found in ==XXX/Depth/== and ==XXX/mapq/== .

==Insert example plot here==

### 2\. TE discovery and masking

Corresponding script: `00_scripts/launch_step05_to_08.sh`

- *De novo* repeat identification using **repeatmodeler** on the input genome(s) to be annotated
- Genome masking using **repeatmasker** with known TE libraries

### 3\. Genome annotation, quality assessment and filtering

Corresponding script: ./00_scripts/06_braker.sh

- Five successive rounds of gene prediction based on protein database, using **braker**
    
- One round of gene prediction using RNA-seq data, using **braker** (only with options *a* and *b*)
    
- Quality assessment and reports production for each round of gene prediction
    

Two tools can be used at this stage for quality assessment:  
\- **Busco** (corresponding script: `00_scripts/07_busco_after_braker.sh`)  
\- **Braker** report on the raw hintsfile  
This report includes number of genes, number of introns per gene, gene support, number of complete genes and various histograms useful for error checking.

- Combination of protein-based and RNA-seq-based gene models using **TSEBRA** (only with options *a* and *b*)

Please read the [TSEBRA manual](https://github.com/Gaius-Augustus/TSEBRA) before running the script.  
The best round of protein-based gene prediction and the RNA-seq-based gene prediction are given as input in TSEBRA.  
Warning: TSEBRA parameters *intron_support* and *stasto_support* are set to 0 in this workflow (default in TSEBRA: 1 and 2 respectively). This means that only overlapping genes between the two gene models will be filtered. You can change this parameter and others to adjust to your desired level of stringency in the TSEBRA config file: `config/default.cfg`

- Final genome annotation reshaping

**For options *a* and *b*:** The final genome annotation  is the output from TSEBRA.  
**For option *c*:** The final genome annotation is the best protein-based braker round, as evaluated with busco.

Corresponding script: `00_scripts/08_braker_reshaping.sh`  
The genes will be renamed to insert the scaffold name for clarity in downstream analyses.  
Because the next steps in the workflow involve single copy ortholog identification, only the longest transcript is kept for each gene.

- Final genome annotation quality assessment

Two more in-depth tools can be used at this stage for quality assessment: (==\+ busco on final genome pred ?==)  
\- **Blast** against **Uniprot**  
If you wish to skip this, comment l.295 of the script `00_scripts/08_braker_reshaping.sh`  
\- **InterProScan** (if option interpro is set to "YES" in the config file and Blast against Uniprot successfully ran)  
This tool is more time-consuming.

# STEP II - Run GeneSpace - Produce synteny plots

## Input of step II

### Parameters set in config ==in yellow parameters to be set==

|     |     |
| --- | --- |
| **option in config** | **description** |
| *scaffolds* | Full path to the list of focal scaffolds (i.e. the scaffolds composing the sex / mating type chromosomes). |
| *==ancestral==* | "chromosome" or "outgroup", whether the sequence used as proxy for the ancestral state is one of the sex / mating type chromosomes or an ougroup provided below. |
| ==\[*ancestral_chromosome_scaffolds*\]== | Compulsory if *ancestral* is set as "chromosome". Full path to the list of scaffolds of the chromosome used as proxy for ancestral state. |
| \[*==outgroup_orthofinder==*\] | Advised if *ancestral* is set as "chromosome". Full path to a list of genomes to be used as outgroups in OrthoFinder only. |
| \[*ancestral_genome*\] | Compulsory if *ancestral* is set as "outgroup". Full path to the genome of the species used as proxy for the ancestral state. |
| \[*ancestral_gff*\] | Compulsory if *ancestral* is set as "outgroup". Full path to the gff (genome annotation) of the species used as proxy for the ancestral state. |
| ==\[*ancestral_outgroup_scaffolds*\]== | Compulsory if *ancestral* is set as "outgroup". Full path to the list of focal scaffolds for the outgroup used as proxy for ancestral state. |

==Give examples of files for scaffolds, ancestral_chromosome_scaffolds, outgroup_orthofinder and ancestral_outgroup_scaffolds ?==

## Operations of step II

### 4a. Minimizer alignment and plots of target region

Corresponding script: `00_scripts/11_run_genesSpace_paml_ideogram.sh`

- Alignment between the two haplotypes using **minimap2**

If you provided as input two haplotypes containing each one of the sex/mating type chromosomes, the whole haplotypes will be aligned.  
If you provided one genome containing both sex/mating type chromosomes, only the corresponding focal scaffolds (as indicated with option *scaffold*) will be aligned.

- Alignment between the two haplotypes and an outgroup genome used as proxy for ancestral state if you have one (option B), using **minimap2**
    
- Construction of whole genome dotplot using **pafR** (only with two haplotypes as input)
    
- Construction of synteny plot on the focal scaffolds using **pafR**
    
ex: minimap based divergence along the mating type chromosomes :

![Fig2.png](https://github.com/QuentinRougemont/genome_annotation/blob/main/pictures/Fig2.png)

ex: minimap based whole genome alignment : 
	
![Fig3.png](https://github.com/QuentinRougemont/genome_annotation/blob/main/pictures/Fig3.png)



### 4b. Ortholog reconstruction

- Launching of **GeneSpace**

In short, this will:  
\- Identify single copy orthologs with **OrthoFinder**  
\- Construct a dotplot and a riparian plot of whole genome (only with two haplotypes as input) \[==GeneSapce==\]  
\- Construct a riparian plot on focal scaffolds \[==GeneSpace==\]  
For more information, consult the [GeneSpace readme](https://github.com/jtlovell/GENESPACE).

ex: Synteny plot from GeneSpace

![Fig4.png](https://github.com/QuentinRougemont/genome_annotation/blob/main/pictures/Fig4.png)


# STEP III - Compute and plot dS - Plot ideogram and rearrangements

## Input of step III

### Parameters set in config

Same as step II (see above).

## Operations of step III

### 5\. Single copy orthologs alignment

Align all coding sequences from the focal scaffolds.

- **TranslatorX**
- **muscle**

### 6\. dS calculation and plotting

- Calculation of dS & dN using **PAML**
    
- Plotting dS values using a custom R script
    

Corresponding script: `00_scripts/Rscripts/03_plot_paml.R`  
dS values are plotted along the focal scaffolds, and, if 2 haplotypes were given as input, along the whole genome.  
The gene order will be that of the genome used as proxy for the ancestral state: either one of the two sex/mating type chromosomes, or an outgroup (see option *ancestral*).  
It is possible to modify the R script to adapt the plotting options to your needs (for instance position and direction of scaffolds).

Ex: Ds plot : 

![Fig5.png](https://github.com/QuentinRougemont/genome_annotation/blob/main/pictures/Fig5.png)


### === - Plot circos (==step III==)

Corresponding script: `00_scripts/Rscripts/05_plot_circos.R`  
Construction of a circos plot of the focal scaffolds, tracing links between their single copy ortholog genes, using **circlize**.  
It is possible to modify the R script to adapt the plotting options to your needs (for instance position and direction of scaffolds).

==insert circos plot here==

The operations automatically performed by the script `master.sh` stop here.

# Step IV
### 1\. Changepoint analyses

Before launching this step, we strongly suggest that you consult the results of the workflow, especially the dS plot. Once you have deciphered clear hypotheses as to whether there are strata on your focal scaffolds, and where they occur, you can use the R script.
`00_scripts/Rscripts/06.MCP_model_comp.R` to perform changepoint analyses on the dS, using **mcp**.
To that end, you can automatically launch the code ```master.sh -o7``` and it will launch 

==insert graph here and how to interpret==

# Options to run part of the workflow

If you wish to perform only part of the workflow or relaunch it from a given step, use option *\-o*

`bash master.sh -o 2` : perform steps I and II  
now if you have a gtf and and genome assembly (either from running this pipeline or any other annotation tools):
`bash master.sh -o 3` ; perform steps II and III (if step I already ran successfully in a previous run)  
`bash master.sh -o 4` : perform step III only (if steps I and II already ran successfully in a previous run)  
`bash master.sh -o 5` : perform step II only (if step I already ran successfully in a previous run)  
`bash master.sh -o 6` : perform step I only
`bash master.sh -o 7`: perform step IV only




# --------------------------------------------------------------------------

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


