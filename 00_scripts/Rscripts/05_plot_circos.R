#!/usr/bin/env Rscript

#Purpose:  Script to plot circos graphs representing synteny between a focus haplotype and its inferred reference state
#Author: LB
#Modified: QR
#Date: 28-11-23 - 12-11-24

#INPUT: 
# 1 - Table of scaffolds of interest (i.e. scaffolds belonging to sex chromosomes) for the haplotype and the reference.
#     2 columns: Scaffold name, Species name
#     ! Note that the order of scaffolds in the file determines the order in which they are displayed on the plot.
#     ! Note that a third, optional column can indicate whether the scaffolds have to be reversed (1) on the plot
#       or not (0).
# 2 - BED file of genes of interest to be displayed on the plot
# 3 - Table of synteny data (will be read automatically if previous steps were sucessfull)
# 4 - .fai index of the two genomes (will be read automatically if previous steps were sucessfull)

####-------------------------- INITIALISATION ------------------------------####

#------------- check that libraries are installed and load them ---------------#
packages <- c('circlize','dplyr','tidyr','wesanderson','magrittr','optparse')

#install.packages(setdiff(packages, rownames(installed.packages())))
install.packages(setdiff(packages, rownames(installed.packages())), repos="https://cloud.r-project.org" )
invisible(lapply(packages, suppressMessages(suppressWarnings(suppressPackageStartupMessages(library))), character.only = TRUE))

#------------- read input from the command line -------------------------------#
option_list <- list(
  make_option(c("-s","--species1"), type="character", default=NULL,
              help="species1 name' [default %default]", 
              ),
  make_option(c("-p","--species2"), type="character", default=NULL,
              help="species2 [default %default]",
              ),
  make_option(c("-c","--chromosome_file"), type="character", default=NULL,
              help="txt file of target chromosomes, a 2 column file with species id\tchromosome id [default %default]",
              ),
  make_option(c("-y","--synteny_table"), type="character", default=NULL,
              help="txt file of synteny table (generated from previous steps) [default %default]",
              dest="synteny_file"),
  make_option(c("-f","--fai_species1"), type="character", default=NULL,
              help="samtools index file from species1 (generated from previous steps) [default %default]",
              dest="fai1"),
  make_option(c("-g","--fai_species2"), type="character", default=NULL,
              help="samtools index file from species2 (generated from previous steps) [default %default]",
              dest="fai2"),
  make_option(c("-i","--gene_species1"), type="character", default=NULL,
              help="bed file of genes for species1 (generated from previous steps) [default %default]",
              dest="g1"),
  make_option(c("-j","--gene_species2"), type="character", default=NULL,
              help="bed file of genes for species2 (generated from previous steps) [default %default]",
              dest="g2"),
  make_option(c("-t","--TE_species1"), type="character", default=NULL,
              help="bed file of TE for species1 (generated from previous steps) [default %default]",
              dest="TE1"),
  make_option(c("-u","--TE_species2"), type="character", default=NULL,
              help="bed file of TE for species2 (generated from previous steps) [default %default]",
              dest="TE2"),
  make_option(c("-l","--links"), type="character", default=NULL,
              help="bed file of regions to highlight (gene/centromere/etc) [default %default]",
              ),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)

options(error=traceback)

parser <- OptionParser(usage = "%prog -s species1 -p species2 -c chromosomes -y synteny_table
                       -f fai_specie1 -g fai_species2 [options]", 
                       option_list = option_list)
opt = parse_args(parser)
#opt = parse_args(OptionParser(option_list=option_list))

if(opt$v){
  cat(paste0("script to perform CIRCOS plot\n\n"))
  cat(paste0("COMPULSORY PARAMETERS :\n\t--species1 (-s): ", opt$species1,"\n"))
  cat(paste0("\t--species2 (-p): ", opt$species2,"\n"))
  cat(paste0("\t--chromosome (-c): ", opt$chromosome_file,"\n"))
  cat(paste0("\t--synteny_table (-y): ", opt$synteny_file,"\n"))
  cat(paste0("\t--fai_species1: index file sp1 (-f): ", opt$fai1,"\n"))
  cat(paste0("\t--fai_species2: index file sp2 (-g): ", opt$fai2,"\n"))
  cat(paste0("optional parameters: \n"))
  cat(paste0("\t--gene_species1 (-i bed file of gene for sp1): ", opt$g1,"\n"))
  cat(paste0("\t--gene_species2 (-j bed file of gene for sp2): ", opt$g2,"\n"))
  cat(paste0("\t--TE_species1 (-t): bed file of TE for sp1", opt$TE1,"\n"))
  cat(paste0("\t--TE_species2 (-u): bed for of TE for sp2", opt$TE2,"\n"))
  cat(paste0("\t--links (-u): bed file of links  to highlight", opt$links,"\n"))

}

#----------- load parameters --------------------------------------------------#
reference <- opt$species1
haplo <- opt$species2 

writeLines("\n~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
writeLines(paste0("\nreference is ", reference,"\n"))
writeLines(paste0("haplo is ", haplo, "\n"))
writeLines("~~~~~~~~~~~~~~~~~~~~~~~~~~\n")


#------------- Import other files ---------------------------------------------#
# import synteny data
writeLines("\n~~~~~~ loading data ~~~~~~~\n")
syn <- read.table(opt$synteny_file, header=T, as.is=T, sep='\t')

# import contig informations from .fai index
index_ref <- read.table(opt$fai1, as.is = T, sep = '\t')[,c(1,2)] %>% 
    set_colnames(., c("chr","end"))

#to fix:
index_hap <- read.table(opt$fai2, as.is = T, sep = '\t')[,c(1,2)] %>% 
    set_colnames(., c("chr","end"))
writeLines("\n~~~~~~ data loaded ~~~~~~~\n")

#import chromosme data: 
chromosomes <- read.table(opt$chromosome_file)
# Check if some contigs have to be inverted
if(ncol(chromosomes)==2) {
  chromosomes$inv=0
}
colnames(chromosomes) <- c("species","chr","inv")



####------------------------ PREPARE CIRCOS DATA ---------------------------####
#------------- Prepare data sets ----------------------------------------------#

# Get list of focus scaffolds for each species
chr_ref <- chromosomes[which(chromosomes$species == reference),]
chr_hap <- chromosomes[which(chromosomes$species == haplo),]
chromosomes <- chromosomes[(chromosomes$species==reference) | (chromosomes$species==haplo),]

# Subset data sets according to the focus contigs
syn <- syn[which(syn$chrom1 %in% chr_ref$chr & syn$chrom2 %in% chr_hap$chr),]
index_ref <- index_ref[match(chr_ref$chr,index_ref$chr),]
index_hap <- index_hap[match(chr_hap$chr,index_hap$chr),]

# Make a unique contig info table
index_ref$species <- reference
index_hap$species <- haplo
contigs <- rbind.data.frame(index_ref,index_hap)

# Make the database for the genomic links
nuc1 <- select(syn, c('chrom1', 'start1','end1'))
nuc2 <- select(syn, c('chrom2', 'start2','end2'))

# Invert contig orientation if needed
to_inv <- chromosomes[which(chromosomes$inv == 1),'chr']

for (contig in to_inv) {
  end=contigs[which(contigs$chr==contig),]$end
  #change the coordinates for the synteny databases
  nuc1[which(nuc1$chrom1==contig),]$start1=end-nuc1[which(nuc1$chrom1==contig),]$start1
  nuc1[which(nuc1$chrom1==contig),]$end1=end-nuc1[which(nuc1$chrom1==contig),]$end1
  nuc2[which(nuc2$chrom2==contig),]$start2=end-nuc2[which(nuc2$chrom2==contig),]$start2
  nuc2[which(nuc2$chrom2==contig),]$end2=end-nuc2[which(nuc2$chrom2==contig),]$end2
}

#make a matrix with the contigs start and end to initialize the circos
nb_contig <- nrow(contigs)
m <- matrix(c(rep(0, nb_contig), c(contigs$end)), ncol=2)

#---------- Optional: getting gene density from bed ---------------------#
if(!(is.null(opt$g1))){
  gedens1 <- read.table(opt$g1)
}

if(!(is.null(opt$g2))){
  gedens2 <- read.table(opt$g2) 
  
  genedensity <- rbind(gedens1,gedens2) %>% 
    select(-V4) %>%
    set_colnames(.,c("chr","start","end")) %>%
    mutate(value = 1) %>%
    filter(chr %in% contigs$chr)

writeLines("\n~~~~~~~~~~~~~~~~~~~~~~~~~~")
cat(paste0("number of gene in sex chr :", nrow(genedensity),"\n"))

}

#---------- Optional: TE density from bed ---------------------#

if(!(is.null(opt$TE1))){
  TE1 <- read.table(opt$TE1)
}

if(!(is.null(opt$TE2))){
  TE2 <- read.table(opt$TE2) 

TEdensity <- rbind(TE1,TE2) %>% 
    select(-V4) %>%
    set_colnames(.,c("chr","start","end")) %>%
    mutate(value = 1) %>%
    filter(chr %in% contigs$chr)

writeLines("\n~~~~~~~~~~~~~~~~~~~~~~~~~~")
cat(paste0("number of TE in sex chr: ", nrow(TEdensity),"\n\n"))

}

#---------- Optional: Prepare table of genes to highlight ---------------------#
#Import the gene positions
if(!(is.null(opt$links))) {
data_genes <- read.table(opt$links, as.is=T, sep='\t') 
  colnames(data_genes)=c("chr","start","end","category")
  #Invert contig orientation if needed
  for (contig in to_inv)
  {
    end=contigs[which(contigs$chr==contig),]$end
    data_genes[which(data_genes$chr==contig),]$start=end-data_genes[which(data_genes$chr==contig),]$start
    data_genes[which(data_genes$chr==contig),]$end=end-data_genes[which(data_genes$chr==contig),]$end
  }

#keep target only: 
data_genes <- data_genes %>% filter(chr %in% contigs$chr)
print(data_genes)
nrow(data_genes)
}
####------------------------ LAUNCH CIRCOS ---------------------------------####
#------------- Define plotting parameters -------------------------------------#

writeLines("\n~~~~~~ preparing colors ~~~~~~~\n")
# Contig colors
col_ref <- "grey"
col_hap <- "grey95"
contig_color <- c(rep(col_ref,nrow(chr_ref)),rep(col_hap,nrow(chr_hap)))

list_cont=chr_ref$chr
rcols=vector(length=nrow(syn))
for(i in 1:nrow(index_ref)) {
  c=list_cont[i]
#  rcols[which(syn$chrom1==c)]=terrain.colors(length(list_cont))[i]
  rcols[which(syn$chrom1==c)]=wes_palette("Zissou1", length(list_cont), type = c("continuous"))[i]
}

#create dir if not present:
if (!dir.exists("02_results/circos")){
  dir.create("02_results/circos")
}


# Output in pdf
pdf(file = paste0('02_results/circos/circos_',haplo,'_on_',reference,'.pdf'))

#------------- Initialize circos ----------------------------------------------#
# Initialization
circos.clear()
circos.par("track.height" = 0.8, 
           "canvas.xlim" = c(-1.1,1.1),
           "canvas.ylim" = c(-1.1,1.1),
           gap.degree = 5, 
           cell.padding = c(0, 0, 0, 0))

circos.initialize(factors=contigs$chr,xlim=m)

# Make contig track
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), max(ylim)+1, chr, cex=0.5, col='black', 
              facing="bending.inside", niceFacing=TRUE)
}, bg.col=contig_color, bg.border="grey40", track.height=0.06)


# Define ticks as megabases
max_size=ceiling(max(c(index_ref$end/10^6,index_hap$end/10^6)))
brk <- seq(0,max_size,1)*10^6

# Trace x axis
circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
  circos.axis(h="top", major.at=brk, labels=round(brk/10^6, 1), labels.cex=0.4,
              col="grey40", labels.col="black", lwd=0.7, labels.facing="clockwise")
}, bg.border=F)

#gene density plots:
if(exists('genedensity')){
  circos.genomicDensity(genedensity, 
                      col="slategray1", 
                      bg.border=contig_color, 
                      #bg.lwd=0.8,
                      window.size = 8000, #force.ylim=FALSE, ylim=c(-0.01,0.99),
                      overlap = FALSE,
                      count_by = "number",  track.height=0.05)
}

#gene density plots:
if(exists('TEdensity')){
  circos.genomicDensity(TEdensity, 
                      col="springgreen", 
                      bg.border=contig_color, 
                      #bg.lwd=0.8,
                      window.size = 8000, #force.ylim=FALSE, ylim=c(-0.01,0.99),
                      overlap = FALSE,
                      count_by = "number",  track.height=0.05)
}

#------------- Plot links -----------------------------------------------------#
#rcols=scales::alpha(ifelse(d$chrom_lag=='MC03',"blue","purple"),alpha=1)
# Color the links according to reference haplotype contigs
# plot links
circos.genomicLink(nuc1, nuc2, col=rcols, border=NA)

#---------- Optional: highlight genes -----------------------------------------#
if(exists('data_genes')){
  writeLines("\nadding links to highlight some regions of interest\n")
  # Make a new track
  circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
    chr=CELL_META$sector.index
    xlim=CELL_META$xlim
    ylim=CELL_META$ylim
  }, bg.col="white", bg.border="grey40", track.height=0.05)
  
  # Plot genes
  list_cat=unique(data_genes$category)

  for(c in 1:length(list_cat)) {
    cat=list_cat[c]
    d_cat=data_genes[which(data_genes$category==cat),]
    col=rainbow(length(list_cat))[c]
    #print(col)
    for(i in d_cat$chr) {
      circos.genomicRect(d_cat[which(d_cat$chr==i),], sector.index=i,
		track.index=1, ytop = 1, ybottom = 0,col=col,border=col)}
    }
}

#
#----------- Write pdf file ---------------------------------------------------#
dev.off()

################################################################################
