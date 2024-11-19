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
packages <- c('circlize','dplyr','tidyr','wesanderson','magrittr')
#install.packages(setdiff(packages, rownames(installed.packages())))
install.packages(setdiff(packages, rownames(installed.packages())), repos="https://cloud.r-project.org" )
#invisible(lapply(packages, library, character.only = TRUE))
invisible(lapply(packages, suppressPackageStartupMessages(library), character.only = TRUE))

#------------- read input from the command line -------------------------------#
args <- commandArgs(T)
#
## test if there are at least 3 arguments: if not, return an error
#if (length(args) < 3) {
#  stop("At least the name of 2 species and a list of focus scaffolds must be supplied", call.=FALSE)
#} else if (length(args)==3) {
#  print("assuming no particular gene to be highlighted")
#} else {
#  print("Genes of interest will be displayed in colors")
#  data_genes=read.table(args[4], as.is=T, sep='\t')
#  # test that bed file has 4 columns
#  if (ncol(data_genes)<4) {
#    stop("Missing one or more columns in BED files. Expected: chr, start, end, state", call.=FALSE)
#  }
#}
## enter variables
haplo <- args[1]
reference <- args[2]
chromosomes <- read.table(args[3])
synt <- args[4] 
fai1 <- args[5]
fai2 <- args[6]

print(paste0("reference is ", reference))
print(paste0("haplo is", haplo))
print(paste0("chromosome are ", args[3]))
print(paste0("synt are", synt))
print(paste0("fai1 is", fai1))
print(paste0("fai2 is", fai2))

####  examples #######
#haplo <- "Mlyc1064a1" #args[1]
#haplo <- "Mlyc1064a2"
#reference <- "ancestral_sp" #args[2]
#reference <- "Mlag129A1" 
#reference <- "Mlyc1064a1"
#chromosomes <- read.table("chromosomes.txt")    #args[3]
#synt <- "synteny_ancestral_sp_Mlyc1064a1.txt"  #args[4]
#synt <- "synteny_ancestral_sp_Mlyc1064a2.txt"
#synt <- "synteny_Mlyc1064a1_Mlyc1064a2.txt"
#fai1 <- "haplo1/03_genome/Mlyc1064a1.fa.fai"  #args[5]
#fai1 <- "ancestral_sp/ancestral_sp.fa.fai"     
#fai2 <- "haplo2/03_genome/Mlyc1064a2.fa.fai"
#fai2 <- "haplo1/03_genome/Mlyc1064a1.fa.fai"   #args[6]

#------------- Import other files ---------------------------------------------#
# import synteny data
syn <- read.table(synt, header=T, as.is=T, sep='\t')

# import contig informations from .fai index
index_ref <- read.table(fai1, as.is = T, sep = '\t')[,c(1,2)] %>% 
    set_colnames(., c("chr","end"))

#to fix:
index_hap <- read.table(fai2, as.is = T, sep = '\t')[,c(1,2)] %>% 
    set_colnames(., c("chr","end"))

####------------------------ PREPARE CIRCOS DATA ---------------------------####
#------------- Prepare data sets ----------------------------------------------#
# Check if some contigs have to be inverted
if(ncol(chromosomes)==2) {
  chromosomes$inv=0
}
colnames(chromosomes) <- c("species","chr","inv")

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

#---------- Optional: Prepare table of genes to highlight ---------------------#
#Import the gene positions
if(length(args) == 7) {
data_genes <- read.table("links.txt", as.is=T, sep='\t') #TMP modif
  colnames(data_genes)=c("chr","start","end","category")
  #Invert contig orientation if needed
  for (contig in to_inv)
  {
    end=contigs[which(contigs$chr==contig),]$end
    data_genes[which(data_genes$chr==contig),]$start=end-data_genes[which(data_genes$chr==contig),]$start
    data_genes[which(data_genes$chr==contig),]$end=end-data_genes[which(data_genes$chr==contig),]$end
  }
print(data_genes)

}
####------------------------ LAUNCH CIRCOS ---------------------------------####
#------------- Define plotting parameters -------------------------------------#
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

#------------- Initialize circos ----------------------------------------------#
# Output in pdf
pdf(file = paste0('circos_',haplo,'_on_',reference,'.pdf'))

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

#------------- Plot links -----------------------------------------------------#
#rcols=scales::alpha(ifelse(d$chrom_lag=='MC03',"blue","purple"),alpha=1)
# Color the links according to reference haplotype contigs
# plot links
circos.genomicLink(nuc1, nuc2, col=rcols, border=NA)

#---------- Optional: highlight genes -----------------------------------------#
if(length(args) == 7) { #TMP
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
    for(i in d_cat$chr) {
      circos.genomicRect(d_cat[which(d_cat$chr==i),], sector.index=i,
		track.index=2, ytop = 1, ybottom = 0,col=col,border=col)}
    }
} #TMP

#----------- Write pdf file ---------------------------------------------------#
dev.off()

################################################################################
