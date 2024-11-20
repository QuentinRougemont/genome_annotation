#!/usr/bin/env Rscript 

if("pafr" %in% rownames(installed.packages()) == FALSE)
{install.packages("pafr", repos="https://cloud.r-project.org") }
if("ggplot2" %in% rownames(installed.packages()) == FALSE)
{install.packages("ggplot2", repos="https://cloud.r-project.org") }
if("ggpubr" %in% rownames(installed.packages()) == FALSE)
{install.packages("ggpubr", repos="https://cloud.r-project.org") }
if("knitr" %in% rownames(installed.packages()) == FALSE)
{install.packages("knitr", repos="https://cloud.r-project.org") }

library(ggplot2)
library(pafr)
library(ggpubr)
library(knitr)

argv <- commandArgs(T)
paf <- argv[1] #paf file
ali <- read_paf(paf)
base <- basename(paf)

#create dir if not present:
if (!dir.exists("02_results/genomeplots")){
  dir.create("02_results/genomeplots")
}

pdf(file=paste0("02_results/genomeplots/divergence", base, ".pdf") )
ggplot(ali, aes(alen, de)) + 
    geom_point(alpha=0.6, colour="steelblue", size=2) + 
    scale_x_continuous("Alignment length (kb)", label =  function(x) x/ 1e3) +
    scale_y_continuous("Per base divergence") + 
    theme_pubr()
dev.off()

by_q <- aggregate(de ~ qname, data=ali, FUN=mean)
knitr::kable(by_q)

prim_alignment <- filter_secondary_alignments(ali)
nrow(ali) - nrow(prim_alignment)

long_ali <- subset(ali, alen > 1e4 & mapq > 40)
long_ali

pdf(file=paste0("02_results/genomeplots/dotplot_prim_ali",base,".pdf") )
dotplot(prim_alignment)
dev.off()

pdf(file=paste0("02_results/genomeplots/dotplot_long_ali",base,".pdf"))
dotplot(long_ali)
dev.off()

pdf(file=paste0("02_results/genomeplots/dotplot_prim_ali_bw", base, ".pdf"))
dotplot(prim_alignment, label_seqs=TRUE, order_by="qstart") + theme_bw()
dev.off()
