#!/usr/bin/env Rscript

#very simple script to plot synteny between target chromosomes using pafR

if("dplyr" %in% rownames(installed.packages()) == FALSE)
{install.packages("dplyr", repos="https://cloud.r-project.org") }
if("ggplot2" %in% rownames(installed.packages()) == FALSE)
{install.packages("ggplot2", repos="https://cloud.r-project.org") }
if("pafr" %in% rownames(installed.packages()) == FALSE)
{install.packages("pafr", repos="https://cloud.r-project.org") }

library(ggplot2)
library(pafr)
library(dplyr)

argv <- commandArgs(T)
paf <- argv[1]     	#paf alignment
keep0 <- read.table(argv[2]) #table of subscaffold  
base <- basename(paf)

#read paf alignment
ali <- read_paf(paf)

#create dir if not present:
if (!dir.exists("02_results/genomeplots")){
  dir.create("02_results/genomeplots")
}


#filter 
prim_alignment <- filter_secondary_alignments(ali)

keep <- keep0 %>% filter(V1> 30) %>% #filter those with less than 30 links
	select(V2,V3) %>% unique()
to_keep <- list(keep$V3, keep$V2)

pdf(file=paste0("02_results/genomeplots/target_scaff_dotplot_",base, ".pdf"))
dotplot(prim_alignment, label_seqs=TRUE, order_by="provided", ordering=to_keep)
dev.off()


### nrow(keep)
for(i in 1:nrow(keep) ){
    pdf(file = paste0("02_results/genomeplots/rcT_synteny_",base,"_", keep$V2[i],"_", keep$V3[i],".pdf"), 10,6 )
    print(plot_synteny(prim_alignment, q_chrom = keep$V3[i], t_chrom =  keep$V2[i], centre = T, rc= T))
    dev.off()
}

for(i in 1:nrow(keep) ){
    pdf(file = paste0("02_results/genomeplots/synteny_",base,"_", keep$V2[i],"_", keep$V3[i],".pdf") , 10, 6)
    print(plot_synteny(prim_alignment, q_chrom = keep$V3[i], t_chrom =  keep$V2[i], centre = T, rc= F))
    dev.off()
}

