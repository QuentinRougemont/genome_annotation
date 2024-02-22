
argv <- commandArgs(T)
genome_id <- argv[1]

library(GENESPACE)

load("results/gsParams.rda")

scaffold <- read.table("scaffold.txt")
roi <- data.frame(scaffold[,-3])
colnames(roi) <-c("genome","chr")
roi$genome <- genome_id



pdf(file = "subset_geneOrder_useRegionF.pdf",8,6)
ripd <- plot_riparian(gsParam=gsParam, highlightBed=roi, backgroundColor=NULL, useRegions = F, minChrLen2plot = 0)
dev.off()

pdf(file = "subset_bp_useRegionF_useOrderF.pdf",8,6)
ripd <- plot_riparian(gsParam=gsParam, highlightBed=roi, backgroundColor=NULL, useRegions = FALSE, useOrder = FALSE, minChrLen2plot = 0)
dev.off()

invchr <- roi

pdf(file = "subset_useRegionF_inverted.pdf",8,6)
ripd <- plot_riparian(gsParam=gsParam, highlightBed=roi,
     invertTheseChrs = invchr,
     backgroundColor=NULL, useRegions = FALSE, useOrder = FALSE, minChrLen2plot = 0)
dev.off()
