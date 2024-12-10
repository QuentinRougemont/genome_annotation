#!/usr/bin/env Rscript

argv <- commandArgs(T)
genome_id <- argv[1]

library(GENESPACE)

load("results/gsParams.rda")

scaffold <- read.table("scaffold.txt")
roi <- data.frame(scaffold[,-3])
colnames(roi) <-c("genome","chr")
roi$genome <- genome_id


pdf(file = "subset_geneOrder_useRegionF.pdf",10,6)
ripd <- plot_riparian(gsParam=gsParam, highlightBed=roi, backgroundColor=NULL, useRegions = F, minChrLen2plot = 0)
dev.off()

pdf(file = "subset_bp_useRegionF_useOrderF.pdf",10,6)
ripd <- plot_riparian(gsParam=gsParam, highlightBed=roi, backgroundColor=NULL, useRegions = FALSE, useOrder = FALSE, minChrLen2plot = 0)
dev.off()

invchr <- roi

pdf(file = "subset_useRegionF_inverted.pdf",10,6)
ripd <- plot_riparian(gsParam=gsParam, highlightBed=roi,
     invertTheseChrs = invchr,
     backgroundColor=NULL, useRegions = FALSE, useOrder = FALSE, minChrLen2plot = 0)
dev.off()


ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  c("darkorange", "skyblue", "darkblue", "purple", "darkred", "salmon"))

#same but using different backrgound :
pdf(file = "white_subset_geneOrder_useRegionF.pdf",12,6)
ripd <- plot_riparian(gsParam=gsParam, 
                      highlightBed=roi, 
                      backgroundColor=NULL, 
                      useRegions = F, 
                      minChrLen2plot = 0,
                      palette = customPal,
                      braidAlpha = .75, 
                      addTheme=ggthemes,
                      chrFill = "lightgrey")
dev.off()

pdf(file = "white_subset_geneOrder_useRegionF_inverted.pdf",12,6)
ripd <- plot_riparian(gsParam=gsParam, 
                      highlightBed=roi, 
                      backgroundColor=NULL, 
                      useRegions = F, 
                      minChrLen2plot = 0,
                      palette = customPal,
                      braidAlpha = .75, 
                      addTheme=ggthemes,
                      invertTheseChrs = invchr,
                      chrFill = "lightgrey")
dev.off()

#add option for more inverted chromosome possibilities ??
