

library(GENESPACE)

load("results/gsParams.rda")
roi <- data.frame(genome=c("Mlag129A1"), chr=c("Mlag129A1_contig_11","Mlag129A1_contig_8"))


pdf(file = "HD_PR_geneOrder_useRegionF.pdf",8,6)
ripd <- plot_riparian(gsParam=gsParam, highlightBed=roi, backgroundColor=NULL, useRegions = F, minChrLen2plot = 0)
dev.off()

pdf(file = "HD_PR_bp_useRegionF.pdf",8,6)
ripd <- plot_riparian(gsParam=gsParam, highlightBed=roi, backgroundColor=NULL, useRegions = FALSE, useOrder = FALSE, minChrLen2plot = 0)
dev.off()

invchr <- roi

pdf(file = "HD_PR_bp_useRegionF.pdf",8,6)
ripd <- plot_riparian(gsParam=gsParam, highlightBed=roi,
		      invertTheseChrs = invchr,
		      backgroundColor=NULL, useRegions = FALSE, useOrder = FALSE, minChrLen2plot = 0)
dev.off()
