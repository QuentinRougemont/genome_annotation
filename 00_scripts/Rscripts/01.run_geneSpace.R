

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jtlovell/GENESPACE")

library(GENESPACE)

gpar <- init_genespace(
  wd = "./",
  path2mcscanx="/home/quentin/.software/MCScanX/",
 path2orthofinder = "/home/quentin/.software/OrthoFinder/orthofinder")
gpar <- run_genespace(gsParam = gpar) 

