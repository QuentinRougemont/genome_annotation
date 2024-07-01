

if("BiocManager" %in% rownames(installed.packages()) == FALSE)            
{install.packages("BiocManager", repos="https://cloud.r-project.org") }   

BiocManager::install("Biostrings")

if("devtools" %in% rownames(installed.packages()) == FALSE)
{install.packages("devtools", repos="https://cloud.r-project.org") }

devtools::install_github("jtlovell/GENESPACE")

library(GENESPACE)

gpar <- init_genespace(path2mcscanx="mcpath", wd = "./" )

gpar <- run_genespace(gsParam = gpar) 

