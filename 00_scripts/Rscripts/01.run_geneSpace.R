#!/usr/bin/env Rscript

#nano script to simply run genespace:
library(GENESPACE)

gpar <- init_genespace(path2mcscanx="mcpath",  wd = "./" )

gpar <- run_genespace(gsParam = gpar) 

