

library(GENESPACE)
gpar <- init_genespace(
  wd = "./",
  path2mcscanx="/home/quentin/.software/MCScanX/",
 path2orthofinder = "/home/quentin/.software/OrthoFinder/orthofinder")
gpar <- run_genespace(gsParam = gpar) 

