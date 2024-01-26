#!/usr/bin/env Rscript

#Purpose:  Script to plot and filter Ds values based on various metadata ##########
#Author: QR
#Date: 12-05-23

#INPUT: 
# 1 - bed files for the ancestral species and other species
# 2 - paml results files (will be read automatically if previous steps were sucessfull)
# 3 - single copy orthologs  (will be read automatically if previous steps were sucessfull)

argv <- commandArgs(T)

if (argv[1]=="-h" || length(argv)==0){
cat("run as:\n./03.plot_paml_micro.R bed1 bed2 bed3 \n" )
}else{

## read data  from the command line ----------------------------------------#

bed1 <- argv[1] #bed for the ancestral gene order - this one is constant for now
bed2 <- argv[2]  #bed for sp1
bed3 <- argv[3]  #bed for sp2

print(paste0("bed file for sp1 is ", bed2 ))
print(paste0("bed file for sp2 is ", bed3 ))

#--------------- check if library are installed -------------------------------#

if("dplyr" %in% rownames(installed.packages()) == FALSE)
{install.packages("dplyr", repos="https://cloud.r-project.org") }
if("ggplot2" %in% rownames(installed.packages()) == FALSE)
{install.packages("ggplot2", repos="https://cloud.r-project.org") }
if("magrittr" %in% rownames(installed.packages()) == FALSE)
{install.packages("magrittr", repos="https://cloud.r-project.org") }
if("cowplot" %in% rownames(installed.packages()) == FALSE)
{install.packages("cowplot", repos="https://cloud.r-project.org") }
if("wesanderson" %in% rownames(installed.packages()) == FALSE)
{install.packages("wesanderson", repos="https://cloud.r-project.org") }


#---------------- load libraries ---------------------------------------------#
libs <- c('dplyr','ggplot2','magrittr','cowplot','wesanderson')
invisible(lapply(libs, library, character.only = TRUE))

#--------------- load thne data -----------------------------------------------#
## Mv-lag-position:
bed <- read.table(bed1.bed) %>% set_colnames(., c("scaff","start","end","Gene"))


## sp1 + sp2
bedSp1 <- read.table(bed2) %>% set_colnames(., c("scaffSp1","startSp1","endSp1","geneX" ))
bedSp2 <- read.table(bed3) %>% set_colnames(., c("scaffSp1","startSp2","endSp2","geneY"))


## yn00 results:
dat <- read.table("paml/results_YN.txt") %>% set_colnames(., c("Ds","SEDs","Dn","SEDn", "geneX", "geneY"))

## single copy orthologs:
single_cp <- read.table("paml/single.copy.orthologs") %>% set_colnames(., c("ortho","Gene","geneX","geneY" ))

## --------------------- generic function ------------------------------------------------- ##

`%nin%` = Negate(`%in%`) #to negate 


## ------------- arrange the data as needed ----------------------------------------------- ##
Ds_table <- merge(dat, single_cp, by.x = "geneX", by.y = "geneX")
Ds_table <- merge(Ds_table, bedLag,    by.x = "Gene", by.y = "Gene") 

Ds_table$status <- ifelse(Ds_table$scaffLag == "Mlag129A1_contig_11", "PR", "HD")


D <- na.omit(Ds_table) #there should not be any!

#this bits is too specific and must be generalized :

#I like to insert HD/PR here:
#D$chr <- paste0(D$scaffLag," (" , D$status, ")" )
#reorder to have PR first!
#D$chr <- factor(D$chr, levels=c("Mlag129A1_contig_8 (HD)", "Mlag129A1_contig_11 (PR)") )
#D <- with(D, D[order(chr),])
#HD <- D %>% 
#    filter(status =="HD") %>%
#    mutate(StartNew = startLag ) %>%
#    arrange(desc(StartNew)) #%>%
#PR <- D %>% 
#    filter(status =="PR") %>%
#    mutate(StartNew = (startLag) ) #%>%

#to rework as well:
all <- rbind(HD, PR) %>%
    mutate(ordre = seq(1:nrow(.))) %>%
    mutate(name = ifelse(mlagGene == "Mlag129A1_contig_8_g7514.t1", "HD1",
    ifelse(mlagGene == "Mlag129A1_contig_8_g7515.t1", "HD2",
    ifelse(mlagGene == "Mlag129A1_contig_11_MP008227.t1", "PR", NA))))

df <- all %>% filter(Ds < 0.3) %>% select(ordre, Ds)

#export the df for model comparison on the cluster:
write.table(df, "dS.modelcomp", quote =F, row.names = F, col.names = T, sep = "\t")


## ------------------ GGPLOT  CUSTOMISATION ------------------------------------------------##
th_plot <-     theme(axis.title.x=element_text(size=14, family="Helvetica",face="bold"),
  axis.text.x=element_text(size=14,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
  axis.title.y=element_text(size=18, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
  axis.text.y=element_text(size=14,family="Helvetica",face="bold"),
  strip.text.x = element_text(size=18),
  panel.grid.major = element_blank())


########################## make plot now using ALL GENES #######################################

mycolor2 <-c("#E69F00",  "#0072B2" ,"#5B1A79",  "#CC79A7", "#D55E00")

Fig1A <- all  %>%   #we plot the D dataframe to obtain the Ds along the order
  filter(Ds < 1.01) %>%
  ggplot(., aes(x = startLag, y = Ds)) +
  geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
  facet_wrap(~chr, scale="free_x") +
  #facet_wrap(~scaffLag, scale="free_x") +
  geom_point( size = 1) + 
  #geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  geom_label_repel(aes(fill = factor(status)), colour = "white", segment.colour = "black",
  #geom_text_repel(
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 0.06,
    direction    = "x",
    angle        = 90,
    hjust        = 0,
    segment.size = 0.4,
    max.iter = 1e4, max.time = 1) +
    #color the gene:
    scale_fill_discrete(type = mycolor2[1:3]) +

  ylim(c(0,1)) +
  xlab("position along chr HD & PR") +
  ylab( expression(italic("Ds"))) +
  th_plot + theme(legend.position = "none") 


#to do: add a trim Y-axis to display PR

#Fig1A
#all_pos
Fig1B <- all %>%   #we plot the D dataframe to obtain the Ds along the order
  filter(Ds < 1) %>%
  ggplot(., aes(x = ordre, y = Ds, colour = chr)) +
  geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
  geom_point( size = 1) + 
  #geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  ylim(c(0,0.6)) +
  xlab("order along M.v. lagerheimii") +
  ylab( expression(italic("Ds"))) +
  th_plot + theme(legend.position = "none") +
  scale_color_manual(values=wes_palette(n=2, name="GrandBudapest1"))   
  
#Fig1B

pdf(file = "plot/Ds.plot.pdf",14,8)
plot_grid(Fig1A, Fig1B, labels="AUTO", ncol = 1)
dev.off()


quit(save = "no") 
}
