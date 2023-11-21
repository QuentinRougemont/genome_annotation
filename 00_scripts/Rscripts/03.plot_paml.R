#!/usr/bin/env Rscript

#Purpose:  Script to plot and filter Ds values based on various metadata ##########
#Author: QR
#Date: 12-05-23

#INPUT: 
# 1 - bed files for the ancestral species and other species
# 2 - paml results files (will be read automatically if previous steps were sucessfull)
# 3 - single copy orthologs  (will be read automatically if previous steps were sucessfull)

argv <- commandArgs(T)

##if (argv[1]=="-h" || length(argv)==0){
#cat("run as:\n./03.plot_paml_micro.R bed1 bed2 bed3 \n" )
#}else{

## read data  from the command line ----------------------------------------#

sp1 <- argv[1]  #bed for the ancestral gene order - this one is constant for now
sp2 <- argv[2]  #bed for sp1
sp3 <- argv[3]  #bed for sp2

#sp1 <- as.character(argv[1])  #ancestral species name
#sp2 <- argv[2]  #species1 
#sp3 <- argv[3]  #species2

print(paste0("bed file for sp1 is ", sp1 ))
print(paste0("bed file for sp2 is ", sp2 ))

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
if("viridis" %in% rownames(installed.packages()) == FALSE)
{install.packages("viridis", repos="https://cloud.r-project.org") }
if("ggrepel" %in% rownames(installed.packages()) == FALSE)
{install.packages("ggrepel", repos="https://cloud.r-project.org") }


#---------------- load libraries ---------------------------------------------#
libs <- c('dplyr','ggplot2','magrittr','cowplot','wesanderson', 'viridis','ggrepel')
invisible(lapply(libs, library, character.only = TRUE))

#--------------- load thne data -----------------------------------------------#
## Mv-lag-position:
#bedLag <- read.table("bed/Mlag129A1.bed") %>% set_colnames(., c("scaffLag","startLag","endLag","mlagGene"))
bedLag <- read.table(paste0("genespace/bed/",sp1, ".bed", sep = "" )) %>% set_colnames(., c("scaffLag","startLag","endLag","mlagGene"))
## sp1 + sp2
bedSp1 <- read.table(paste0("genespace/bed/",sp2, ".bed", sep = "" )) %>% set_colnames(., c("scaffSp1","startSp1","endSp1","geneX" ))
bedSp2 <- read.table(paste0("genespace/bed/",sp3, ".bed", sep = "" )) %>% set_colnames(., c("scaffSp1","startSp2","endSp2","geneY"))


## yn00 results:
dat <- read.table("paml/results_YN.txt") %>% set_colnames(., c("Ds","SEDs","Dn","SEDn", "geneX", "geneY"))

## single copy orthologs:
single_cp <- read.table("paml/single.copy.orthologs") %>% set_colnames(., c("ortho","mlagGene","geneX","geneY" ))

## --------------------- generic function ------------------------------------------------- ##

`%nin%` = Negate(`%in%`) #to negate 


## ------------- arrange the data as needed ----------------------------------------------- ##
Ds_table <- merge(dat, single_cp, by.x = "geneX", by.y = "geneX")
Ds_table <- merge(Ds_table, bedLag,    by.x = "mlagGene", by.y = "mlagGene") 

Ds_table$status <- ifelse(Ds_table$scaffLag == "Mlag129A1_contig_11", "PR", "HD")


D <- na.omit(Ds_table) #there should not be any!
#I like to insert HD/PR here:
D$chr <- paste0(D$scaffLag," (" , D$status, ")" )

#reorder to have PR first!
D$chr <- factor(D$chr, levels=c("Mlag129A1_contig_8 (HD)", "Mlag129A1_contig_11 (PR)") )
D <- with(D, D[order(chr),])

HD <- D %>% 
    filter(status =="HD") %>%
    mutate(StartNew = startLag ) %>%
    arrange(desc(StartNew)) #%>%

PR <- D %>% 
    filter(status =="PR") %>%
    mutate(StartNew = (startLag) ) #%>%


all <- rbind(HD, PR) %>%
    mutate(ordre = seq(1:nrow(.))) %>%
    mutate(name = ifelse(mlagGene == "Mlag129A1_contig_8_g7514.t1", "HD1",
    ifelse(mlagGene == "Mlag129A1_contig_8_g7515.t1", "HD2",
    ifelse(mlagGene == "Mlag129A1_contig_11_MP008227.t1", "PR", NA))))

df <- all %>% filter(Ds < 0.3) %>% select(ordre, Ds)

#export the df for model comparison on the cluster:
write.table(df, "dS.modelcomp", quote =F, row.names = F, col.names = T, sep = "\t")

#D <- D %>% 
#    mutate(startNew = ifelse(status =="HD", rev(startLag), startLag )) %>%
#    mutate(ordre = seq(1:nrow(.)))


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
  ggplot(., aes(x = startLag, y = Ds, label = name)) +
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

pdf(file = "plots/Ds.pdf",14,8)
plot_grid(Fig1A, Fig1B, labels="AUTO", ncol = 1)
dev.off()



print("-------------------------------------------------------")
print("------- constructing graph with gene order-------------")

# --- now add the order
ordSp1<- merge(all, bedSp1, by.x = "geneX", by.y = "geneX", sort = F) %>%
  group_by(scaffSp1) %>%
  filter(n()>2) %>%
  ungroup() %>%
  mutate(rankA1 = dense_rank(startSp1))

print("-------------------------------------------------------")
head(bedSp2)
print("-------------------------------------------------------")
head(all)

ordSp2 <-  merge(all, bedSp2, by.x = "geneY.y", by.y = "geneY")%>%
  group_by(scaffSp1) %>%
  filter(n()>2) %>%
  ungroup() %>%
  mutate(rankA1 = dense_rank(startSp2))


pordSp1 <- ggplot(ordSp1, aes(x = ordre, y = rankA1, colour = scaffSp1 )) +
  geom_point( size = 2) + 
  #geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  #ylim(c(0,0.5)) +
  xlab("order along M.v. lagerheimii") +
  ylab( expression(italic("gene rank in A1"))) +
  th_plot + theme(legend.position = "none") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) +
  scale_color_viridis(discrete=TRUE) 
  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest1"))   

pordSp2 <- ggplot(ordSp2, aes(x = ordre, y = rankA1, colour = scaffSp1 )) +
  geom_point( size = 2) + 
  #geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  #ylim(c(0,0.5)) +
  xlab("order along M.v. lagerheimii") +
  ylab( expression(italic("gene rank in A2"))) +
    th_plot + theme(legend.position = "none") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() 
  ) +
  scale_color_viridis(discrete=TRUE) 
  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest1"))   

pdf(file = "plots/Ds_and_arrangements.pdf",18,20)
plot_grid(Fig1A, Fig1B, pordSp1, pordSp2, labels="AUTO", ncol = 1, rel_heights = c(1,1,0.9,0.9))
dev.off()

#quit(save = "no") 
#}

