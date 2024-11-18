#!/usr/bin/env Rscript

#Purpose:  Script to plot and filter Ds values based on various metadata ##########
#Author: QR
#Date: 12-05-23

#INPUTs -- 3 arguments needed:  
# 1 - bed files for the species1
# 2 - bed files for the species2
# 3 - a string: "N" or "R" for "Normal" or Reversed (R): in the case were the focal region is spread on two scaffold, this string should state wether the second scaffold should be reversed or not. 
#this will not work for more than two scaffold

# optional:
# 4 - bed file for the ancestral species 

# - paml results files will be read automatically if previous steps were sucessfull

argv <- commandArgs(T)

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

## --------------------- generic function ------------------------------------------------- ##

`%nin%` = Negate(`%in%`) #to negate 

#--------------- load thne data -----------------------------------------------#

#- common results
# yn00 results:
dat <- read.table("paml/results_YN.txt") %>% 
	set_colnames(., c("Ds","SEDs","Dn","SEDn", "geneX", "geneY"))



if (length(argv)<3) {
	  stop("At least the name of 2 species to compare and a txt file containing the name and order of scaffold must be supplied.n", call.=FALSE)
} else if (length(argv)==3) {
	print("assuming no ancestral species was used")
	sp1 <- argv[1]     # only the basename is needed !
	sp2 <- argv[2]     # only the basename is needed !
	chr <- argv[3]     # table with chr\tstatus [Reversed or Not]

	scaf <- read.table(chr, sep="\t") %>% set_colnames(., c("haplo","chr","order"))

	#orthofinder single copy orthologs:
	single_cp <- read.table("paml/single.copy.orthologs", sep = "\t") %>% 
		     set_colnames(., c("ortho","geneX","geneY" ))

} else {
	print("assuming an ancestral species exist")
	sp1 <- argv[1]     #only the basename is needed !
	sp2 <- argv[2]     #only the basename is needed !
	chr <- argv[3]     # table with chr\tstatus [Reversed or Not]
	#optional 
	sp3 <- argv[4]     #the basename of the ancestral species !
	
	print("load scaffold info")
	scaf <- read.table(chr, sep ="\t") %>% set_colnames(., c("haplo","chr","order"))

	#orthofinder single copy orthologs:
	print("load single copy info")
	single_cp <- read.table("paml/single.copy.orthologs", sep = "\t") %>% 
		     set_colnames(., c("ortho","gene","geneX","geneY" ))

	#link <- argv[6] 
	#links <- read.table(link, stringsAsFactors = T) %>% set_colnames(.,c("gene1", "gene2","status"))	
	#we will create a vector of color according to the number of status
	
	## read Ancestral species :
	print("load ancestral species info")
	bedAnc <- read.table(paste0("genespace/bed/",sp3, ".bed", sep = "" )) %>% 
		set_colnames(., c("scaff","start","end","gene"))

}

## sp1 + sp2
bedSp1 <- read.table(paste0("genespace/bed/",sp1, ".bed", sep = "" )) %>% 
	set_colnames(., c("scaff","start","end","gene" ))
bedSp2 <- read.table(paste0("genespace/bed/",sp2, ".bed", sep = "" )) %>% 
	set_colnames(., c("scaffSp2","startSp2","endSp2","geneY"))

## ------------- arrange the data as needed ----------------------------------------------- ##
Ds_table <- merge(dat, single_cp, by.x = "geneX", by.y = "geneX")

#merge with the coordinate of the reference sequence (either Sp1 or ancestral species):
if (length(argv)==3) {
    Ds_table <- merge(Ds_table, bedSp1,    by.x = "geneX", by.y = "gene") 
} else {
    Ds_table <- merge(Ds_table, bedAnc,    by.x = "gene", by.y = "gene") 

}

#now we must: 
    #1 - reorder according to the scaffold orientation
    #2 - create an incremenantial gene order accordingly:
if (exists(sp3) {
    #assuming ancestral species was provided
    all <- merge(bedAnc, scaf, by.x = "scaff", by.y = "chr") %>%
        left_join(., Ds_table, by=join_by(gene == gene) ) %>%
        arrange(scaff, start) %>%
        group_by(scaff) %>%
        mutate(St = ifelse(order == "N", start, rev(start) )) %>% 
        arrange(St, .by_group = TRUE) %>%
        ungroup() %>%
        mutate(orderchp = seq(1:nrow(.)))
}
else {
    #assuming non ancestral species 
    #plotting along the X:
    all <- merge(bedSp1, scaf, by.x = "scaff", by.y = "chr") %>%
        left_join(., Ds_table, by=join_by(gene == gene) ) %>%
        arrange(scaff, start) %>%
        group_by(scaff) %>%
        mutate(St = ifelse(order == "N", start, rev(start) )) %>% 
        arrange(St, .by_group = TRUE) %>%
        ungroup() %>%
        mutate(orderchp = seq(1:nrow(.)))
}

#Ds values above 0.3 will be considered as pseudo-genes for the changepoint analyses. 
#df <- all %>% filter(Ds < 0.3) %>% select(order, Ds)
df <- all %>% 
    filter(Ds < 0.3) %>% 
    select(scaff, start, order, orderchp, St, Ds) %>%
    na.omit()

#export the df for model comparison on the cluster:
write.table(df, "dS.values.forchanepoint.txt", quote =F, row.names = F, col.names = T, sep = "\t")

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
  ggplot(., aes(x = start, y = Ds )) +
  geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
  facet_wrap(~scaff, scale="free_x") +
  geom_point( size = 1) + 
  #geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  #geom_label_repel(aes(fill = factor(status)), colour = "white", segment.colour = "black",
  #geom_text_repel(
  #  force_pull   = 0, # do not pull toward data points
  #  nudge_y      = 0.06,
  #  direction    = "x",
  #  angle        = 90,
  #  hjust        = 0,
  #  segment.size = 0.4,
  #  max.iter = 1e4, max.time = 1) +
    #color the gene:
  #  scale_fill_discrete(type = mycolor2[1:3]) +

  ylim(c(0,1)) +
  xlab("position along chr") +
  ylab( expression(italic("Ds"))) +
  th_plot + theme(legend.position = "none") 

#to do: add a trim Y-axis to display high Ds genes

Fig1B <- all %>%   #we plot the D dataframe to obtain the Ds along the order
  filter(Ds < 1) %>%
  ggplot(., aes(x = order, y = Ds, colour = scaff)) +
  geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
  geom_point( size = 1) + 
  theme_classic() +
  ylim(c(0,0.6)) +
  xlab("order along reference") +
  ylab( expression(italic("Ds"))) +
  th_plot + theme(legend.position = "none") +
  scale_color_manual(values=wes_palette(n=2, name="GrandBudapest1"))   
  
#Fig1B

pdf(file = "plots/Ds.pdf",14,8)
plot_grid(Fig1A, Fig1B, labels="AUTO", ncol = 1)
dev.off()


print("-------------------------------------------------------")
print("------- constructing graph with gene order-------------")

#only if we have an ancestral reference, otherwise it is a bit meaningless

if(length(argv)==4){
	# --- now add the order
	colnames(bedSp1) <- c("scaffSp1","startSp1","endSp1","geneX") 
	ordSp1<- merge(all, bedSp1, by.x = "geneX", by.y = "geneX", sort = F) %>%
  		group_by(scaff) %>%
  		filter(n()>2) %>%
	  ungroup() %>%
	  mutate(rankA1 = dense_rank(startSp1))
	
	ordSp2 <-  merge(all, bedSp2, by.x = "geneY.y", by.y = "geneY")%>%
	  group_by(scaff) %>%
	  filter(n()>2) %>%
	  ungroup() %>%
	  mutate(rankA1 = dense_rank(startSp2))
	
	
	pordSp1 <- ggplot(ordSp1, aes(x = order, y = rankA1, colour = scaffSp1 )) +
	  geom_point( size = 2) + 
	  #geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
	  theme_classic() +
	  #ylim(c(0,0.5)) +
	  xlab("order along reference") +
	  ylab( expression(italic("gene rank in A1"))) +
	  th_plot + theme(legend.position = "none") +
	  theme(axis.text.y=element_blank(),
	        axis.ticks.y=element_blank() 
	  ) +
	  scale_color_viridis(discrete=TRUE) 
	  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest1"))   
	
	pordSp2 <- ggplot(ordSp2, aes(x = order, y = rankA1, colour = scaffSp2 )) +
	  geom_point( size = 2) + 
	  #geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
	  theme_classic() +
	  #ylim(c(0,0.5)) +
	  xlab("order along reference") +
	  ylab( expression(italic("gene rank in A2"))) +
	    th_plot + theme(legend.position = "none") +
	  theme(axis.text.y=element_blank(),
	        axis.ticks.y=element_blank() 
	  ) +
	  scale_color_viridis(discrete=TRUE) 
	  #scale_color_manual(values=wes_palette(n=4, name="GrandBudapest1"))   
	
	pdf(file = "plots/Ds_and_arrangements.pdf",18,20)
	print(plot_grid(Fig1A, Fig1B, pordSp1, pordSp2, labels="AUTO", ncol = 1, rel_heights = c(1,1,0.9,0.9)) )
	dev.off()
}

