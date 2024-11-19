#!/usr/bin/env Rscript

#Purpose:
#script to run rideogram plot 
#input required:  
   # 1 - a list of single copy orthologs for the target species pairs
   # 2 - the two bed files for the target species 
   # 3 - the genomes indexes
   # 4 - optionally: a link file of the type "gene1\tgene2\tstatus" with the gene being the single copy orthologs and a status that will be used for colors
	#to do: add option for position based links 

#--------------- check if library are installed -------------------------------#
if("dplyr" %in% rownames(installed.packages()) == FALSE)
{install.packages("dplyr", repos="https://cloud.r-project.org") }
if("RIdeogram" %in% rownames(installed.packages()) == FALSE)
{install.packages("RIdeogram", repos="https://cloud.r-project.org") }
if("magrittr" %in% rownames(installed.packages()) == FALSE)
{install.packages("magrittr", repos="https://cloud.r-project.org") }
if("data.table" %in% rownames(installed.packages()) == FALSE)
{install.packages("data.table", repos="https://cloud.r-project.org") }

#---------------- load libraries ---------------------------------------------#
libs <- c('dplyr','RIdeogram','magrittr','data.table')
invisible(lapply(libs, suppressMessages(library), character.only = TRUE))

#read species name from the 
argv <- commandArgs(T)

single_copy_ortho <- argv[1]
sco <- read.table(single_copy_ortho, sep = "\t") %>% select(-V1) 

# test if there is at least one argument: if not, return an error
if (length(argv)<=4) {
	  stop("At least the bed and index names of 2 species to compare must be supplied.n", call.=FALSE)
} else if (length(argv)==5) {
	print("assuming no particular link to be highlighted")
	bedA <- argv[2]   #full path to bed file for sp1
	bedB <- argv[3]   #full path to bed file for sp2
	indexA <- argv[4] #full path to index file for sp1
	indexB <- argv[5] #full path to index file for sp2

} else {
	print("links provided in the links file will be displayed in colors")
	print("link file must contain name of gene for species1, name of ortholog for species2 and a status that will be used for coloring the gene")
	#to do: add option to provide a coordinate file with status instead of gene file 

	bedA <- argv[2]   #full path to bed file for sp1
	bedB <- argv[3]   #full path to bed file for sp2
	indexA <- argv[4] #full path to index file for sp1
	indexB <- argv[5] #full path to index file for sp2

       	link <- argv[6] 
	links <- read.table(link, stringsAsFactors = T) %>% set_colnames(.,c("gene1", "gene2","status"))	
	#we will create a vector of color according to the number of status
}

sp1 <- basename(gsub(".bed", "",bedA)) 
sp2 <- basename(gsub(".bed", "",bedB)) 

#read bed files
#they will be use to create the jointed file:
bed1 <- read.table(bedA) %>% 
    merge(sco, ., by.x="V2", by.y = "V4", sort = F ) %>% 
    select(-V3.x) %>%
    set_colnames(., c("gene1", "contig1", "Start_1", "End_1") ) %>%
    mutate(species1 = sp1 ) %>%
    select(species1, gene1, contig1, Start_1, End_1) 

bed2 <- read.table(bedB) %>%
    merge(sco, ., by.x="V3", by.y = "V4", sort = F ) %>% 
    select(-V2.x) %>%
    set_colnames(., c("gene2", "contig2", "Start_2", "End_2") ) %>%
    mutate(species2  = sp2 ) %>%
    select(species2, gene2, contig2, Start_2, End_2) 


n1 <- nrow(bed1)
n2 <- nrow(bed2)
print(paste0("there is ", n1, " single copy gene in bed1"))
print(paste0("there is ", n2, " single copy gene in bed2"))

#to do: add a check and exit with error if n1 != n2 

#-------------- merging bed1 and bed2 - fill colors - create rank to match RIdeogram weird requirement 
#---- rename the rank as species
#---- select wanted columns to match RIdeogram requirements: 
#we will merge the bed1 and bed2 


if (length(argv)==5) {
all <- cbind(bed1, bed2) %>% group_by(contig1) %>% 
    filter(n()>4) %>% group_by(contig2) %>% filter(n()>4) %>%
    mutate(fill = 'cccccc') %>%
    as.data.frame() %>%  mutate(Species_1 = dense_rank(contig1)) %>%
    mutate(Species_2 = dense_rank(contig2)) %>% 
    #select(Species_1,Start_1,End_1,Species_2,Start_2,End_2,fill) %>%
    as.data.frame(.)
} else {
    #assumming links were provided
    all <- cbind(bed1, bed2) %>% group_by(contig1) %>% 
    filter(n()>4) %>% group_by(contig2) %>% filter(n()>4) %>%
    mutate(fill = 'cccccc') %>%
    as.data.frame() %>%  mutate(Species_1 = dense_rank(contig1)) %>%
    mutate(Species_2 = dense_rank(contig2)) %>% 
    #select(Species_1,Start_1,End_1,Species_2,Start_2,End_2,fill) %>%
    as.data.frame(.)

    #assuming we have a link file that is provided
    #some cols:
    colS <- c("f1bb7b", "fd6467","5b1a18","5b1a88","d67236")
    col_pal <- c("#2b8cbe","#de2d26", "#fc9272", "#fee0d2" , "#edf8b1" ,"#636363" )
    links$fill <- rep(colS[1:length(levels(links$status))], c(data.frame(table(links$status))[,2]))
    #finally we use matching of gene to have it all together:
    all$fill[match(links$gene1,all$gene1)] <- links$fill
}

#export the joint bed:
write.table(all, paste0("joint_", sp1, "_" , sp2, ".bed" ) , sep="\t", row.names =F, col.names=T, quote = F)

#here it would be important to check that the order of the genes in one or the two species is identical
#to the order of the genes in the sco! 


#read index to filter chromosome and create a pseudobed file :
#structure of the pseudobed:
#Chr     Start   End      fill    species size    color
#Chr01   0       101369167 969696  Sconica 12      252525

#create pseudo-index1/
index1 <- read.table(indexA)  %>% 
    select(V1,V2) %>%
    filter(V1 %in% all$contig1) %>% 
    mutate(Chr = V1, Start = 0, End = V2, fill = 969696,species = sp1, size = 12, color = 252525) %>%
    select(-V1, -V2)
  

#create pseudo-index2
index2 <- read.table(indexB)  %>% 
    select(V1,V2) %>%
    filter(V1 %in% all$contig2) %>% 
    mutate(Chr = V1, Start = 0,  End = V2, fill = 969696, species = sp2, size = 12, color = 252525) %>%
    select(-V1, -V2)


#combine contig1 and contig 2
karyo <- bind_rows(index1, index2)
#karyo

#----  RIdeogram plot genomes of the same size
#---- this is a probelm if genome size if different. 
#----  so we create a false chromosome that a the size of the difference between the 2 genomes
gapsize <- karyo %>% group_by(species) %>% summarise(size = sum(End) ) %>% summarise(diff = max(size) - min(size))
sp <- karyo %>% group_by(species) %>% summarise(size = sum(End) ) %>%  filter(size == min(size)) %>% select(species)

small  <- data.frame(Chr = "none",  
                     Start = 0, 
                     End = gapsize$diff, 
                     fill = "#0000FF00",  #"#FF0000CC", # "fffffff", #empty fill
                     species = sp$species, #name of 
                     size = 12, 
                     color = 25252525) 
#small
print("index1 is :")
index1 

print("-----------------")
print("index2 is: ")
index2

print("-----------------")
print("karyo : ")
karyo 

print("-----small is: ")
small

if(small$species==index1$species[1]) {
karyo <- rbind(index1, small, index2)
} else { karyo <- rbind(index1, index2, small) }


#/!\/!\
#bits of code to rework depending on wether species 1 or species 2 is the one with the smallest genome size 
#if not ordered properly the script will fail
#karyo <- karyo[c(1,3,2),]
#karyo <- karyo[order(karyo$species),]

all %<>% select(Species_1,Start_1,End_1,Species_2,Start_2,End_2,fill) 

ideogram(karyotype = karyo, synteny = all, output=paste0('plots/', sp1,sp2,'.svg'))
convertSVG(paste0('plots/', sp1,sp2,'.svg', sep=''), file = paste0('plots/', sp1,sp2,'.pdf'), device = "pdf")
