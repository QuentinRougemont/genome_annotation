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

#bed1 <- argv[1] #bed for the ancestral gene order - this one is constant for now
bed2 <- argv[1]  #bed for sp1
bed3 <- argv[2]  #bed for sp2

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
bedLag <- read.table("bed/Mlag129A1.bed") %>% set_colnames(., c("scaffLag","startLag","endLag","mlagGene"))


## sp1 + sp2
bedSp1 <- read.table(bed2) %>% set_colnames(., c("scaffSp1","startSp1","endSp1","geneX" ))
bedSp2 <- read.table(bed3) %>% set_colnames(., c("scaffSp1","startSp2","endSp2","geneY"))


## yn00 results:
dat <- read.table("results_YN.txt") %>% set_colnames(., c("Ds","SEDs","Dn","SEDn", "geneX", "geneY"))

## single copy orthologs:
single_cp <- read.table("single.copy.orthologs") %>% set_colnames(., c("ortho","mlagGene","geneX","geneY" ))

## --------------------- generic function ------------------------------------------------- ##

`%nin%` = Negate(`%in%`) #to negate 


## ------------- arrange the data as needed ----------------------------------------------- ##
Ds_table <- merge(dat, single_cp, by.x = "geneX", by.y = "geneX")
Ds_table <- merge(Ds_table, bedLag,    by.x = "mlagGene", by.y = "mlagGene") 

Ds_table$status <- ifelse(Ds_table$scaffLag == "Mlag129A1_contig_11", "PR", "HD")


D <- na.omit(Ds_table) #there should not be any!
head(D)
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

pdf(file = "Ds.plot.pdf",14,8)
plot_grid(Fig1A, Fig1B, labels="AUTO", ncol = 1)
dev.off()


quit(save = "no") 
}

#------------------------------------------------------------------------------#
#                   perform the changepoint analyis here: 
#------------------------------------------------------------------------------#

#define the model we want to test:
model4st = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1)
model3strata = list(Ds ~ 1, 1~ 1, 1 ~ 1)
model2strata = list(Ds ~ 1, 1~ 1)
modelsimple = list(Ds ~ 1 + ordre)
model5st = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1, 1 ~ 1)

df <- D %>% filter(Ds < 0.3) %>% select(ordre, Ds)

#export the df for model comparison on the cluster:
write.table(df, "dS.modelcomp", quote =F, row.names = F, col.names = T, sep = "\t")

fit_2st = mcp(model2strata,  data = df, par_x = "ordre",   iter = 1e5, adapt = 15e3,  chains = 5, cores = 5 )
fit_3st = mcp(model3strata,  data = df, par_x = "ordre",   iter = 1e5, adapt = 15e3,  chains = 5, cores = 5 )
fit_1st = mcp(modelsimple,  data = df, par_x = "ordre",   iter = 1e5, adapt = 15e3,  chains = 5, cores = 5 )
fit_4st = mcp(model4st,  data = df, par_x = "ordre",   iter = 1e5, adapt = 15e3,  chains = 5, cores = 5 )
fit_5st = mcp(model5st,  data = df, par_x = "ordre",   iter = 1e5, adapt = 15e3,  chains = 5, cores = 5 )

#test which model is the best here: 

#summarise and look at point change:
summary(fit_4st)
summary(fit_3st)
summary(fit_2st)

#test hypothesis relative to the change;
hypothesis(fit_4st, c("int_1 < int_2", "int_2 < int_3", "int_3 > int_4"))
hypothesis(fit_3st, c("int_1 < int_2", "int_2 > int_3" ))

#store it for later
xl <- expression(paste("order along ", italic("M. lagerheimeii"), " mating chromsome"))
FigMCP4Strata  <- plot(fit_4st, q_fit = TRUE) + ggplot2::ggtitle("Posterior fit - 4 strata") + theme_classic() + th_plot +
  geom_point(color = "darkblue", size = 0.1) +
  xlab(xl) + ylab(expression(italic(d[s])))

FigMCP3Strata  <- plot(fit_3st, q_fit = TRUE) + ggplot2::ggtitle("Posterior fit - 3 strata") + theme_classic() + th_plot +
  geom_point(color = "darkblue", size = 0.21) +
  xlab(xl) + ylab(expression(italic(d[s])))

FigMCP2Strata  <- plot(fit_2st, q_fit = TRUE) + ggplot2::ggtitle("Posterior fit - 2 strata") + theme_classic() + th_plot +
  geom_point(color = "darkblue",  size = 0.21) +
  xlab(xl) + ylab(expression(italic(d[s])))

pdf(file = "original_posterior_4strata_401Ds.pdf",12,6)
FigMCP4Strata
dev.off()
pdf(file = "original_posterior_3strata_401Ds.pdf",12,6)
FigMCP3Strata
dev.off()
pdf(file = "original_posterior_2strata_401Ds.pdf",12,6)
FigMCP2Strata
dev.off()

#plot themall
pdf(file = "Strata_compMnot.pdf", 12,14)
plot_grid(FigMCP2Strata, FigMCP3Strata, FigMCP4Strata, labels = "AUTO", ncol = 1)
dev.off()

#statistical test of the difference among strata:

#extraite les cp des estimations de changepoints: 
location_cp1 <-
location_cp2 <-
location_cp3 <-
  

dat3 <- dat3 %>% 
  mutate(strata = ifelse(start < location_cp1, "strata1",ifelse(start > location_cp2, "strata3", "strata2")) )





################################################################################
#Final figure for Mnot: 
#Fig1A: Ds along genome
#Fig1B: Ds along order
#Fig1C: changepoint
#Fig1D: rang in a1 and a2

#Fig2: statisical comparison of Ds and of Dn/Ds



#Fig3: 


#Fig0: geneSpace M-lag and target species only

#Fig0b: minimap whole genome D/T, minimap LagA1 et D, et zoom sur chr

#
################################################################################






awk '(length($7)>=3 || length($8)>=3) && ((gsub(/a/,"&",$7)>=length($7)-1|| gsub(/c/,"&",$7)>=length($7)-1 || gsub(/t/,"&",$7)>=length($7)-1 || gsub(/g/,"&",$7)>=length($7)-1) || (gsub(/a/,"&",$8)>=length($8)-1|| gsub(/c/,"&",$8)>=length($8)-1 || gsub(/t/,"&",$8)>=length($8)-1 || gsub(/g/,"&",$8)>=length($8)-1))  ' out.var |grep 'MC03\|MC12\|MC16' |awk '(length($7)<50 && length($8)<50) {print $2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > potential.homopolymers.txt 

homo <- read.table("potential.homopolymers.txt") %>% select(V4,V5,V6,V7,V8) %>% 
  set_colnames(. , c("h1","h2","scaff","start","end")) %>% data.table(.)

setkey(homo, scaff, start, end)
homo$homo <- "homopolymers"
Ds_table <- data.table(Ds_table)
setkey(Ds_table, scaff, start, end)

t2 <- foverlaps(Ds_table,homo, type ="any")
t2$homo <- t2$homo %>% replace_na("not_homo")
t2$length = t2$i.end - t2$i.start
filter(t2, Ds >0.15)

all_pos <- t2 %>% 
  filter(length < 3000) %>%
  #ggplot(., aes(x = i.start, y = Ds, color = homo, group = homo)) +
  ggplot(., aes(x = i.start, y = Ds, color = length)) + #, group = homo)) +
  geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
  facet_wrap(~scaff, scale="free_x") +
  geom_point( size = 1) + 
  #geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  ylim(c(0,1)) +
  xlab("position along chr HD & PR") +
  ylab( expression(italic("Ds"))) +
   th_plot #+ theme(legend.position = "none") 

all_pos 

pdf(file = "Ds_with_homop.pdf", 12,9)
all_pos
dev.off()

################## more quality check:
#awk '$3=="transcript" || $3=="mRNA" ' ../../10.genespace/Mlag12901A1.best/myFile_lociMerged_longestIsoform.gff |cut -f 1-9 |\
#  sed -e 's/ID=//g' |cut -d ";" -f 1 |column -t  > info.mlag1290A1

#awk '$3=="transcript" || $3=="mRNA" ' Mvlag1253A1.gff | grep "MC03\|MC12\|MC16" |\
#  cut -f 1-9 |sed -e 's/ID=//g' |cut -d ";" -f 1 |column -t  > info.mlag1253A1

mlag125 <- read.table("info.mlag1253A1") %>% 
  set_colnames(., c("scaf","tool","inf","start","end","prob","strand","other","geneX")) %>%
  select(-inf, -other, -strand)

mlag129 <- read.table("info.mlag1290A1") %>% 
  set_colnames(., c("scaf129","tool","inf","start129","end129","prob129","strand","other","geneY")) %>%
  select(-inf, -other, -strand)

head(mlag125)
head(t2)
t3 <- merge(mlag125, t2, by = "geneX") %>% select(-start.x, -end.x)
t4 <- merge(mlag129, t3, by = "geneY") %>% mutate(length129 = end129 - start129)

head(t4)
filter(t4, Ds > 0.2)
t4$prob2 <- (as.numeric(t4$prob)) 
t4$prob129.2 <- (as.numeric(t4$prob129)) 

head(t4)
sub <- select(t4, prob2, prob129.2) %>% na.omit()
cor(sub$prob2, sub$prob129.2)

#### look at "true single copy" 
grep -Ff ../orthofinder/Results_Jul19/Orthogroups/Orthogroups_SingleCopyOrthologues.txt MC12_MC16.PR.txt > true_single_copy.PR.txt
grep -Ff ../orthofinder/Results_Jul19/Orthogroups/Orthogroups_SingleCopyOrthologues.txt MC03.txt > true_single_copy.HD.txt
cat true_single_copy.HD.txt true_single_copy.PR.txt > true_single_copy.HD.PR.txt

single <- read.table("true_single_copy.HD.PR.txt") %>% set_colnames(.,c("OG","geneX", "geneY"))
head(single)
#t5 <- merge(t4, single)
#head(Ds_table)
t5 <- merge(Ds_table, single) 
head(t5)
filter(t5, Ds > 0.1)

all_pos2 <- t5 %>% 
  #filter(length < 3000) %>%
  ggplot(., aes(x = start, y = Ds) ) +  #, color = homo, group = homo)) +
  #ggplot(., aes(x = i.start, y = Ds, color = length)) + #, group = homo)) +
  geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
  facet_wrap(~scaff, scale="free_x") +
  geom_point( size = 1) + 
  #geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  ylim(c(0,1)) +
  xlab("position along chr HD & PR") +
  ylab( expression(italic("Ds"))) +
  th_plot #+ theme(legend.position = "none") 

all_pos2

pdf(file = "Ds_mvlag1253_A1_vs_mvlag1253_A2.pdf", 16, 12)
plot_grid(all_pos, all_pos2, labels = c("A: single copy A1-A2 orthologues","B single copy A1-A2-intermedium-lag129A1-Robab"), nrow = 2)
dev.off()
t5 %>% group_by(homo) %>% summarise(mean = mean(Ds))

pdf(file = "Ds_corrected_true_single_copy.pdf",16,9)
all_pos
dev.off()



all_order <- all %>%  
  ggplot(., aes(x = order, y = Ds, color = group, label = gene_name)) +
  geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
  
  geom_point( size = 1) + 
  geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  ylim(c(0,0.9)) +
  xlab("order along chr X") +
  ylab( expression(italic("Ds"))) +
  th_plot + theme(legend.position = "none") +
  scale_color_gradient(low="blue", high="red")

all_order 

pdf(file = "all_genes_26_06.pdf", 12,9)
plot_grid(all_order, all_pos, plottedGenes, nrow = 3)
dev.off()
################################################################################
#         SYNTENY BASED PLOT USING CONICA
################################################################################

#read single copy orthologue between conica and latifolia:
commands <- paste0("awk 'NF==4' Orthogroups/Orthogroups.txt | awk '$3 !~/FUN/ && $3 !~/Chr/ && $3 !~/Scon/ && $4 ~/scaff/ ' |\
                   grep 'chr12' |cut -d " " -f 3-4 |awk '$1 ~/chr/ {print $1"\t"$2"\t"$1"_"$2}' > cp3.txt " )
system(commands)

getwd()
singleCP <- read.table(("/home/quentin/01.analysis/04.silene_latifolia/SILAT/Sconica/orthofinder_26_06/data/cp4.txt"))

#head(singleCP)
#mCP <- read.table("~/03.COLLAB/SILAT/Sconica/orthofinder_26_06/data/Orthmp.txt")
#### creation d'un jeu de données avec plusieurs match possibles ######################################
#awk 'NF<=6 && NF>2 ' Orthogroups/Orthogroups.txt  |    grep "Scon\|FUN\|Chr" |    grep "chr12\|scaff" |\
#  awk '$3 !~/Scon/ && $3 !~/FUN/ && $3 !~/Chr/ ' |\
#  awk 'NF>3 && $3 ~/chr12/ || $4 ~/scaff/ || $5 ~/scaff/ || $6 ~/scaff/ ' |\
#  awk '$3 !~/scaff/ ' | tr -s " " |\
#  sed 's/ /\t/g'  | \
#  awk '{if(NF==4 && $3 ~/chr12/ && $4 ~/scaf/ ) print 
#        else if ( (NF==5 && $5 ~/scaff/) && 
#                    ($3 ~/chr12/  || $4 ~/chr12/ ) ) print 
#            else if ( (NF==6 && $5 ~/scaff/ || $6 ~/scaff/ ) &&
#                     ( $3 ~/chr12/ || $4 ~/chr12/ ) ) print }' |cut -f 3- |perl -pe 's/\t/\n/g' > mp.txt

############################################################################################################
##################################################################################
#### ----------- optional analysis ---- command out to not use it ----------------
all$V3 = paste0(all$geneX, "_", all$geneY)
toto <- all %>%
  filter(V3 %in% singleCP$V4)

#just a check
#toto <- base::merge(all, singleCP, by.x="V3", by.y="V4" )
#head(toto)

#toto1 <- all %>%
#  filter(geneX %in% singleCP$V1 ) %>%
#  filter(geneY %in% singleCP$V1)
na.omit(all$gene_name)
#E707"   "X4"     "X7"     "X3"     "E807"   "DD44X"  "E330"   "Y1"     "SlCypX" "X9"     "SlAP3X" "E162"   "E316"   "X6b" 

na.omit(toto$gene_name)
#"X4"   "X7"   "X3"   "E330" "X9"   "E316"
#we need to recover the 6 other missing gene: 
missing_gene <- c("E707", "E807","DD44X", "Y1", "SlCypX","SlAP3X", "X6b", "E162")
#colnames(toto2)
wanted <- all %>% filter(gene_name %in% missing_gene) 
toto <- rbind(toto, wanted) %>%
  filter(evalue.x < 1e-50 & evalue.y < 1e-50) %>%
  filter(identity.x >75) 


#order the toto to make groups:
toto <- toto[order(i.start), ]
toto$order <- seq(1: nrow(toto))
wanted_l <- trunc(nrow(toto)/20)
part1 <- rep(1: wanted_l ,  each = 20) 
l = length(part1)
part2 <- rep(32, nrow(toto) - l )
toto$group <-  c(part1, part2)

colnames(toto)
################################################################################
################################################################################
#
#                       changepoint analysis    
################################################################################
#https://lindeloev.github.io/mcp/articles/packages.html
library(mcp)

#we reverse the position and order to have it as in previous paper:
toto$i.start <- rev(toto$i.start)
toto$order <- rev(toto$order)

#test based on all data (736 points:)
#model = list(Ds ~ 1, 1~ 1, 1 ~ 1)
#df = select(all, order, Ds) %>% filter(Ds <0.6 ) %>% set_colnames(. , c("chrX", "Ds"))
#fit_mcp = mcp(model, data = df, par_x = "chrX")

#plot(fit_mcp) + ylab(expression(Ds)) + theme_bw() + 
#  geom_point(size = 1, shape = 3) + #,  color = "springgreen") + 
#  plot_pars(fit_mcp, pars = c("cp_1", "cp_2"), type = "dens_overlay") +
#  plot_layout(widths = c(2, 1))

#summary(fit_mcp)

#based only on the 372 genes:
model = list(Ds ~ 1, 1~ 1, 1 ~ 1)

df = toto %>% filter(Ds <1 ) %>%
  #filter(evalue.x < 1e-100 & evalue.y < 1e-100) %>%
  #filter(identity.x > 80) %>% 
  select(order, Ds) %>%
  set_colnames(. , c("chrX", "Ds")) #%>%
#mutate(chrX = rev(chrX))

fit_mcp3 = mcp(model, data = df, par_x = "chrX")
summary(fit_mcp3)

#test the significance of the difference in mean using the bayesian approach:
hypothesis(fit_mcp3, c("int_1 < int_2",
                       "int_2 < int_3"))

pdf(file = "chpt11good_order_run2.pdf", 12, 8)
plot(fit_mcp3, q_fit = TRUE) + 
  ylab(expression(italic(d[S]))) + 
  xlab("order along chr X") + 
  theme_bw() + 
  geom_point(size = 1, shape = 3) + th_plot  #,  color = "springgreen") + 
dev.off()

plot_pars(fit_mcp3, pars = c("cp_1","cp_2")) + 
  plot_pars(fit_mcp3, pars =c("cp_1","cp_2"), type = "hex" )

fitted(fit_mcp3)
predict(fit_mcp3)

toto$strata <- ifelse(toto$order < 97.519, "strata1", 
                      ifelse(toto$order > 334, "strata3","strata2" ))

toto$strata <- ifelse(toto$order < 93.67, "strata1", 
                      ifelse(toto$order > 335, "strata3","strata2" ))

#test significance of difference using tukeyHSD:
summary(fm1 <- aov(Ds ~ strata, data = toto))
TukeyHSD(fm1,  ordered = TRUE)
plot(TukeyHSD(fm1))

table(toto$strata)

########################## make plot now #######################################
library(wesanderson) #to have some customcolors

#look at the single extrem value:
toto %>% filter(Ds >0.45)

toto_order <- toto %>%  
  #filter(evalue.x < 1e-100 & evalue.y < 1e-100) %>%
  #filter(identity.x >80) %>%
  #ggplot(., aes(x = order, y = Ds, color = group, label = gene_name)) +
  ggplot(., aes(x = order, y = Ds, color = strata, label = gene_name)) +
    geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
  
  geom_point( size = 1) + 
  geom_text(hjust = 0, vjust = 0, size = 3, color="black") +
  geom_vline(xintercept =  70, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 156, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 327, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 339, linetype = "dotted", color = "black") +
  
  theme_classic() +
  ylim(c(0,0.9)) +
  xlab("order along chr X") +
  ylab( expression(italic(d[S]))) +
  th_plot + theme(legend.position = "none") +
  scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1"))
  #scale_color_gradient(low="blue", high="red")

toto_order 

#get the value of order to obtain the corresponding position in bp and abline:
toto %>% filter(order == 74 | order == 180 | order == 334 | order == 339)
toto %>% filter(order == 70 | order == 156 | order == 327 | order == 339)

#plot this based on position along the X chromosome now:
toto_pos <- toto %>%  
  ggplot(., aes(x = i.start, y = Ds, color = strata, label = gene_name)) +
  geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
  geom_vline(xintercept =  7284777, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 15161963, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 308912380, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 312879472, linetype = "dotted", color = "black") +
  
  #add the PAR:
  geom_rect(aes(xmin = 0, xmax = 20e6, ymin = 0, ymax = 0.75), fill = "white", alpha = 0, color = "blue", linetype ="dotted"  ) +

  geom_point( size = 1) + 
  #geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  ylim(c(0,0.9)) +
  xlab("position along chr X") +
  ylab( expression(italic(d[S]))) +
  th_plot + theme(legend.position = "none") +
  scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1"))

toto_pos 

pdf(file = "single_copy_ortho_11_07_PAR_v3.pdf", 12,6)
plot_grid(toto_pos, toto_order, nrow = 2, labels = "AUTO")
dev.off()

head(toto)

############### coloring by Y ##########################"
toto_pos2 <- toto %>%  
  tidyr::separate(geneY, into = c(NA, "scaff","id"),  sep = "_") %>%
  mutate(Yname = paste0(scaff,"_", id)) %>%
  ggplot(., aes(x = i.start, y = Ds, color = Yname, label = gene_name)) +
  geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
  geom_vline(xintercept =  7284777, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 15161963, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 308912380, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 312879472, linetype = "dotted", color = "black") +
  
  #add the PAR:
  #geom_rect(aes(xmin = 0, xmax = 20e6, ymin = 0, ymax = 0.75), fill = "white", alpha = 0, color = "blue", linetype ="dotted"  ) +
  
  geom_point( size = 1) + 
  #geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  ylim(c(0,0.9)) +
  xlab("position along chr X") +
  ylab( expression(italic(d[S]))) +
  th_plot + 
  #theme(legend.position = "none") +
  #scale_color_manual(values=wes_palette(n=5, name="Darjeeling1"))
  scale_color_manual(values=wes_palette(n=5, name="Darjeeling2"))
  #scale_color_manual(values=wes_palette(n=5, name="BottleRocket2"))

toto_pos2 
pdf(file="Ds_colored_onY.pdf", 16,7)
toto_pos2
dev.off()


#################################################################################
tannier <- select(toto, geneX, geneY, V1, i.start, i.end, gene_name, length, strata)
head(tannier)
na.omit(tannier)
tannier <- select(toto, geneX, geneY, V1, i.start, i.end, length, strata)

write.table(tannier, "single_copy_orthologue_with_tentative_strata.txt", quote = F, row.names = F)

############## Dn and Dn/Ds plot for the X #####################################
########################## make plot now #######################################
cor(toto$Dn, toto$Ds)
summary(lm(toto$Dn ~ toto$Ds))
head(toto)
toto %>% group_by(strata) %>% summarise(meanDn = mean(Dn), meanSE = mean(SEDn))

summary(fm2 <- aov(Dn ~ strata, data = toto))
TukeyHSD(fm2,  ordered = TRUE)
plot(TukeyHSD(fm1))


dn_pos <- toto %>%  
  ggplot(., aes(x = i.start, y = Dn, color = strata, label = gene_name)) +
  geom_errorbar(aes(ymin = Dn-SEDn, ymax = Dn + SEDn), width = .1) +
  
  geom_point( size = 1) + 
  geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  #ylim(c(0,0.9)) +
  xlab("position along chr X") +
  ylab( expression(italic("Dn"))) +
  th_plot + theme(legend.position = "none")  +
  scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1"))


dn_pos 

dn_order <- toto %>%  
  ggplot(., aes(x = order, y = Dn, color = strata, label = gene_name)) +
  geom_errorbar(aes(ymin = Dn-SEDn, ymax = Dn + SEDn), width = .1) +
  
  geom_point( size = 1) + 
  geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  #ylim(c(0,0.9)) +
  xlab("order along chr X") +
  ylab( expression(italic("Dn"))) +
  th_plot + theme(legend.position = "none") +
  scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1"))
  #scale_color_gradient(low="blue", high="red")

dn_order 

dnds_order <- toto %>%  
  ggplot(., aes(x = order, y = Dn/Ds, color = strata, label = gene_name)) +
  #geom_errorbar(aes(ymin = Ds-SEDn, ymax = Ds + SEDn), width = .1) +
  
  geom_point( size = 1) + 
  geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  ylim(c(0,1)) +
  xlab("order along chr X") +
  ylab( expression(italic("Dn/Ds"))) +
  th_plot + theme(legend.position = "none") +
  scale_color_manual(values=wes_palette(n=3, name="GrandBudapest1"))
  #scale_color_gradient(low="blue", high="red")

dnds_order 

pdf(file = "DnDs_v2.pdf", 12,8)
plot_grid(dn_pos, dn_order, dnds_order, nrow = 3, labels = "AUTO")
dev.off()

















################################## Deprecated #####################################
########## ----- TEsT on MULTICOPY ---- #######################################
toto <- filter(tmp, Ds < 1000) %>% 
  #filter(prob > 0.15) %>%
  select(geneY, geneX,  i.start, i.end, gene_name, Ds, SEDs, evalue.x, evalue.y, identity.x, length) %>% 
  distinct(.) %>% 
  mutate(Diff = Ds-SEDs) %>% 
  mutate(gene_name = str_replace_all(gene_name, "Gene.", NA_character_ ) ) %>%
  distinct(.) %>%
  filter(geneX %in% mCP$V1 | geneY %in% mCP$V1)

#order the toto to make groups:
toto <- toto[order(i.start), ]
toto$order <- seq(1: nrow(toto))
wanted_l <- trunc(nrow(toto)/20)
part1 <- rep(1: wanted_l ,  each = 20) 
#part1 <- rep(1:32, each = 20) 
l = length(part1)
part2 <- rep(32, nrow(toto) - l )
toto$group <-  c(part1, part2)

########################## make plot now #######################################
toto_pos <- toto %>%  
  ggplot(., aes(x = i.start, y = Ds, label = gene_name)) +
  geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
  
  geom_point( size = 1) + 
  geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  ylim(c(0,0.9)) +
  xlab("position along chr X") +
  ylab( expression(italic("Ds"))) +
  th_plot + theme(legend.position = "none") #+
#scale_color_gradient(low="blue", high="red")

toto_pos 

toto_order <- toto %>%  
  ggplot(., aes(x = order, y = Ds, color = group, label = gene_name)) +
  geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
  
  geom_point( size = 1) + 
  geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  ylim(c(0,0.9)) +
  xlab("order along chr X") +
  ylab( expression(italic("Ds"))) +
  th_plot + theme(legend.position = "none") +
  scale_color_gradient(low="blue", high="red")

toto_order 


pdf(file = "Ds_2-4copy_ortho_26_06.pdf", 12,6)
plot_grid(toto_order, toto_pos,  nrow = 2)
dev.off()

#########################################################################################################################################################


gameto <- read.table("gametologs", h =T)
head(gameto)
head(Ds_table)

#### plot based on the gametolog ##############"
#toto_pos <- toto %>%  
D <- D %>% tidyr::separate_wider_delim(cols = geneX, names = c("geneX",NA), delim = ".")
mCP2 <- mCP %>%  tidyr::separate_wider_delim(cols = V1, names = c("geneX",NA), delim = ".")   

totoG <- Ds_table %>%  
  tidyr::separate_wider_delim(cols = geneX, names = c("geneX",NA), delim = ".") %>%
  filter(geneX %in% gameto$name_gene_X) %>% 
  rbind(., D) %>%
  distinct(.) %>%
  filter(geneX %in% mCP2$geneX)


totoG <- totoG[order(totoG$i.start), ]
totoG$order <- seq(1: nrow(totoG))
wanted_l <- trunc(nrow(totoG)/20)
part1 <- rep(1: wanted_l ,  each = 20) 
#part1 <- rep(1:32, each = 20) 
l = length(part1)
part2 <- rep(32, nrow(totoG) - l )
totoG$group <-  c(part1, part2)

#https://stackoverflow.com/questions/53289315/r-calculating-rolling-average-with-window-based-on-value-not-number-of-rows-or

toto_gameto <- filter(totoG, length > 1) %>%  
  ggplot(., aes(x = order, y = Ds, color = group, label = name)) +
  geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
  
  geom_point( size = 1) + 
  geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  ylim(c(0,0.9)) +
  xlab("order along chr X") +
  ylab( expression(italic("Ds"))) +
  th_plot + theme(legend.position = "none") +
  scale_color_gradient(low="blue", high="red")

toto_gameto 

toto_gameto_pos <- filter(totoG, length > 1) %>%  
  ggplot(., aes(x = i.start, y = Ds, label = name)) +
  geom_errorbar(aes(ymin = Ds-SEDs, ymax = Ds + SEDs), width = .1) +
  
  geom_point( size = 1) + 
  geom_text(hjust = 0, vjust = 0, size = 5, color="black") +
  theme_classic() +
  ylim(c(0,0.9)) +
  xlab("position along chr X") +
  ylab( expression(italic("Ds"))) +
  th_plot + theme(legend.position = "none") #+
  #scale_color_gradient(low="blue", high="red")
toto_gameto_pos


pdf(file = "Ds_MultyCopy_gametolog_and_known_genes_order_26_06.pdf", 12,8)
plot_grid(toto_gameto, toto_gameto_pos, plottedGenes, nrow = 3)
dev.off()


pdf(file = "Ds_gametolog_and_known_genes_order_26_06.pdf", 12,6)
plot_grid(toto_order, plottedGenes, nrow = 2)
dev.off()

