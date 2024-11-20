
################################################################################
######  Script to plot Ds values and perform changepoint analyses     ##########
#Author: QR
#Date: 26-01-24
################################################################################
### load librariers ############################################################

# perform the usual verifications :
if("mcp" %in% rownames(installed.packages()) == FALSE)
{install.packages("mcp", repos="https://cloud.r-project.org") }
if("dplyr" %in% rownames(installed.packages()) == FALSE)
{install.packages("dplyr", repos="https://cloud.r-project.org") }
if("ggplot2" %in% rownames(installed.packages()) == FALSE)
{install.packages("ggplot2", repos="https://cloud.r-project.org") }
if("cowplot" %in% rownames(installed.packages()) == FALSE)
{install.packages("cowplot", repos="https://cloud.r-project.org") }


library(mcp)
library(ggplot2)
library(cowplot)
library(dplyr)

#---- load data ---- # 
df <- read.table("02_results/dS.values.forchangepoint.txt", h = T) #a table with two columns : Ds and order 
df <- filter(df, Ds < 0.3) 
#this input is normally produced from the script (03.plot_paml.R)


################################################################################
#                   perform the changepoint analyis here: 
################################################################################
#define the model we want to test:
model8strata = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1, 1 ~ 1 , 1 ~ 1 , 1 ~ 1, 1  ~ 1 ) 
model7strata = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1, 1 ~ 1 , 1 ~ 1 , 1 ~ 1)
model6strata = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1, 1 ~ 1 , 1 ~ 1)
model5strata = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1, 1 ~ 1)
model4strata = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1)
model3strata = list(Ds ~ 1, 1~ 1, 1 ~ 1)
model2strata = list(Ds ~ 1, 1~ 1)
modelsimple = list(Ds ~ 1 + orderchp)

fit_1st = mcp(modelsimple,  data = df, par_x = "orderchp",   iter = 8e3, adapt = 1.5e3,  chains = 5, cores = 5 )
fit_2st = mcp(model2strata, data = df, par_x = "orderchp",   iter = 8e3, adapt = 1.5e3,  chains = 5, cores = 5 )
fit_3st = mcp(model3strata, data = df, par_x = "orderchp",   iter = 8e3, adapt = 1.5e3,  chains = 5, cores = 5 )
fit_4st = mcp(model4strata, data = df, par_x = "orderchp",   iter = 8e3, adapt = 1.5e3,  chains = 5, cores = 5 )
fit_5st = mcp(model5strata, data = df, par_x = "orderchp",   iter = 8e3, adapt = 1.5e3,  chains = 5, cores = 5 )
fit_6st = mcp(model6strata, data = df, par_x = "orderchp",   iter = 8e3, adapt = 1.5e3,  chains = 5, cores = 5 )
fit_7st = mcp(model7strata, data = df, par_x = "orderchp",   iter = 8e3, adapt = 1.5e3,  chains = 5, cores = 5 )
fit_8st = mcp(model8strata, data = df, par_x = "orderchp",   iter = 8e3, adapt = 1.5e3,  chains = 5, cores = 5 )


#use my usual theme:
## ------------------ GGPLOT  CUSTOMISATION ------------------------------------------------##
th_plot <-  theme(axis.title.x=element_text(size=14, family="Helvetica",face="bold"),
    axis.text.x=element_text(size=14,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
    axis.title.y=element_text(size=18, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
    axis.text.y=element_text(size=14,family="Helvetica",face="bold"),
    strip.text.x = element_text(size=18),panel.grid.major = element_blank())


#make simple plot of model and fits:
xl <- expression(paste("order along ", italic("ancestral species"), " mating chromsome"))
FigMCP4Strata  <- plot(fit_4st, q_fit = TRUE) + ggplot2::ggtitle("Posterior fit - 4 strata") + theme_classic() + th_plot +
  geom_point(color = "darkblue", size = 0.1) +
  xlab(xl) + ylab(expression(italic(d[s])))

FigMCP3Strata  <- plot(fit_3st, q_fit = TRUE) + ggplot2::ggtitle("Posterior fit - 3 strata") + theme_classic() + th_plot +
  geom_point(color = "darkblue", size = 0.21) +
  xlab(xl) + ylab(expression(italic(d[s])))

FigMCP2Strata  <- plot(fit_2st, q_fit = TRUE) + ggplot2::ggtitle("Posterior fit - 2 strata") + theme_classic() + th_plot +
  geom_point(color = "darkblue",  size = 0.21) +
  xlab(xl) + ylab(expression(italic(d[s])))

FigMCP5Strata  <- plot(fit_5st, q_fit = TRUE) + ggplot2::ggtitle("Posterior fit - 5 strata") + theme_classic() + th_plot +
  geom_point(color = "darkblue", size = 0.1) +
  xlab(xl) + ylab(expression(italic(d[s])))

FigMCP6Strata  <- plot(fit_6st, q_fit = TRUE) + ggplot2::ggtitle("Posterior fit - 6 strata") + theme_classic() + th_plot +
  geom_point(color = "darkblue", size = 0.1) +
  xlab(xl) + ylab(expression(italic(d[s])))

FigMCP7Strata  <- plot(fit_7st, q_fit = TRUE) + ggplot2::ggtitle("Posterior fit - 7 strata") + theme_classic() + th_plot +
  geom_point(color = "darkblue", size = 0.1) +
  xlab(xl) + ylab(expression(italic(d[s])))

FigMCP8Strata  <- plot(fit_8st, q_fit = TRUE) + ggplot2::ggtitle("Posterior fit - 8 strata") + theme_classic() + th_plot +
  geom_point(color = "darkblue", size = 0.1) +
  xlab(xl) + ylab(expression(italic(d[s])))


#create dir if not present:
if (!dir.exists("02_results/modelcomp")){
  dir.create("02_results/modelcomp")
}


#plot them all
pdf(file = "02_results/modelcomp/Strata_comp_XX.pdf", 12,18)
plot_grid(FigMCP2Strata, FigMCP3Strata, FigMCP4Strata, 
          FigMCP5Strata, FigMCP6Strata, 
          FigMCP7Strata, FigMCP8Strata,
          labels = "AUTO", ncol = 1)
dev.off()

pdf(file = "02_results/modelcomp/Strata_comp_short.pdf", 12,18)
plot_grid(
          FigMCP5Strata, FigMCP6Strata, 
          FigMCP7Strata, FigMCP8Strata,
          labels = "AUTO", ncol = 1)
dev.off()

#summarise and look at point change:
m1 <- summary(modelsimple)
m2 <- summary(fit_2st)
m3 <- summary(fit_3st)
m4 <- summary(fit_4st)
m5 <- summary(fit_5st)
m6 <- summary(fit_6st)
m7 <- summary(fit_7st)
m8 <- summary(fit_8st)

## make some plot based on mean changes:
#getting mean-value along the X of the first changepoint 
#2-changepoint model:
cp1m2 <- m2$mean[1]

#3-changepoint :
cp1m3 <- m3$mean[1]
cp2m3 <- m3$mean[2]

#4 - changepoint model :
cp1m4 <- m4$mean[1]
cp2m4 <- m4$mean[2]
cp3m4 <- m4$mean[3]

#5 - changepoint model :
cp1m5 <- m5$mean[1] ##alternatively: summary(fit_5st)[2][[1]][1]
cp2m5 <- m5$mean[2]
cp3m5 <- m5$mean[3]
cp4m5 <- m5$mean[4]

#5 - changepoint model :
cp1m6 <- m6$mean[1] ##alternatively: summary(fit_5st)[2][[1]][1]
cp2m6 <- m6$mean[2]
cp3m6 <- m6$mean[3]
cp4m6 <- m6$mean[4]
cp5m6 <- m6$mean[5]


#5 - changepoint model :
cp1m7 <- m7$mean[1] ##alternatively: summary(fit_5st)[2][[1]][1]
cp2m7 <- m7$mean[2]
cp3m7 <- m7$mean[3]
cp4m7 <- m7$mean[4]
cp5m7 <- m7$mean[5]
cp6m7 <- m7$mean[5]


#5 - changepoint model :
cp1m8 <- m8$mean[1] ##alternatively: summary(fit_5st)[2][[1]][1]
cp2m8 <- m8$mean[2]
cp3m8 <- m8$mean[3]
cp4m8 <- m8$mean[4]
cp5m8 <- m8$mean[5]
cp6m8 <- m8$mean[5]
cp7m8 <- m8$mean[5]

#filter data and visualise differences in mean Ds values:
#note: we could do better using the BF and hypothesis function instead of working with means ....

df$two_strata <- ifelse(df$orderchp < cp1m2, "strata1", "strata2")
df$three_strata <- ifelse(df$orderchp < cp1m3, "strata1", 
                    ifelse(df$orderchp > cp2m3, "strata3", "strata2"))

df$four_strata <- ifelse(df$orderchp < cp1m4, "strata1", 
                    ifelse(df$orderchp > cp3m4, "strata4",
                    ifelse(df$orderchp > cp2m4 & df$orderchp < cp3m4, "strata3", "strata2")))


df$five_strata <- ifelse(df$orderchp < cp1m5, "strata1", 
                    ifelse(df$orderchp > cp4m5, "strata5",
                    ifelse(df$orderchp > cp1m5 & df$orderchp < cp2m5, "strata2",
                    ifelse(df$orderchp > cp2m5 & df$orderchp < cp3m5, "strata3","strata4"))))


df$six_strata <- ifelse(df$orderchp < cp1m6, "strata1", 
                    ifelse(df$orderchp > cp5m6, "strata6",
                    ifelse(df$orderchp > cp1m6 & df$orderchp < cp2m6, "strata2",
                    ifelse(df$orderchp > cp2m6 & df$orderchp < cp3m6, "strata3",
                    ifelse(df$orderchp > cp3m6 & df$orderchp < cp4m6, "strata4", "strata5" )))))


df$set_strata <- ifelse(df$orderchp < cp1m7, "strata1", 
                    ifelse(df$orderchp > cp6m7, "strata7",
                    ifelse(df$orderchp > cp1m7 & df$orderchp < cp2m7, "strata2",
                    ifelse(df$orderchp > cp2m7 & df$orderchp < cp3m7, "strata3",
                    ifelse(df$orderchp > cp3m7 & df$orderchp < cp4m7, "strata4",
                    ifelse(df$orderchp > cp4m7 & df$orderchp < cp5m7, "strata5","strata6"))))))
#                    ifelse(df$orderchp > cp5m7 & df$orderchp < cp6m7, "strata6",
#                           "strata7" )))))))

df$huit_strata <- ifelse(df$orderchp < cp1m8, "strata1", 
                    ifelse(df$orderchp > cp7m8, "strata8",
                    ifelse(df$orderchp > cp1m8 & df$orderchp < cp2m8, "strata2",
                    ifelse(df$orderchp > cp2m8 & df$orderchp < cp3m8, "strata3",
                    ifelse(df$orderchp > cp3m8 & df$orderchp < cp4m8, "strata4",
                    ifelse(df$orderchp > cp4m8 & df$orderchp < cp5m8, "strata5",
                    ifelse(df$orderchp > cp5m8 & df$orderchp < cp6m8, "strata6",
                           "strata7" )))))))

write.table(df,"02_results/modelcomp/df.txt",quote=F,row.names=F,col.names=T,sep="\t")

mycolor <-c("#E69F00",  "#0072B2" ,"#5B1A79",  "#CC79A7", "#D55E00", "red", "black","yellow")

plot2cp <- ggplot(df, aes(x = two_strata, y = Ds, fill = two_strata)) + 
    geom_violin(trim = FALSE) + 
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    theme_classic() + 
    th_plot  + 
    ylab("dS") +
    scale_fill_manual(values=mycolor[1:2])  + theme(legend.position="none")

plot3cp <- ggplot(df, aes(x = three_strata, y = Ds, fill = three_strata)) + 
    geom_violin(trim = FALSE) + 
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    theme_classic() + 
    th_plot  + 
    ylab("dS") +
    scale_fill_manual(values=mycolor[1:3])  + theme(legend.position="none")


plot4cp <- ggplot(df, aes(x = four_strata, y = Ds, fill = four_strata)) + 
    geom_violin(trim = FALSE) + 
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    theme_classic() + 
    th_plot  + 
    ylab("dS") +
    scale_fill_manual(values=mycolor[1:4])  + theme(legend.position="none")


plot5cp <- ggplot(df, aes(x = five_strata, y = Ds, fill = five_strata)) + 
    geom_violin(trim = FALSE) + 
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    theme_classic() + 
    th_plot  + 
    ylab("dS") +
    scale_fill_manual(values=mycolor[1:5])  + theme(legend.position="none")


plot6cp <- ggplot(df, aes(x = six_strata, y = Ds, fill = six_strata)) + 
    geom_violin(trim = FALSE) + 
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    theme_classic() + 
    th_plot  + 
    ylab("dS") +
    scale_fill_manual(values=mycolor[1:6])  + theme(legend.position="none")


plot7cp <- ggplot(df, aes(x = set_strata, y = Ds, fill = set_strata)) + 
    geom_violin(trim = FALSE) + 
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    theme_classic() + 
    th_plot  + 
    ylab("dS") +
    scale_fill_manual(values=mycolor[1:7])  + theme(legend.position="none")


plot8cp <- ggplot(df, aes(x = huit_strata, y = Ds, fill = huit_strata)) + 
    geom_violin(trim = FALSE) + 
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    theme_classic() + 
    th_plot  + 
    ylab("dS") +
    scale_fill_manual(values=mycolor[1:8])  + theme(legend.position="none")

#plot them all
pdf(file = "02_results/modelcomp/Strata_comp4.pdf", 12,18)
plot_grid(plot2cp, plot3cp, plot4cp, plot5cp , plot6cp , 
          plot7cp, plot8cp, labels = "AUTO", ncol = 1)
dev.off()

#----------------  a bit of model choice --------------------------------------#
#to be considered precautionously!
fit_1st$loo <- loo(fit_1st)
fit_2st$loo <- loo(fit_2st)
fit_3st$loo <- loo(fit_3st)
fit_4st$loo <- loo(fit_4st)
fit_5st$loo <- loo(fit_5st)
fit_6st$loo <- loo(fit_6st)
fit_7st$loo <- loo(fit_7st)
fit_8st$loo <- loo(fit_8st)

m.choice <- loo::loo_compare(fit_1st$loo, fit_2st$loo, fit_3st$loo, 
                             fit_4st$loo, fit_5st$loo, fit_6st$loo,
                             fit_7st$loo, fit_8st$loo)

# ----- rule of thumb - to be used as a helper for model choice ---- : #
# you should read the manual (https://lindeloev.github.io/mcp/articles/comparison.html)
# and associated paper to help in choosing a model.

#finally compute the weight:
loo_list = list(fit_1st$loo, fit_2st$loo, fit_3st$loo, 
                fit_4st$loo, fit_5st$loo, fit_6st$loo)

weights <- loo::loo_model_weights(loo_list, method="pseudobma")

write.table(weights,"02_results/modelcomp/model_weights",quote=F, col.names=("weights"))


# ----- some hypothesis testing regarding differences among intervals -------- #
#testing hypothesis
# interval1 lower than interval2:
h2st <- hypothesis(fit_2st, c("int_1 < int_2"))
#the reverse:
h2strev <- hypothesis(fit_2st, c("int_2 < int_1"))

#3-model test:
#fitting traditionnal expectation of a gradual loss :
h3st1 <-  hypothesis(fit_3st, c("int_1 > int_2", "int_2 > int_3" ))
#and on the other way around:
h3st2 <-hypothesis(fit_3st, c("int_1 < int_2", "int_2 < int_3" ))
#more convoluted possibility:
h3st3 <-hypothesis(fit_3st, c("int_1 > int_2", "int_2 < int_3" ))
h3st4 <-hypothesis(fit_3st, c("int_1 < int_2", "int_2 > int_3" ))

#this could be done with 4 and 5 strata. 
#However we invite the user to perform is own investigation.

#see more here: https://lindeloev.github.io/mcp/articles/comparison.html


#----------------save the data -------------------------------------------------#
#summaries:
write.table(m1, "02_results/modelcomp/summary1strata.txt", quote = F)
write.table(m2, "02_results/modelcomp/summary2strata.txt", quote = F)
write.table(m3, "02_results/modelcomp/summary3strata.txt", quote = F)
write.table(m4, "02_results/modelcomp/summary4strata.txt", quote = F)
write.table(m5, "02_results/modelcomp/summary5strata.txt", quote = F)
write.table(m6, "02_results/modelcomp/summary6strata.txt", quote = F)
write.table(m7, "02_results/modelcomp/summary7strata.txt", quote = F)
write.table(m8, "02_results/modelcomp/summary8strata.txt", quote = F)

#model choice attempt: 
write.table(m.choice ,"02_results/modelcomp/model.choice.txt", quote = F)

#hypothesis testing attempt: 
write.table(h2st,  "02_results/modelcomp/hypothesis2strata", quote= F)
write.table(h2st,  "02_results/modelcomp/hypothesis2strata.rev.txt", quote= F)
write.table(h3st1, "02_results/modelcomp/hypothesis3strata.1.txt", quote= F)
write.table(h3st2, "02_results/modelcomp/hypothesis3strata.2.txt", quote= F)
write.table(h3st3, "02_results/modelcomp/hypothesis3strata.3.txt", quote= F)
write.table(h3st4, "02_results/modelcomp/hypothesis3strata.4.txt", quote= F)

save.image( file = "02_results/modelcomp/changepoint_analysis.RData")
#save(fit_1st, fit_2st,fit_3st, fit_4st, fit_5st, file = "02_results/modelcomp/changepoint_analysis.RData")


th_plot2 <-  theme(axis.title.x=element_text(size=8, family="Helvetica",face="bold"),
    axis.text.x=element_text(size=8,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
    axis.title.y=element_text(size=8, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
    axis.text.y=element_text(size=8,family="Helvetica",face="bold"),
    strip.text.x = element_text(size=7),panel.grid.major = element_blank())

dplot <- function(df_of_ds, nstrata, columnstrata) {
    columnstrata=sym(columnstrata)
    ggplot(df_of_ds, aes(x = start, y = Ds, colour = !!columnstrata)) + 
    geom_point( size = 1) + 
    facet_wrap(~scaff, scale="free_x") +
    theme_classic() +
    ylim(c(0,0.5)) +
    xlab("position along chr (bp)") +
    ylab( expression(italic("dS"))) +
    th_plot2 + 
    theme(legend.position = "none") + 
    #scale_colour_manual(values=mycolor[1:4])  
    scale_colour_manual(values=mycolor[1:nstrata])  

}

ds3 <- dplot(df, nstrata=3, "three_strata")
ds4 <- dplot(df, nstrata=4, "four_strata")
ds5 <- dplot(df, nstrata=5, "five_strata")
ds6 <- dplot(df, nstrata=6, "six_strata")
ds7 <- dplot(df, nstrata=7, "set_strata")
ds8 <- dplot(df, nstrata=8, "huit_strata")


pdf(file="02_results/modelcomp/plot_all_ds.pdf",8,12)
plot_grid(ds3, ds4, ds5, ds6, ds7, ds8, 
          labels = c("A - three changepoint",
                     "B - four changepoint" ,
                     "C - five changepoint" ,
                     "D - six changepoint"  ,
                     "E - seven changepoint" ,
                     "F - eigth changepoint") , 
         label_size = 7,
         hjust = -0.5, vjust = -0.5,
         ncol = 1)
dev.off()


#finally: 
s3.anc.h1 <- select(df, gene, geneX, three_strata)
s4.anc.h1 <- select(df, gene, geneX, four_strata)
s5.anc.h1 <- select(df, gene, geneX, five_strata)
s6.anc.h1 <- select(df, gene, geneX, five_strata)
s6.anc.h1 <- select(df, gene, geneX, six_strata)
s7.anc.h1 <- select(df, gene, geneX, set_strata)
s8.anc.h1 <- select(df, gene, geneX, huit_strata)

s3.h1.h2 <- select(df, geneX, geneY.x, three_strata)
s4.h1.h2 <- select(df, geneX, geneY.x, four_strata)
s5.h1.h2 <- select(df, geneX, geneY.x, five_strata)
s6.h1.h2 <- select(df, geneX, geneY.x, five_strata)
s6.h1.h2 <- select(df, geneX, geneY.x, six_strata)
s7.h1.h2 <- select(df, geneX, geneY.x, set_strata)
s8.h1.h2 <- select(df, geneX, geneY.x, huit_strata)


write.table(s3.anc.h1,"02_results/modelcomp.s3.anc.h1",quote=F,row.names=F,col.names=F,sep="\t")
write.table(s4.anc.h1,"02_results/modelcomp.s4.anc.h1",quote=F,row.names=F,col.names=F,sep="\t")
write.table(s5.anc.h1,"02_results/modelcomp.s5.anc.h1",quote=F,row.names=F,col.names=F,sep="\t")
write.table(s6.anc.h1,"02_results/modelcomp.s6.anc.h1",quote=F,row.names=F,col.names=F,sep="\t")
write.table(s7.anc.h1,"02_results/modelcomp.s7.anc.h1",quote=F,row.names=F,col.names=F,sep="\t")
write.table(s8.anc.h1,"02_results/modelcomp.s8.anc.h1",quote=F,row.names=F,col.names=F,sep="\t")

write.table(s3.h1.h2,"02_results/modelcomp.s3.h1.h2",quote=F,row.names=F,col.names=F,sep="\t")
write.table(s4.h1.h2,"02_results/modelcomp.s4.h1.h2",quote=F,row.names=F,col.names=F,sep="\t")
write.table(s5.h1.h2,"02_results/modelcomp.s5.h1.h2",quote=F,row.names=F,col.names=F,sep="\t")
write.table(s6.h1.h2,"02_results/modelcomp.s6.h1.h2",quote=F,row.names=F,col.names=F,sep="\t")
write.table(s7.h1.h2,"02_results/modelcomp.s7.h1.h2",quote=F,row.names=F,col.names=F,sep="\t")
write.table(s8.h1.h2,"02_results/modelcomp.s8.h1.h2",quote=F,row.names=F,col.names=F,sep="\t")

