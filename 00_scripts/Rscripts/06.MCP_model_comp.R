
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
df <- read.table("dS.values.forchanepoint.txt", h = T) #a table with two columns : Ds and ordre 

#this input is normally produced from the script (03.plot_paml.R)


################################################################################
#                   perform the changepoint analyis here: 
################################################################################
#define the model we want to test:
model5strata = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1, 1 ~ 1)
model4strata = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1)
model3strata = list(Ds ~ 1, 1~ 1, 1 ~ 1)
model2strata = list(Ds ~ 1, 1~ 1)
modelsimple = list(Ds ~ 1 + ordre)

fit_1st = mcp(modelsimple,  data = df, par_x = "ordre",   iter = 8e3, adapt = 1.5e3,  chains = 5, cores = 5 )
fit_2st = mcp(model2strata, data = df, par_x = "ordre",   iter = 8e3, adapt = 1.5e3,  chains = 5, cores = 5 )
fit_3st = mcp(model3strata, data = df, par_x = "ordre",   iter = 8e3, adapt = 1.5e3,  chains = 5, cores = 5 )
fit_4st = mcp(model4strata, data = df, par_x = "ordre",   iter = 8e3, adapt = 1.5e3,  chains = 5, cores = 5 )
fit_5st = mcp(model5strata, data = df, par_x = "ordre",   iter = 8e3, adapt = 1.5e3,  chains = 5, cores = 5 )


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

#plot them all
pdf(file = "Strata_comp.pdf", 12,14)
plot_grid(FigMCP2Strata, FigMCP3Strata, FigMCP4Strata, FigMCP5Strata, labels = "AUTO", ncol = 1)
dev.off()

#summarise and look at point change:
m1 <- summary(modelsimple)
m2 <- summary(fit_2st)
m3 <- summary(fit_3st)
m4 <- summary(fit_4st)
m5 <- summary(fit_5st)

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
cp4m5 <- m5$mean[3]

#filter data and visualise differences in mean Ds values:
#note: we could do better using the BF and hypothesis function instead of working with means ....

df$two_strata <- ifelse(df$ordre < cp1m2, "strata1", "strata2")
df$three_strata <- ifelse(df$ordre < cp1m3, "strata1", 
                    ifelse(df$ordre > cp2m3, "strata3", "strata2"))

df$four_strata <- ifelse(df$ordre < cp1m4, "strata1", 
                    ifelse(df$ordre > cp3m4, "strata4",
                    ifelse(df$ordre > cp2m4 & df$ordre < cp3m4, "strata3", "strata2")))


df$five_strata <- ifelse(df$ordre < cp1m5, "strata1", 
                    ifelse(df$ordre > cp4m5, "strata5",
                    ifelse(df$ordre > cp1m5 & df$ordre < cp2m5, "strata2",
                    ifelse(df$ordre > cp2m5 & df$ordre < cp3m5, "strata3","strata4"))))


mycolor <-c("#E69F00",  "#0072B2" ,"#5B1A79",  "#CC79A7", "#D55E00")

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


#----------------  a bit of model choice --------------------------------------#
#to be considered precautionously!
fit_2st$loo <- loo(fit_2st)
fit_3st$loo <- loo(fit_3st)
fit_4st$loo <- loo(fit_4st)
fit_5st$loo <- loo(fit_5st)

m.choice <- loo::loo_compare(fit_2st$loo, fit_3st$loo, fit_4st$loo fit_5st$loo)

# ----- rule of thumb - to be used as a helper for model choice ---- : #
# you should read the manual (https://lindeloev.github.io/mcp/articles/comparison.html)
# and associated paper to help in choosing a model.

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
write.table(m1, "summary_1strata.txt", quote = F)
write.table(m2, "summary_2strata.txt", quote = F)
write.table(m3, "summary_3strata.txt", quote = F)
write.table(m4, "summary_4strata.txt", quote = F)
write.table(m5, "summary_5strata.txt", quote = F)

#model choice attempt: 
write.table(m.choice ,"model.choice.txt", quote = F)

#hypothesis testing attempt: 
write.table(h2st,"hypothesis2strata", quote= F)
write.table(h2st,"hypothesis2strata.rev.txt", quote= F)
write.table(h3st1,"hypothesis3strata.1.txt", quote= F)
write.table(h3st2,"hypothesis3strata.2.txt", quote= F)
write.table(h3st3,"hypothesis3strata.3.txt", quote= F)
write.table(h3st4,"hypothesis3strata.4.txt", quote= F)

save(fit_1st, fit_2st,fit_3st, fit_4st, fit_5st, file = "changepoint_analysis.RData")
