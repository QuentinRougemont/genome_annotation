
################################################################################
######  Script to plot and filter Ds values based on various metadata ##########
#Author: QR
#Date: 12-05-23
###############################################################################
### load library ################################################################
#library(dplyr)
#library(magrittr)
#library(data.table)
library(mcp)


#---- load data ---- # 
df <- read.table("dS.modelcomp", h = T) #a table with two columns : Ds and ordre #hey should be named like this!

################################################################################
#                   perform the changepoint analyis here: 
################################################################################
#define the model we want to test:
model4strata = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1)
model5strata = list(Ds ~ 1, 1~ 1, 1 ~ 1, 1 ~ 1)
model3strata = list(Ds ~ 1, 1~ 1, 1 ~ 1)
model2strata = list(Ds ~ 1, 1~ 1)
modelsimple = list(Ds ~ 1 + ordre)

fit_2st = mcp(model2strata, data = df, par_x = "ordre",   iter = 8e4, adapt = 15e3,  chains = 5, cores = 5 )
fit_3st = mcp(model3strata, data = df, par_x = "ordre",   iter = 8e4, adapt = 15e3,  chains = 5, cores = 5 )
fit_1st = mcp(modelsimple,  data = df, par_x = "ordre",   iter = 8e4, adapt = 15e3,  chains = 5, cores = 5 )
fit_4st = mcp(model4strata, data = df, par_x = "ordre",   iter = 8e4, adapt = 15e3,  chains = 5, cores = 5 )

#test which model is the best here: 

#summarise and look at point change:
summary(fit_4st)
summary(fit_3st)
summary(fit_2st)
summary(modelsimple)
#test hypothesis relative to the change;
hypothesis(fit_4st, c("int_1 < int_2", "int_2 < int_3", "int_3 > int_4"))
hypothesis(fit_3st, c("int_1 < int_2", "int_2 > int_3" ))
hypothesis(fit_2st, c("int_1 < int_2"))

fit_2st$loo <- loo(fit_2st)
fit_3st$loo <- loo(fit_3st)
fit_4st$loo <- loo(fit_4st)

fit_5st = mcp(model5strata, data = df, par_x = "ordre",   iter = 1e5, adapt = 15e3,  chains = 5, cores = 5 )


#stor eit for later
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

#plot themall
pdf(file = "Strata_comp.pdf", 12,14)
plot_grid(FigMCP2Strata, FigMCP3Strata, FigMCP4Strata, FigMCP5Strata, labels = "AUTO", ncol = 1)
dev.off()


loo::loo_compare(fit_2st, fit_3st, fit_4st)
