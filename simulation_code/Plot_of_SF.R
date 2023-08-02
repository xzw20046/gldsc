library(ggplot2)
library("gridExtra")
library("cowplot")
library('ggpubr')
library(scales)
library(svglite)
location='/Users/zewei/Desktop/HKUgrad/G-LDSC/plots/figure/'
ggsave(paste0(location,"figure2.svg"),plot = main.plot.2, width = 20, height = 14, units = "cm")
ggsave(paste0(location,"figure3.svg"),plot = main.plot.6, width = 14, height = 14, units = "cm")
ggsave(paste0(location,"figure4.svg"),plot = main.plot.4, width = 20, height = 16, units = "cm")
ggsave(paste0(location,"figure5.svg"),plot = main.plot.5, width = 20, height = 14, units = "cm")
ggsave(paste0(location,"figure1.pdf"),plot = main.plot.1, width = 114, height = 185, units = "mm")
ggsave(paste0(location,"figure2.pdf"),plot = main.plot.2, width = 174, height = 185, units = "mm")
ggsave(paste0(location,"figure3.pdf"),plot = main.plot.6, width = 114, height = 114, units = "mm")
ggsave(paste0(location,"figure4.pdf"),plot = main.plot.4, width = 174, height = 185, units = "mm")
ggsave(paste0(location,"figure5.pdf"),plot = main.plot.5, width = 174, height = 150, units = "mm")
#plot Figure 3 
source('~/Desktop/HKUgrad/G-LDSC/gldsc/gldsc/simulation_code/plot_of_N_mis.R')

ggsave(paste0(location,"figure3.pdf"),plot = boxplot.mis1, width = 174, height = 185, units = "mm")
#ggsave(paste0(location,"figureS4.svg"),plot = boxplot.mis2, width = 174, height = 185, units = "mm")

#Figure S1 S2 S3
source('~/Desktop/HKUgrad/G-LDSC/gldsc/gldsc/simulation_code/plot of N.R')
ggsave(paste0(location,"figureS1.svg"),plot = plot.h2.N, width = 174, height = 185, units = "mm")
ggsave(paste0(location,"figureS2.svg"),plot = plot.h2.N.004, width = 174, height = 185, units = "mm")
ggsave(paste0(location,"figureS3.svg"),plot = plot.e.N.004, width = 174, height = 185, units = "mm")

#Figure S4 S5 S6 S7 S8
source('~/Desktop/HKUgrad/G-LDSC/gldsc/gldsc/simulation_code/plot_confound_v2.R')
ggsave(paste0(location,"figureS4.svg"),plot = plot.h2.confound, width = 174, height = 185, units = "mm")
ggsave(paste0(location,"figureS5.svg"),plot = plot.e.confound, width = 174, height = 185, units = "mm")
ggsave(paste0(location,"figureS6.svg"),plot = plot.bias.1, width = 174, height = 185, units = "mm")
ggsave(paste0(location,"figureS7.svg"),plot = plot.bias.1.05, width = 174, height = 185, units = "mm")
ggsave(paste0(location,"figureS8.svg"),plot = plot.bias.1.1, width = 174, height = 185, units = "mm")

#Figure S9 S10
source('/Users/zewei/Desktop/HKUgrad/G-LDSC/plot codes/UKBB plots.v2.R')
ggsave(paste0(location,"figureS9.svg"),plot = plot.applica.comp.h2, width = 14, height = 14, units = "cm")
ggsave(paste0(location,"figureS10.svg"),plot = plot.average.e.rest, width = 20, height = 16, units = "cm")



