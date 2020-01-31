#!/usr/bin/env Rscript

library(ggplot2)
# UASg: start, midpoint, end positions
UASg.df <- data.frame(x=c(278568, 278627, 278685), y=25)

folder <- '~/Desktop/NGS/Project_6686/Sample_CS1001-009_IGO_06686_9/'
file <- file.path(folder, "gal4-UASg-vplot.csv", fsep = .Platform$file.sep)
binding <- read.csv(file)
colnames(binding) <- c('midpoint', 'length')
ggplot(binding, aes(x=midpoint, y=length)) + geom_point(alpha = 0.5, color="blue", size=.1) +
    geom_line(data=UASg.df , aes(x, y)) +
    geom_point(data=UASg.df, aes(x=x,y=y))
