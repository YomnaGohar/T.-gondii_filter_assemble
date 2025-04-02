#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
Genome_size <- as.integer(args[3])
# Load necessary libraries
library(ggplot2)
library(scales)

dat <- read.table(file =input_file, header = FALSE,sep="",col.names = c('Chr','Position','Coverage'))
# Plot bar-plot
unique.levels <- sort(unique(dat$Coverage))
count <- table(dat$Coverage)
count.df <- data.frame(unique.levels, count/Genome_size*100)
p=ggplot(count.df, aes(factor(unique.levels),Freq, color=factor(unique.levels), fill =factor(unique.levels)))+
  geom_bar(stat="identity") + 
  labs(title="Per-base coverage of T.gondii genome by the assembly",
       subtitle="",
       y="percent %", x="coverage per base") + 
  theme(legend.position="none") +
  scale_y_continuous(labels = comma, breaks=seq(0,max(count.df$Freq), 10), limits = c(0,max(count.df$Freq))) #to suppress the scientfic notation on the y axis (ex,5e+00 to represent 5)
ggsave(output_file, plot = p, width = 10, height = 10)
