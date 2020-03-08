#! /urs/bin/env Rscript
# Rank Abundance

library(BiodiversityR)
library(vegan)

infile <- '/Users/congliu/prs/ruijin_xirou/combined_otu_table_m2_std.txt'
map_file <- '/Users/congliu/prs/ruijin_xirou/map.txt'
outpath <- '/Users/congliu/prs/ruijin_xirou/'

otu <- read.delim(infile, row.names = 1, sep = '\t', 
                  stringsAsFactors = FALSE, check.names = FALSE)
mapfile <- read_delim(map_file, delim = '\t') %>% select(Type, Description)
otu <- subset(otu, select = -taxonomy)

otu_relative <- otu / rowSums(otu) #转化为相对丰度
rank_dat <- data.frame()
for (i in rownames(otu_relative)) {
  rank_dat_i <- data.frame(rankabundance(subset(otu_relative, rownames(otu_relative) == i), digits = 6))[1:2]
  rank_dat_i$sample <- i
  rank_dat <- rbind(rank_dat, rank_dat_i)
}
rank_dat <- subset(rank_dat, abundance != 0)

#ggplot2 作图，更好的可视化效果，还请自行修改 ggplot2 作图细节
library(ggplot2)

ggplot(rank_dat, aes(rank, log(abundance, 10), color = sample)) +
  geom_line() +
  scale_colour_manual(limits = c('a1', 'a2', 'a3', 'a4', 'a5', 'a6'), 
                      values = c('orange', 'purple', 'green', 'blue', 'red', 'gray40')) +
  labs(x = 'OTUs rank', y = 'Relative adundance (%)', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.key = element_rect(fill = 'transparent')) +
  scale_y_continuous(breaks = 0:-5, labels = c('100', '10', '1', '0.1', '0.01', '0.001'), limits = c(-5, 0))


