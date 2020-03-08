#! /urs/bin/env Rscript
##manhadun

library(qqman)
library(tidyverse)


infile <- '/Users/congliu/prs/ruijin_xirou/lefse_input_polyp_stool_unpolyp_stool.res'
map_file <- '/Users/congliu/prs/ruijin_xirou/map.txt'
outpath <- '/Users/congliu/prs/ruijin_xirou/'
one_group <- c('polyp_stool', 'unpolyp_stool')
group1 <- one_group[1]
group2 <- one_group[2]

res <- read.delim(infile,  sep = '\t', 
                  stringsAsFactors = FALSE, check.names = FALSE, header = F)
mapfile <- read.delim(map_file, sep = '\t') %>% dplyr::select(Type, Description)
res <- subset(res, V4>0)
diff <- read.delim('/Users/congliu/prs/ruijin_xirou/polyp_stool_vs_unpolyp_stool_meta_diff.txt', sep = '\t',
                   header = T)
diff <- subset(diff, select=c(OTU, oddsRatio, lower, upper, pvalues, adjPvalues, taxonomy))
diff[which(diff$pvalues < 0.05), 'sig'] <- 'sign'
diff[which(diff$pvalues >= 0.05), 'sig'] <- 'no-sign'
###按照目排序
diff_o <- tidyr::separate(diff, taxonomy, into=c('a', 'b'), sep='o__') %>% dplyr::select(-a) %>%
  tidyr::separate(b, into=c('order', 'y'), sep='; f__') %>% dplyr::select(-y)
otu_stat <- diff_o[order(diff_o$order), ]
otu_stat$otu_sort <- 1:nrow(otu_stat)

p <- ggplot(otu_stat, aes(otu_sort, -log(pvalues, 10))) +
  geom_point(aes(size = pvalues, color = order, shape = sig)) +
  scale_size(range = c(1, 2))+
  scale_shape_manual(limits = c('sign', 'no-sign'), values = c(16, 1)) +     
  labs(x = NULL, y = '-log10(P)', size = 'relative abundance (%)', shape = 'enriched') +
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_rect(fill = 'transparent'), legend.key = element_rect(fill = 'transparent')) +
  guides(color = 'none') +
  geom_hline(yintercept = -log10(0.01), color = 'gray', linetype = 2, size = 1)
p
phylum_num <- length(otu_stat$order)
phylum_range <- c(0, as.numeric(phylum_num[[1]]))
phylum_name <- as.numeric(phylum_num[[1]]) / 2
for (i in 2:length(phylum_num)) {
  phylum_range[i+1] <- phylum_range[i] + phylum_num[i]
  phylum_name[i] <- phylum_range[i] + phylum_num[i] / 2
}
p <- p +
  scale_x_continuous(breaks = phylum_name, labels = names(phylum_num), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
