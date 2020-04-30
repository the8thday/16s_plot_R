#! /urs/bin/env Rscript
##旭日图
library(tidyverse)
#library(plotly)
library(sunburstR)
#library(htmltools)
#library(d3r)
library(ggrepel)

#一个简单的旭日图，2 层的圆环状结构，参考自
#https://stackoverflow.com/questions/26748069/ggplot2-pie-and-donut-chart-on-same-plot
# https://github.com/marbl/Krona/wiki
#https://stackoverflow.com/questions/50004058/multiple-dependent-level-sunburst-doughnut-chart-using-ggplot2

# data tidy ---------------------------------------------------------------

infile <- '/Users/congliu/prs/ruijin_xirou/combined_otu_table_m2_std.txt'


otu <- read_delim(infile, delim = '\t')
otu_one <- otu %>% dplyr::select(c(PRS029180041, taxonomy)) %>% 
  arrange(desc(PRS029180041)) %>% filter(taxonomy != 'Unassigned') %>% 
  slice(1:50) %>% separate(col = taxonomy, 
                           into = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'),
                           sep = "; ") %>% pivot_longer(-PRS029180041)

otu_one$name <- factor(otu_one$name, levels = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'))

ggplot(otu_one, aes(x = name, y = PRS029180041, fill = value, alpha = name)) +
  
  geom_col(width = 1, color = 'gray90', size = 0.25, position = position_stack()) +
  geom_text(aes(label = value), size = 2.5, position = position_stack(vjust = 0.5), stat = 'identity', check_overlap = T) +
  #geom_text(aes(label = value), size = 2.5, stat = 'identity',check_overlap = TRUE, angle= 90, hjust = -.2) +
  coord_polar(theta = 'y') +
  scale_alpha_manual(values = c('kingdom' = 1, 'phylum' = 1, 'class' = 1, 'order'=1, 
                                'family'=1, 'genus'=1, 'specise'=1), guide = F) +
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  #scale_fill_brewer(palette = 'Dark2', na.translate = F) +
  labs(x = NULL, y = NULL) +
  theme_minimal() + guides(fill=FALSE)


# sunburstR ----------------------------------------------------------------



