#! /urs/bin/env Rscript
# Y轴截断

# 加载包 ---------------------------------------------------------------------

library(tidyverse)
library(plotrix)
library(patchwork)

# tidy data ---------------------------------------------------------------

infile <- '/Users/congliu/prs/ZJDXGW/combined_otu_table_m2_std.txt'
map_file <- '/Users/congliu/prs/ZJDXGW/map.txt'
outpath <- '/Users/congliu/prs/ZJDXGW/'

otu <- read.delim(infile, row.names = 1, sep = '\t', 
                  stringsAsFactors = FALSE, check.names = FALSE)
mapfile <- read.delim(map_file, sep = '\t') %>% dplyr::select(Type, Description)
otu <- subset(otu, taxonomy != 'Unassigned')
otu <- subset(otu, select = -taxonomy)
otu <- data.frame(t(otu))
#otu2 <- log(otu + 1)
otu2 <- otu
otu2$Description <- row.names(otu2)
otu3 <- merge(otu2, mapfile, by.x = 'Description', by.y = 'Description')
row.names(otu3) <- otu3$Description
otu3 <- subset(otu3, Type %in% one_group) ##挑选出一组样本，提供分组信息
otu4 <- subset(otu3, select = -c(Description, Type))
# plot --------------------------------------------------------------------

ab <- read.delim('/Users/congliu/prs/ZJDXGW/A_vs_D_deseq2_diff.txt', sep = '\t', 
                 stringsAsFactors = FALSE, check.names = FALSE)
ab <- subset(ab, padj < 0.05, select= -taxonomy)[1:20, ]
ab['group'] <- ifelse(ab$log2FoldChange > 0, 'a', 'n')
filter_ <- c(as.vector(ab$OTU),c('Type'))
final_otu <- otu3[filter_] %>% 
  group_by(Type) %>% summarise_all(funs(mean(., na.rm=TRUE), sd)) %>% 
  pivot_longer(-Type, names_to = c('OTU', 'fun'), names_sep = '_') %>% 
  pivot_wider(names_from = fun, values_from = value)
#%>% select(-Description)
#aggregate(final_otu[, 1:(length(final_otu)-1)], list(final_otu$Type), mean)

group <- unique(final_otu$Type)

# gap.barplot(ab$value.mean, gap=c(0,4,20), col=ab$taxonomy, 
#             xlab="index", ylab="value",
#             horiz = F)
# axis.break(2, from, breakcol="snow", style="gap")
# axis.break(2, from*(1+0.02), breakcol="black", style="slash")
# axis.break(4, from*(1+0.02), breakcol="black", style="slash")
# axis(2, at=from)

# plot by ggplot2 ---------------------------------------------------------

bar_plot <- ggplot(final_otu, aes(OTU, mean, fill = Type)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
               width = 0.5, position=position_dodge(0.8)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())

split1 <- bar_plot + coord_cartesian(ylim = c(0, 1000)) + 
  theme(legend.position='none')
split2 <- bar_plot + coord_cartesian(ylim = c(20000, 40000)) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none')
split3 <- bar_plot + coord_cartesian(ylim = c(50000, 200000)) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        #legend.position = c(0.95, 0.7)
        )

split3 / split2/ split1









