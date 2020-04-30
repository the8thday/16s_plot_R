#! /urs/bin/env Rscript
# 加上PERMANOVA分析


args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
infile <- args[1]
map_file <- args[2]
outpath <- args[3]
group_str <- args[4]
one_group <- eval(parse(text = (group_str)))
print(one_group)


#library(plyr)
library(tidyverse)
library(vegan)
library(ropls)
library(mixOmics)
library(dendextend)
library(RColorBrewer)
library(ggrepel)

tax_color <- c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', 
               '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', 'gray')
infile <- '/Users/congliu/prs/WXRM/Reporter/left_with_right/combined_otu_table_m2_std.txt'
map_file <- '/Users/congliu/prs/WXRM/Reporter/left_with_right/map.txt'
outpath <- '/Users/congliu/prs/WXRM/Reporter/left_with_right'
one_group <- c('left', 'right') #, 'Polyp', 'Carcinoma','Post')
group1 <- one_group[1]
group2 <- one_group[2]
# group1 <- 'all'
# group2 <- 'group'


input_data <- function(){
  otu <- read.delim(infile, row.names = 1, sep = '\t', 
                  stringsAsFactors = FALSE, check.names = FALSE)
  mapfile <- read.delim(map_file, sep = '\t') %>% dplyr::select(Type, Description)
  otu <- subset(otu, taxonomy != 'Unassigned')
  otu <- subset(otu, select = -taxonomy)
  otu <- data.frame(t(otu))
  otu2 <- log(otu + 1)
  #otu2 <- otu
  otu2$Description <- row.names(otu2)
  otu3 <- merge(otu2, mapfile, by.x = 'Description', by.y = 'Description')
  row.names(otu3) <- otu3$Description
  otu3 <- subset(otu3, Type %in% one_group) ##挑选出一组样本，提供分组信息
  otu4 <- subset(otu3, select = -c(Description, Type)) #为所有计算所用的otu数值
  return(otu4)
}
otu4 <- input_data()
#############PCoA#################
# I have used the weighted UniFrac metric to determine the distance between samples and PCoA to visualise the data
# 首先是常用的bray距离
pcoa_bray <- function(){
  distance <- vegan::vegdist(otu4, method = 'bray')
  pcoa <- cmdscale(distance, k = (nrow(otu4) - 1), eig = TRUE)
  #write.csv(as.matrix(distance), 'distance.csv', quote = F)
  #ordiplot(scores(pcoa)[ ,c(1, 2)], type = 't')
  species <- wascores(pcoa$points[,1:5], otu4)
  pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
  sample_site <- data.frame({pcoa$point})[1:2]
  sample_site$names <- rownames(sample_site)
  names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
  sample_site <- merge(sample_site, mapfile, by.x = 'names', by.y = 'Description', all.x = TRUE)
  group_border <- plyr::ddply(sample_site, 'Type', function(df) df[chull(df[[2]], df[[3]]), ])
  
  p_pcoa <- ggplot(sample_site, aes(PCoA1, PCoA2, group = Type)) +
    theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
          panel.background = element_rect(color = 'black', fill = 'transparent'), 
          legend.key = element_rect(fill = 'transparent')) +
    geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
    geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
    geom_polygon(data = group_border, aes(fill = Type)) +
    geom_point(aes(color = Type), size = 1.5, alpha = 0.8) +
    scale_shape_manual(values = c(17, 16)) +
    scale_color_manual(values = c('orange', 'blue2', 'red4', '#C673FF2E', 'green2')) +
    scale_fill_manual(values = c('#C673FF2E', '#73D5FF2E', '#49C35A2E', '#FF985C2E', '#C673FF2E')) + #可在这里修改区块的颜色
    guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2), color = guide_legend(order = 3)) +
    labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig[1], 2), '%'), 
         y = paste('PCoA axis2: ', round(100 * pcoa_eig[2], 2), '%'))
  #可通过修改下面四句中的点坐标、大小、颜色等，修改“A、B、C、D”标签
  # annotate('text', label = 'A', x = -0.31, y = -0.15, size = 5, colour = '#C673FF') +
  # annotate('text', label = 'B', x = -0.1, y = 0.3, size = 5, colour = '#73D5FF') +
  # annotate('text', label = 'C', x = 0.1, y = 0.15, size = 5, colour = '#49C35A') +
  # annotate('text', label = 'D', x = 0.35, y = 0, size = 5, colour = '#FF985C')
  ggsave(p_pcoa, filename = file.path(outpath, paste('pcoa_2_', group1, '_', group2, '.pdf', sep = '')), height = 8, width = 8)
}
pcoa_bray()

#断棍模型评估各轴特征值
pcoa_eig <- pcoa$eig
n <- length(pcoa_eig)
bsm <- data.frame(j=seq(1:n), p = 0)
bsm$p[1] <- 1/n
for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
bsm$p <- 100*bsm$p/n

barplot(t(cbind(100 * pcoa_eig/sum(pcoa_eig), bsm$p[n:1])), beside = TRUE, main = '% change', col = c('orange', 'bisque'), las = 2)
abline(h = mean(100 * pcoa_eig/sum(pcoa_eig)), col = 'red')
legend('topright', c('% 特征值', '断棍模型', '平均特征值'), pch = c(15, 15, NA), col = c('orange', 'bisque', 'red'), lwd = c(NA, NA, 1), bty = 'n')
# 由断棍模型可以看到用前两个轴作为代表显得勉强，而最后面的负特征根在末端的几个轴中，不会带来很大的影响

# 将物种信息加入pcoa中
ordiplot(pcoa, type = 'text', main = 'PCoA with Species')
points(species[ ,1:2], pch = 3, cex = 0.7, col = 'pink')
abundance <- apply(otu4, 2, sum)
abundance_top10 <- names(abundance[order(abundance, decreasing = TRUE)][1:10])


# weighted UniFrac --------------------------------------------------------

pcoa_unifrac <- function(){
  distance <- phyloseq::distance()
}
