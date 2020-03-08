#! /urs/bin/env Rscript
## ---------------------------
##
## Script name: 
##
## Purpose of script:计算基于otu表的nmds
# 似乎需要加上anosim分析
##
## Author: liuc
##
## Date Created: 2020-02-05
##
## Copyright (c) 
## Email:
##
## ---------------------------
##
## Notes:
## Usage: Rscript NMDS.R m2.txt map.txt outpath   
##
## ---------------------------

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
infile <- args[1]
map_file <- args[2]
outpath <- args[3]

library(tidyverse)
#library(ape)

if (!require(vegan)){
  install.packages('vegan')
} else {
  library(vegan)
}

options(scipen = 200)
##############################
# infile <- '/Users/congliu/prs/ruijin_xirou/combined_otu_table_m2_std.txt'
# map_file <- '/Users/congliu/prs/ruijin_xirou/map.txt'
# outpath <- '/Users/congliu/prs/ruijin_xirou/'
one_group <- c('polyp_stool', 'unpolyp_stool')
group1 <- one_group[1]
group2 <- one_group[2]

otu <- read.delim(infile, row.names = 1, sep = '\t', 
                  stringsAsFactors = FALSE, check.names = FALSE)
otu <- subset(otu, select = -taxonomy)
otu$sum <- rowSums(otu)
otu <- otu[order(otu$sum, decreasing = TRUE), ]
top10_otu <- rownames(otu[1:10, ])
otu <- subset(otu, select = -sum)
otu <- data.frame(t(otu))
#otu <- otu[!row.names(otu) %in% c('taxonomy'), ]
mapfile <- read_delim(map_file, delim = '\t') %>% select(Type, Description)
#mapfile <- subset(read.delim('~/prs/ruijin_pang/map.txt', sep = '\t'), select = c(Type, Description))
###anosim分析
otu2 <- otu
otu2$Description <- row.names(otu2)
otu3 <- merge(otu2, mapfile, by.x = 'Description', by.y = 'Description')
row.names(otu3) <- otu3$Description
#otu3 <- subset(otu3, Type %in% one_group)
otu4 <- subset(otu3, select = -c(Description, Type))


# anosim分析 ----------------------------------------------------------------

anosim_result_otu <- anosim(otu4, otu3$Type, permutations = 999, distance = 'bray')
(anosim_all_plot <- plot(anosim_result_otu, col=c('red4', 'blue3', 'yellow')))
ggsave(anosim_all_plot, filename = file.path(outpath, 'anosim_all_boxplot.png'), 
       height = 8, width = 8, dpi = 800)
Rvalue <- anosim_result_otu$statistic
Pvalue <- anosim_result_otu$signif

group_name <- unique(otu3$Type) #或者是指定的分组
anosim_result_two <- NULL
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij <- subset(otu3, Type %in% c(group_name[i], group_name[j]))
    otu_ij <- otu4[group_ij$Description, ]
    anosim_result_otu_ij <- anosim(otu_ij, group_ij$Type, permutations = 999, distance = 'bray')
    if (anosim_result_otu_ij$signif <= 0.001) Sig <- '***'
    else if (anosim_result_otu_ij$signif <= 0.01) Sig <- '**'
    else if (anosim_result_otu_ij$signif <= 0.05) Sig <- '*'
    else Sig <- NA
    anosim_result_two <- rbind(anosim_result_two, c(paste(group_name[i], group_name[j], sep = '/'), 
                                                    'Bray-Curtis', anosim_result_otu_ij$statistic, 
                                                    anosim_result_otu_ij$signif, Sig))
    png(file.path(outpath, paste('anosim_', group_name[i], '_', group_name[j], '.png', sep = '')), 
        width = 800, height = 800)
    plot(anosim_result_otu_ij, col = c('red4', 'blue3', 'yellow'))
    dev.off()
  }
}
anosim_result_two <- data.frame(anosim_result_two, stringsAsFactors = FALSE)
names(anosim_result_two) <- c('group', 'distance', 'R', 'P_value', 'Sig')
write.table(format(anosim_result_two, scientific = FALSE), file = file.path(outpath, 'ANOSIM.result_two.txt'), row.names = FALSE, 
            sep = '\t', quote = FALSE, na = '')
#PERMANOVA分析
adonis_result_otu <- adonis(otu4~Type, otu3, permutations = 999, distance = 'bray')
otuput <- data.frame(adonis_result_otu$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
write.table(otuput, file = file.path(outpath, 'PERMANOVA.result_all.txt'), row.names = FALSE, 
            sep = '\t', quote = FALSE, na = '')
adonis_result_two <- NULL
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij <- subset(otu3, Type %in% c(group_name[i], group_name[j]))
    otu_ij <- otu4[group_ij$Description, ]
    adonis_result_otu_ij <- adonis(otu_ij~Type, group_ij, permutations = 999, distance = 'bray')
    adonis_result_two <- rbind(adonis_result_two, c(paste(group_name[i], group_name[j], sep = '/'), 
                                                    'Bray-Curtis', 
                                                    unlist(data.frame(adonis_result_otu_ij$aov.tab, check.names = FALSE)[1, ])))
  }
}
adonis_result_two <- data.frame(adonis_result_two, stringsAsFactors = FALSE)
names(adonis_result_two) <- c('group', 'distance', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
for (i in 1:nrow(adonis_result_two)) {
  if (adonis_result_two[i, 'Pr (>F)'] <= 0.001) adonis_result_two[i, 'Sig'] <- '***'
  else if (adonis_result_two[i, 'Pr (>F)'] <= 0.01) adonis_result_two[i, 'Sig'] <- '**'
  else if (adonis_result_two[i, 'Pr (>F)'] <= 0.05) adonis_result_two[i, 'Sig'] <- '*'
}
write.table(adonis_result_two, file = file.path(outpath, 'PERMANOVA.result_two.txt'), 
            row.names = FALSE, sep = '\t', quote = FALSE, na = '')


####NMDS分析
Color <- c('magenta2', 'blue3', 'green', 'darkorange', 'yellow', 
           'hotpink', 'grey', 'skyblue', 'violet', 'cadetblue', '#0072B2')
#排序，预设 n个排序轴
nmds1 <- metaMDS(otu, distance = 'bray', k = 2)
#提取应力函数值（stress）
nmds1.stress <- nmds1$stress
#提取样本排序坐标
nmds1.point <- data.frame(nmds1$point)
#提取物种（OTU）排序坐标
nmds1.species <- data.frame(nmds1$species)
#############全部样本的NMDS分析################
plot(nmds1) #样本坐标+OTU坐标
stressplot(nmds1)
nmds1.point$Description = row.names(nmds1.point)
input_data <- left_join(nmds1.point, mapfile, by = c('Description' = 'Description'))
names(input_data)[1:2] <- c('NMDS1', 'NMDS2')
species_site <-{nmds1.species[top10_otu, ]}[1:2]
species_site$group <- rownames(species_site)
names(species_site)[1:2] <- c('NMDS1', 'NMDS2')

sample_point <- ggplot(input_data, aes(x=NMDS1, y=NMDS2, colour=Type)) +
  geom_point(aes(shape=Type)) +
  stat_ellipse(linetype = 2, type = 't') +
  theme_bw() +
  #scale_shape_manual(values = c(16, 17)) +
  scale_color_manual(values = Color[1:length(group_name)]) +
  labs(title = paste('Stress =', round(nmds1$stress, 3))) +
  geom_text(aes(label = group), data = species_site, color = 'green4', size = 2) + 
  theme(plot.title = element_text(hjust = 0.5)
  )

ggsave(sample_point, filename = file.path(outpath, paste('NMDS_all', '.pdf', sep = '')), height = 8, width = 8)


##################分组，应该怎么分组###################
for (i in 1:(length(group_name)-1)) {
  for(j in (i+1):length(group_name)){
    input_data_ij <- subset(input_data, Type %in% c(group_name[i], group_name[j]))
    sample_point_ij <- ggplot(input_data_ij, aes(x=NMDS1, y=NMDS2, colour=Type)) +
      geom_point(aes(shape=Type)) +
      stat_ellipse(linetype = 2, type = 't') +
      theme_bw() +
      scale_shape_manual(values = c(16, 17)) +
      scale_color_manual(values = c('magenta2', 'blue3')) +
      labs(title = paste('Stress =', round(nmds1$stress, 3))) +
      geom_text(aes(label = group), data = species_site, color = 'green4', size = 2) + 
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(sample_point_ij, 
           filename = file.path(outpath, 
                                paste('NMDS', group_name[i], group_name[j], '.png',sep = '_')),
           height = 8, width = 8, dpi = 800)
  }
}
#######################分组样本的NMDS分析##########################
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij <- subset(otu3, Type %in% c(group_name[i], group_name[j]))
    otu_ij <- otu4[group_ij$Description, ]
    otu_top <- data.frame(t(otu_ij))
    otu_top$sum <- rowSums(otu_top)
    otu_top <- otu_top[order(otu_top$sum, decreasing = TRUE), ]
    top10_otu_ij <- rownames(otu_top[1:10, ])
    nmds_ij <- metaMDS(otu_ij, distance = 'bray', k = 2)
    nmds_ij_stress <- nmds_ij$stress
    nmds_ij_point <- data.frame(nmds_ij$point)
    nmds_ij_species <- data.frame(nmds_ij$species)
    nmds_ij_point$Description = row.names(nmds_ij_point)
    input_data_ij <- left_join(nmds_ij_point, mapfile, by = c('Description' = 'Description'))
    names(input_data_ij)[1:2] <- c('NMDS1', 'NMDS2')
    species_site_ij <- {nmds_ij_species[top10_otu_ij, ]}[1:2]
    species_site_ij$group <- rownames(species_site_ij)
    names(species_site_ij)[1:2] <- c('NMDS1', 'NMDS2')
    sample_point_ij <- ggplot(input_data_ij, aes(x=NMDS1, y=NMDS2, colour=Type)) +
      geom_point(aes(shape=Type)) +
      stat_ellipse(linetype = 2, type = 't') +
      theme_bw() +
      #scale_shape_manual(values = c(16, 17)) +
      scale_color_manual(values = Color[1:length(group_name)]) +
      labs(title = paste('Stress =', round(nmds_ij_stress, 3))) +
      geom_text(aes(label = group), data = species_site_ij, color = 'green4', size = 2) + 
      theme(plot.title = element_text(hjust = 0.5)
      )
    ggsave(sample_point_ij, filename = file.path(outpath, 
                                                 paste('NMDS', group_name[i], group_name[j], '.png',sep = '_')), 
           height = 8, width = 8)
  }
}
