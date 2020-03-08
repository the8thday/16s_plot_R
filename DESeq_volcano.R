#! /urs/bin/env Rscript
# DESeq2 火山图

#options(BioC_mirror='http://mirrors.ustc.edu.cn/bioc/')
library(DESeq2)
library(tidyverse)
library(patchwork)

infile <- '/Users/congliu/prs/ruijin_xirou/combined_otu_table_m2_std.txt'
map_file <- '/Users/congliu/prs/ruijin_xirou/map.txt'
outpath <- '/Users/congliu/prs/ruijin_xirou/'
one_group <- c('polyp_saliva', 'unpolyp_saliva')
group1 <- one_group[1]
group2 <- one_group[2]

otu <- read.delim(infile, row.names = 1, sep = '\t', 
                  stringsAsFactors = FALSE, check.names = FALSE)
mapfile <- read.delim(map_file, sep = '\t') %>% dplyr::select(Type, Description)
row.names(mapfile) <- mapfile$Description
group_file <- subset(mapfile, select = -c(Description))
otu <- subset(otu, select = -taxonomy)
otu <- data.frame(t(otu))
otu2 <- otu
otu2$Description <- row.names(otu2)
otu3 <- merge(otu2, mapfile, by.x = 'Description', by.y = 'Description')
row.names(otu3) <- otu3$Description
otu3 <- subset(otu3, Type %in% one_group) ##挑选出一组样本，提供分组信息
otu4 <- subset(otu3, select = -c(Description, Type))
otu5 <- t(otu4)
otu5[rowSums(otu5)>10, ]
group_file <- subset(otu3, select=c(Type))


dds <- DESeqDataSetFromMatrix(countData = otu5, colData = group_file, design = ~Type)
dds <- dds[rowSums(counts(dds)) >= 10,]
dds <- DESeq(dds)
suppressMessages(dds)
res <- results(dds, contrast=c('Type', group1, group2))
##结果提取
deseq_res <- as.data.frame(res[order(res$padj), ])
deseq_res$otu_id <- rownames(deseq_res)
#可选使用下命令将“baseMean”这一列转化为相对丰度数值
#deseq_res$baseMean <- deseq_res$baseMean / sum(deseq_res$baseMean)
#write.table(deseq_res[c(7, 1:6)], 'DESeq2.txt', row.names = FALSE, sep = '\t', quote = FALSE)
deseq_res <- deseq_res[complete.cases(deseq_res[,2]),]
for (i in 1:nrow(deseq_res)) {
  if (abs(deseq_res[i,'log2FoldChange']) >= 1) deseq_res[i,'select_change'] <- 'y' else deseq_res[i,'select_change'] <- 'n'
  if (deseq_res[i,'padj'] %in% NA | abs(deseq_res[i,'padj']) >= 0.05) deseq_res[i,'select_pvalue'] <- 'n' else deseq_res[i,'select_pvalue'] <- 'y'
  deseq_res[i,'select'] <- paste(deseq_res[i,'select_change'], deseq_res[i,'select_pvalue'], sep = '')
}

deseq_res$select <- factor(deseq_res$select, levels = c('nn', 'ny', 'yn', 'yy'), 
                           labels = c('p >= 0.05, FC < 2', 'p < 0.05, FC < 2', 'p >= 0.05, FC >= 2', 'p < 0.05, FC >= 2'))
write.table(deseq_res, 'DESeq2_2.txt', row.names = FALSE, sep = '\t', quote = FALSE)
#纵轴为显著性 p 值
volcano_plot_pvalue <- ggplot(deseq_res, aes(log2FoldChange, -log(padj, 10))) +
  geom_point(aes(color = select), alpha = 0.6) +
  scale_color_manual(values = c('gray30', 'green4', 'red2', 'blue2')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-0.5, 0.5), color = 'gray', size = 0.5) + 
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.5) +
  labs(x = 'log2 Fold Change', y = '-log10 p-value')

volcano_plot_abundance <- ggplot(deseq_res, aes(log2FoldChange, 100 * baseMean / sum(deseq_res$baseMean))) +
  geom_point(aes(color = select), alpha = 0.6) +
  scale_color_manual(values = c('gray30', 'green4', 'red2', 'blue2')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-0.5, 0.5), color = 'gray', size = 0.5) + 
  labs(x = 'log2 Fold Change', y = 'Abundance (%)')
volcano_plot_abundance / volcano_plot_pvalue
