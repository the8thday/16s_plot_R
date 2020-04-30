#! /urs/bin/env Rscript
## ---------------------------
##
## Script name: PCA_PCoA_PLSDA_UPGMA.R
##
## Purpose of script:计算基于otu表的PCA\PCoA\PLSDA\CA\UPGMA
##
## Author: liuc
##
## Date Created: 2020-02-21
##
## Copyright (c) 
## Email:
##
## ---------------------------
##
## Notes: 每次只传入一组数据，所以要对map表进行筛选
## Usage: Rscript   
##
## ---------------------------

options(BioC_mirror='http://mirrors.ustc.edu.cn/bioc/')
#options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
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
# infile <- '/Users/congliu/prs/ZJDXGW/combined_otu_table_m2_std.txt'
# map_file <- '/Users/congliu/prs/ZJDXGW/map.txt'
# outpath <- '/Users/congliu/prs/ZJDXGW/'
# one_group <- c('A', 'B', 'C', 'D','E')
group1 <- one_group[1]
group2 <- one_group[2]
# group1 <- 'all'
# group2 <- 'group'


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
otu4 <- subset(otu3, select = -c(Description, Type)) #为所有计算所用的otu数值

#############PCoA#################
# I have used the weighted UniFrac metric to determine the distance between samples and PCoA to visualise the data
# 但是下文不是
distance <- vegan::vegdist(otu4, method = 'bray')
pcoa <- cmdscale(distance, k = (nrow(otu4) - 1), eig = TRUE)
#write.csv(as.matrix(distance), 'distance.csv', quote = F)
#ordiplot(scores(pcoa)[ ,c(1, 2)], type = 't')
species <- wascores(pcoa$points[,1:2], otu4)
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
ggsave(p_pcoa, filename = file.path(outpath, paste('pcoa_', group1, '_', group2, '.pdf', sep = '')), height = 8, width = 8)

#####################PCA####################
#Transformation-based Principal Component Analysis，tb-PCA
#因为PCA计算的欧几里得距离的原因, 而物种数据0较多，可选择进行hellinger，hellinger转化本身会降低高丰度物种的权重
otu5 <- otu
otu5$Description <- row.names(otu5)
otu6 <- merge(otu5, mapfile, by.x = 'Description', by.y = 'Description')
row.names(otu6) <- otu6$Description
otu6 <- subset(otu6, Type %in% one_group)
otu7 <- subset(otu6, select = -c(Description, Type))

otu_hel <- decostand(otu7, method = 'hellinger')
# rda scale为TRUE时排序对象是相关矩阵、scale为FASLE时排序对象是离差矩阵
pca_sp <- rda(otu_hel, scale = FALSE)
#ordiplot(pca_sp, scaling = 1, display = 'site', type = 'points')
pca_eig <- pca_sp$CA$eig
pca_exp <- pca_sp$CA$eig / sum(pca_sp$CA$eig)
#site.scaling1 <- scores(pca_exp, choices = 1:2, scaling = 1, display = 'site')
#env.scaling2 <- scores(pca_exp, choices = 1:2, scaling = 2, display = 'sp')
site.scaling1 <- data.frame(summary(pca_sp, scaling = 1)$sites[ ,1:2])
otu.scaling1 <- data.frame(summary(pca_sp, scaling = 1)$species[ ,1:2])
site.scaling1$sample <- rownames(site.scaling1)
site.scaling1 <- merge(site.scaling1, mapfile, by.x = 'sample', by.y = 'Description')
otu.scaling1$group <- rownames(otu.scaling1)
otu.scaling1 <- otu.scaling1[1:5,]
pc1 <- paste('PC1:', round(pca_exp[1]*100, 2), '%')
pc2 <- paste('PC2:',round(pca_exp[2]*100, 2), '%')

p_pca <- ggplot(site.scaling1, aes(PC1, PC2)) +
  geom_point(aes(color = Type)) +
  stat_ellipse(aes(PC1, PC2, color = Type), show.legend = FALSE, linetype = 2) +
  scale_color_manual(values = c('red4', 'green3', '#C673FF2E', 'orange', 'blue2')) +
  scale_shape_manual(values = c(16, 17)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title = element_blank(), legend.key = element_rect(fill = 'transparent')) + 
  labs(x = pc1, y = pc2) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_segment(data = otu.scaling1, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, 'cm')), size = 0.3, color = 'blue') +
  geom_text(data = otu.scaling1, aes(PC1 * 1.1, PC2 * 1.1, label = group), color = 'blue', size = 3)

ggsave(p_pca, filename = file.path(outpath, paste('pca_', group1, '_', group2, '.pdf', sep = '')), 
       height = 8, width = 8)

##################PLSDA######################
# otu.plsda <- opls(x = otu4, y = factor(otu3$Type), orthoI = 0)
# head(otu.plsda@scoreMN) 
# plot(otu.plsda, typeVc = 'overview', parAsColFcVn = factor(otu3$Type))

plsda_result <-plsda(otu4, otu3$Type, ncomp = 2)
#plotIndiv(plsda_result, ind.names = TRUE, style = 'ggplot2')
sample_site <- data.frame(plsda_result$variates)[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('plsda1', 'plsda2')
sample_site <- merge(sample_site, mapfile, by.x = 'names', by.y = 'Description', all.x = TRUE)
plsda_result_eig <- {plsda_result$explained_variance$X}[1:2]

p_plsda <- ggplot(sample_site, aes(plsda1, plsda2, color = Type, label = names)) +
  geom_point(size = 1.5, alpha = 0.6) + 
  stat_ellipse(show.legend = FALSE, linetype = 2) +
  scale_color_manual(values = c('#1D7ACC', '#F67433', '#00815F','#C673FF2E', 'blue2')) +
  scale_shape_manual(values = c(16, 17)) +
  theme(panel.grid = element_line(color = 'grey50'), 
        panel.background = element_rect(color = 'black', fill = 'transparent')) + 
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent')) +
  labs(x = paste('PLS-DA axis1 ( explained variance ', round(100 * plsda_result_eig[1], 2), '% )', sep = ''), 
       y = paste('PLS-DA axis2 ( explained variance ', round(100 * plsda_result_eig[2], 2), '% )', sep = ''))
ggsave(p_plsda, filename = file.path(outpath, paste('plsda_', group1, '_', group2, '.pdf', sep = '')), height = 8, width = 8)
# 
###################CA/DCA/CCA/RDA####################
# CA对应分析，展示的为样方之间的卡方距离，卡方距离不受零值的影响，可以直接输入物种矩阵数据
# 选择用RDA还是CCA分析？先用“样本-物种”文件做DCA分析！DCA分析后大于4一般选择单峰模型，否则选择线性模型
decorana(sp)
otu5 <- otu
otu5$Description <- row.names(otu5)
otu6 <- merge(otu5, mapfile, by.x = 'Description', by.y = 'Description')
row.names(otu6) <- otu6$Description
otu6 <- subset(otu6, Type %in% one_group) ##挑选出一组样本，提供分组信息
otu7 <- subset(otu6, select = -c(Description, Type)) #为所有计算所用的otu数值
ca <- vegan::cca(otu7)
ca_eig <- ca$CA$eig
ca_exp <- ca$CA$eig / sum(ca$CA$eig)
ca.scaling1 <- data.frame(summary(ca, scaling = 1)$sites[ ,1:2])
otu.scaling1 <- data.frame(summary(ca, scaling = 1)$species[ ,1:2])
ca.scaling1$sample <- rownames(ca.scaling1)
ca.scaling1 <- merge(ca.scaling1, mapfile, by.x = 'sample', by.y = 'Description', all.x = TRUE)
otu.scaling1$group <- rownames(otu.scaling1)
#只保留 top10 丰度的细菌门
top10_otu <- names(colSums(otu7)[order(colSums(otu7), decreasing = TRUE)][1:10])
top10_otu.scaling1 <- otu.scaling1[top10_otu, ]

p_ca <- ggplot(ca.scaling1, aes(CA1, CA2)) +
  geom_point(aes(color = Type)) +
  scale_color_manual(values = c('red', 'blue2')) +
  scale_shape_manual(values = c(16, 17)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'), plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_point(data = top10_otu.scaling1, color = 'green') +
  geom_text_repel(data = top10_otu.scaling1, aes(label = group), 
                 size = 3, box.padding = unit(0.5, 'lines'), color = 'blue')
  #labs(x = ca1, y = ca2, title = 'I 型标尺，仅展示 top10 丰度种群', color = NULL)

ggsave(p_ca, filename = file.path(outpath, paste('ca_', group1, '_',  group2, '.pdf', sep = '')), height = 8, width = 8)

##################UPGMA##################
distance <- vegan::vegdist(otu4, method = 'bray')
  
dis_ch <- vegdist(decostand(otu4, 'normalize'), method = 'euc') #弦距离，权重的？
#使用中值的非权重成对组法，UPGMA
clust_average <- hclust(distance, method = 'average')
#使用质心的非权重成对组法，UPGMC
clust_centroid <- hclust(distance, method = 'centroid')
#使用中值的权重成对组法，WPGMA
clust_mcquitty <- hclust(distance, method = 'mcquitty')
#使用质心的权重成对组法，WPGMC
clust_median <- hclust(distance, method = 'median')
##Ward 最小方差聚类（包含两种不同的 Ward 聚类算法，详情 ?hclust）
#该方法中只能使用欧式几何属性的距离测度，因此更换为平方根转化后的 Bray-curtis 距离测度
dis_bray_euc <- distance^0.5
#ward.D
clust_ward <- hclust(dis_bray_euc, method = 'ward.D')
#ward.D2
clust_ward2 <- hclust(dis_bray_euc, method = 'ward.D2')

clust_average <- hclust(distance, method = 'average')
#plot(clust_average, main = 'UPGMA\nBray-curtis', sub = '', xlab = 'Sample', ylab = 'Height')
#rect.hclust(clust_average, k = 2, border = c('red', 'blue'))
hcd=as.dendrogram(clust_average)
labelColors=brewer.pal(n=4, name="Set1")
clusMember <- otu3$Type
names(clusMember) <- otu3$Description
#clusMember=cutree(hcd, 4)
hcd=hcd %>% set("labels_cex", 1.5) %>% set("branches_lwd", 2) %>% 
  set("branches_k_color", k=4) %>% set("branches_k_lty", k=4)
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, 'nodePar') <- c(a$nodePar, lab.col = labCol)
  }
  n
}
clusDendro=dendrapply(hcd, colLab)
# pdf(file.path(outpath, paste('upgma_', group1, '_', group2, '.pdf', sep = '')), 
#     width = 800, height = 800)
#plot(clusDendro, main="UPGMA Tree", type="rectangle", horiz=FALSE, ylab = 'Height', cex=0.1)
plot(clusDendro, main = 'UPGMA\n(Bray-curtis distance)', sub = '', ylab = 'Height')
#axis(side=1, cex.lab=0.1)
dev.off()

















