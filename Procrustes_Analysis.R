#! /urs/bin/env Rscript

library(tidyverse)
library(vegan)

#http://john-quensen.com/tutorials/procrustes-analysis/
# tidy data ---------------------------------------------------------------

env <- read.delim('/Users/congliu/OneDrive/普瑞森/200328-普鲁克分析（Procrustes Analysis）评估物种-环境或功能关联度的一个示例/env_table.txt',
                  sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

otu <- read.delim('/Users/congliu/OneDrive/普瑞森/200328-普鲁克分析（Procrustes Analysis）评估物种-环境或功能关联度的一个示例/spe_table.txt',
                  sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

group <- read.delim('/Users/congliu/OneDrive/普瑞森/200328-普鲁克分析（Procrustes Analysis）评估物种-环境或功能关联度的一个示例/group.txt',
                    sep = '\t', stringsAsFactors = FALSE)

# 在进行比较之前因数据属性不同，直接比较两数据集不合适，因此先做PCA处理(或其他降维的数据PCoA,NMDS)
env.std <- decostand(env, method = "standardize")
env_pca <- rda(env, scale = TRUE) #scale为TRUE,数据正态分布
otu_hell <- decostand(otu, method = 'hellinger')
otu_pca <- rda(otu_hell, scale = FALSE)
# biplot(env_pca, choices = c(1,2), scaling = 1, col = c('red2', 'blue2'))
# biplot(otu_pca, choices = c(1,2), scaling = 1, col = c('red2', 'blue2'))
# 从示例数据看二者间数据分布相似，在实际数据中要仔细分析数据集的分布

# Procrustes 分析 -----------------------------------------------------------

#??procrustes
# X:target matrix, Y: Matrix to be rotated, symmetric为FALSE时m^2会随XY转换而不同，symmetric为TRUE时m^2 is same
#提取两个 PCA 中的样方排序坐标，均以 I 型标尺为例；PCA维度不一致。
site_env <- summary(env_pca, scaling = 1)$site
site_otu <- summary(otu_pca, scaling = 1)$site
proc <- procrustes(X=site_env, Y=site_otu, symmetric = TRUE)
# proc$X：分析后的X坐标
# proc$Yrot: 分析后的Y坐标, 可以看到坐标名称发生了变化, 即是是按照X为基础 进行的转换
# proc$ss  偏差平方和 M2 统计量
# proc$rotation  通过该值可获得旋转轴的坐标位置
plot(proc, kind = 1)
plot(proc, kind = 2)
# residuals(proc) #残差值该图显示了每个配对样方的残差，从底部到顶部的水平线是残差的25％（虚线），对应点坐标之间的偏差称为矢量残差（vector residuals）
# 50％（实线）和75％（虚线）分位数。如果样方中环境和物种组成更为一致，则旋转图中两个匹配的点距离越近，
# 残差较小，反之环境和物种更无规律，则残差更大。

# M2统计量的显著性检验（PROTEST）
?protest()
#注：protest() 中执行的是对称 Procrustes 分析，X 和 Y 的分配调换不影响 M2 统计量的计算
set.seed(42)
prot <- protest(X = env_pca, Y = otu_pca, permutations = how(nperm = 999))
# prot$signif  #p 值
# prot$ss  #偏差平方和 M2 统计量
signif <- prot$signif
m2 <- prot$ss

# ggplot2 绘图 --------------------------------------------------------------
#提取 Procrustes 分析的坐标, 只绘制两个数据集中前两个维度
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)
Y$sample <- rownames(Y)
Y <- merge(Y, group, by.x = 'sample', by.y = 'samples')
p <- ggplot(data = Y) +
  geom_point(aes(x = X1, y = X2, colour = groups), size = 1.5, shape = 16) +
  geom_point(aes(x = PC1, y = PC2, color = groups), size = 1.5, shape = 2) +
  geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.05, 'cm'))) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(x = 'Dimension 1', y = 'Dimension 2') +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
  annotate('text', label = sprintf(paste('M^2 == ', signif)),
           x = -0.21, y = 0.42, size = 3, parse = TRUE) +
  annotate('text', label = paste('P < ', round(m2, 3)),
           x = -0.21, y = 0.38, size = 3, parse = TRUE)

p








