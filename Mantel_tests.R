#! /urs/bin/env Rscript

library(tidyverse)
library(vegan)

# Mantel tests是确定两组距离测度矩阵（而非两组变量矩阵）之间相关性的相关性测试方法，
# 用于判断一个矩阵中的样本距离与另一矩阵中的样本距离是否相关。
# Qiime 似乎提供了脚本
# 相关不等于因果，相关只是表明有关系
# https://jkzorz.github.io/2019/07/08/mantel-test.html

env <- read.delim('/Users/congliu/prs/ZKZL_analysis/saliva/BMI/BMI.txt',
                  sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

otu <- read.delim('/Users/congliu/prs/ZKZL_analysis/saliva/BMI/combined_otu_table_m2_std.txt',
                  sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

group <- read.delim('/Users/congliu/prs/ZKZL_analysis/saliva/BMI/map.txt',
                    sep = '\t', stringsAsFactors = FALSE)

outpath <- '/Users/congliu/prs/ZKZL_analysis/saliva/BMI/'

#env <- env %>% select(c(PFS))
env <- subset(env, select = BMI)
#env$type <- ifelse(env$Type == 'nPR', 0, 1)
otu <- otu %>% select(-c(taxonomy))
otu <- t(otu)
otu <- log(otu + 1)

otu.dist <- vegdist(otu, method = 'bray')
env.scale <- scale(env$PFS, center = TRUE, scale = TRUE)

# 对于环境因素其距离该如何选择？不同的量纲，此处采用欧式距离
env.dist <- dist(env$BMI, method = 'euclidean')


# Mantel Test && Plot -----------------------------------------------------
# ??mantel
# permutation可视具体情况进行调整
mantel.result <- mantel(otu.dist, env.dist, method = 'spearman', permutations = 9999, na.rm = TRUE)
mantel.result
# 结果中 Mantel statistic r:  0.76； Significance: 0.001，otu矩阵和环境距离之间有很强的相关性，P值拒绝H0假设（距离矩阵间无相关性）
pvalue <- round(mantel.result$signif,3)
r2 <- round(mantel.result$statistic,3)
otu_distance <- 'bray'
env_distance <- 'euclidean'
mantel_res <- as.data.frame(cbind(otu_distance, env_distance, pvalue, r2))
write.table(mantel_res, file = file.path(outpath, 'Mantel_result.txt'), row.names = FALSE, 
            sep = '\t', quote = FALSE, na = '')


# 重要的是分析，比如环境因素如何影响丰度的变化，需要具体结合
df_p <- cbind(env, otu)
p <- ggplot(df_p, aes(x = env1, y = species1)) +
  geom_smooth(stat = 'smooth', method = 'lm') +
  geom_point(aes(colour = Latitude), size = 4) +
  labs(y = 'Pelagibacteraceae (OTU 307744) (%)', x = 'Temperature (C)') + 
  theme( axis.text.x = element_text(face = 'bold',colour = 'black', size = 12), 
         axis.text.y = element_text(face = 'bold', size = 11, colour = 'black'), 
         axis.title= element_text(face = 'bold', size = 14, colour = 'black'), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = 'black'), 
         legend.title = element_text(size =12, face = 'bold', colour = 'black'),
         legend.text = element_text(size = 10, face = 'bold', colour = 'black')) +
  scale_colour_continuous(high = 'navy', low = 'salmon')




















