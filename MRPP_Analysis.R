#! /urs/bin/env Rscript

# MRPP（Multi Response Permutation Procedure）分析，是一种用于分析高维度数据组间相似性的统计方法。
# 它以距离矩阵为基础，基于组内和组间变异程度的置换检验，用于判断两组或两组以上实体之间的差异
# 例如，在群落分析中，常使用MRPP结合排序分析一起，评估不同环境的群落组成结构差异是否显著

library(tidyverse)
library(vegan)

infile <- '/Users/congliu/prs/IBD_project/combined_otu_table_m2_std.txt'
map_file <- '/Users/congliu/prs/IBD_project/map.txt'
outpath <- '/Users/congliu/prs/IBD_project/'
one_group <- c('IBD', 'Normal')

otu <- read.delim(infile, row.names = 1, sep = '\t', 
                  stringsAsFactors = FALSE, check.names = FALSE)
mapfile <- read.delim(map_file, sep = '\t') %>% dplyr::select(Type, Description)
otu <- subset(otu, taxonomy != 'Unassigned')
otu <- subset(otu, select = -taxonomy)
otu <- data.frame(t(otu))

??mrpp
# 如下是另外一批数据的测试，相同样本不同的丰度和环境数据进行一致性对比
mrpp_result <- mrpp(otu, group$site, distance = 'bray', permutations = 999)
mrpp_result
# Significance of delta，即显著性p值，小于0.05说明差异显著
# Chance corrected within-group agreement A，即A值，大于0说明组间差异大于组内差异，小于0说明组内差异大于组间差异；
# observe delta值越小说明组内差异小，expect delta值越大说明组间差异大

write.table(data.frame(Group = 'all', Distance = 'Bray-Curtis', A = mrpp_result$A,
                       Observe_delta = mrpp_result$delta, Expect_delta = mrpp_result$E.delta, P_value = mrpp_result$Pvalue),
            file = 'MRPP.result_all.txt', row.names = FALSE, sep = '\t', quote = FALSE)


# 以上对所有的分组进行了分析；可以分别挑选出相应的分组进行分析






