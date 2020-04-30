#! /urs/bin/env Rscript

# mRMR algorithm

# library(sideChannelAttack)

library(tidyverse)
library(mRMRe)

m2 <- '/Users/congliu/prs/lasso_data/crc2_4_45_46_47_48_49_mv34_qiime_silva_v2_m2_std.txt'
map_file <- '/Users/congliu/prs/lasso_data/map_20200108_add_ref_n_a_c.txt'
one_group <- c('c', 'n')

otu <- read.delim(m2, row.names = 1, sep = '\t', 
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
otu3$Type <- as.character(otu3$Type)
otu3 <- subset(otu3, Type %in% one_group)
otu3$Type <- ordered(as.factor(otu3$Type))
otu4 <- subset(otu3, select = -c(Description))
otu4$Type <- c('n'=0, 'c'=1)[as.character(otu4$Type)]
otu4$Type <- ordered(as.factor(otu4$Type)) #在将Type列作为factor的时候，出现了问题，难道只能用于连续性的数值？
#sapply(otu4, class)
#otu4 <- transform(otu4, X.Output. = as.numeric(X.Output.))
#write_delim(otu4, '~/hehe.txt', delim = '\t')

otu5 <- subset(otu3, select = -c(Description, Type))
# Begin -------------------------------------------------------------------

set.thread.count(2)

data(cgps)
data.annot <- data.frame(cgps.annot)
data.cgps <- data.frame(cgps.ic50, cgps.ge)
head(data.cgps)[, 1:5]
# 数据类型还能是因子 或 Surv
# dd <- mRMR.data(data = data.cgps)
# 真实数据
dd <- mRMR.data(data = otu4)
# correlate(X=dd[,1], Y=otu4$Type, method = 'cindex')
ddd <- subsetData(dd, 1:10, 1:10)
# Uses Spearman as correlation estimator
spearman_mim <- mim(ddd, continuous_estimator = "spearman")
# Uses Pearson as correlation estimator
pearson_mim <- mim(ddd, continuous_estimator = "pearson")

# 模型运行
mrmr_classic <- mRMR.classic(data = dd, target_indices = c(3136), feature_count = 30)

# ensemble mRMR
# For mRMR.classic-like results
mRMR.ensemble(data = dd, target_indices = c(1), solution_count = 1, feature_count = 30)
# For mRMR.ensemble-like results
ensemble <- mRMR.ensemble(data = dd, target_indices = c(3136), solution_count = 5, feature_count = 30)
causality(ensemble) #the minimum result is kept
featureCount(dd)
featureCount(mrmr_classic)
featureData(dd) # 好像有点问题
# 选择特征，如何设置solution的参数, 如何得到mi_threshold，
selected.features.mrmre <- solutions(mrmr_classic, 
                                     #mi_threshold = 0.01
                                     )
selected_feature <- mrmr_classic@feature_names[unlist(selected.features.mrmre)]
foo <- str_replace_all(selected_feature, '\\.\\.', '; ')
readr::write_lines(foo, 'selected_feature_by_mrmr.txt')

# network
feature_data <- mRMR.data(data =  data.frame(cgps.ge))
mim(feature_data)
network <- new('mRMRe.Network', data = dd, target_indices = c(1), 
               levels = c(3,1), layers =2)
network <- mRMR.network(data = feature_data, target_indices = c(1,2), levels = c(2,1), layers = 1)
mRMRe::visualize(network)
# adjacencyMatrix(network)



















