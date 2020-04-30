#！/urs/bin/env Rscript
# ALDEx2是一种常用的微生物组差异分析方法，
# 1.用原始输入数据生成每个分类单元的后验概率分布；然后将该分布进行中心对数变换
# 2.将变换后的值，用参数或非参数检验进行单变量统计检验，并返回 p 值和 Benjamini-Hochberg 校正后的 p 值


if (!requireNamespace("ALDEx2", quietly = TRUE))
  install.packages("ALDEx2")

library(ALDEx2)

