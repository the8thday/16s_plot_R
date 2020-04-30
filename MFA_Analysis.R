#! /urs/bin/env Rscript


if (!require(FactoMineR) & !require(factoextra)){
  install.packages(c("FactoMineR", "factoextra"))
} else {
  library("FactoMineR")
  library("factoextra")
}
library(tidyverse)
#library(ade4)
#library("corrplot")
library(ggsci)


# MFA是一种探索性的对称排序方法，不是因果关系的假设检验，MFA对数据输入格式比较宽松