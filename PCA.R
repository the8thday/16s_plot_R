#! /urs/bin/env Rscript
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials

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


infile <- '/Users/congliu/prs/ruijin_xirou/combined_otu_table_m2_std.txt'
map_file <- '/Users/congliu/prs/ruijin_xirou/map.txt'
outpath <- '/Users/congliu/prs/ruijin_xirou/'
one_group <- c('polyp_stool', 'unpolyp_stool')
group1 <- one_group[1]
group2 <- one_group[2]

otu <- read.delim(infile, row.names = 1, sep = '\t', 
                  stringsAsFactors = FALSE, check.names = FALSE)
mapfile <- read.delim(map_file, sep = '\t') %>% dplyr::select(Type, Description)
otu <- subset(otu, select = -taxonomy)
otu <- data.frame(t(otu))
# PCA处理数据
otu2 <- log(otu + 1)
#otu2 <- otu
otu2$Description <- row.names(otu2)
otu3 <- merge(otu2, mapfile, by.x = 'Description', by.y = 'Description')
row.names(otu3) <- otu3$Description
otu3 <- subset(otu3, Type %in% one_group)
otu4 <- subset(otu3, select = -c(Description, Type))
otu5 <- scale(otu)

res.pca <- PCA(otu4, scale.unit = FALSE, graph = FALSE)
eig.val <- get_eigenvalue(res.pca)
#fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 30))
var.values <- get_pca_var(res.pca)
var.values$contrib #通过PCA的方式选择变量
#fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 10)
#fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
res.desc <- dimdesc(res.pca, axes = c(1,2), proba = 0.05)


# Graph of individuals ----------------------------------------------------

ind <- get_pca_ind(res.pca)
ind.plot <- fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = otu3$Type, # color by groups
             fill.ind = otu3$Type,
             #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             palette = 'npg',
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "convex",
             ellipse.linetype = 2,
             legend.title = "Groups"
)

ggpubr::ggpar(ind.plot,
              title = "Principal Component Analysis",
              caption = "Source: prs",
              xlab = "PC1", ylab = "PC2",
              legend.title = "Group", legend.position = "top",
              ggtheme = theme_gray(), palette = "npg"
)


