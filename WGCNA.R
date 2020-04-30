#！/urs/bin/env Rscript
# WGCNA

# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# http://pages.stat.wisc.edu/~yandell/statgen/ucla/WGCNA/wgcna.html
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
library(WGCNA)
# WGCNA 的聚类具有生物学意义(主要体现在其对相关性的加权分析)，而非常规的聚类算法；module为一组具有类似表达的基因
# 计算两个基因间的相关性，可人为定义相似性的阈值，为了减免人为所致的因素，WGCNA采用软阈值的方法规避了这一问题

# 本例的目的在于提取出与样本特征相关的基因共表达模块，比如和肠癌相关的物种菌属的模块


# The following setting is important, do not omit
options(stringsAsFactors = FALSE)
#enableWGCNAThreads() #似乎容易报错, 目前只支持R自带GUI

otu <- read.delim('/Users/congliu/Work/some_knowlege/WGCNA/ExpData.txt', header = T,
                  stringsAsFactors = F, 
                  row.names = 1)
mapfile <- read.delim('/Users/congliu/Work/some_knowlege/WGCNA/Sam_info.txt',
                      header = T,
                      stringsAsFactors = F,
                      row.names = 1)
# 数据为行为样本、列为特征
otu <- t(otu)
# 筛选一定比值的数据, 数据处理前的处理
otu.sums = apply(otu, 2, function(x){sum(x != 0)})
#otu <- otu[which(otu.sums > 3), ]
otu <- otu[, which(otu.sums > 3)]

notus = ncol(otu)
nSamples = nrow(otu)

# 样本聚类检查离群值
gsg = goodSamplesGenes(otu, verbose = 3)
if (gsg$allOK) {
  print("all sample test OK")
} else {
  print("your data need attention!")
  #stop("Sample test failed")
}

# 看一下样本的UPGMA聚类的分布
sampleTree = hclust(dist(otu), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers"
     , sub="", xlab="")

# 如果需要去除离散值, 离散样本
# clust = cutreeStatic(sampleTree, cutHeight = 20000, minSize = 10)
# table(clust)
# keepSamples = (clust==1)
# datExpr = datExpr[keepSamples, ]
# nGenes = ncol(datExpr)
# nSamples = nrow(datExpr)

# 软阈值筛选
#powers = c(c(1:10), seq(from = 12, to=20, by=2))
powers = 1:20
# 计算无尺度分布拓扑矩阵
sft = pickSoftThreshold(otu, powerVector = powers, verbose = 5)
sft$fitIndices # 不同power值所对应的参数
sft$powerEstimate # 最佳阈值
# 绘制拟合指数和power散点图；平均连通性和power散点图
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1, col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# sft$powerEstimate 值NA，两批次微生物组的数据都显示为NA，是不是数据需要什么处理
# 一个合适的power值的选择
power = 13
# power = sft$powerEstimate
##一步法网络构建：One-step network construction and module detection##
otu_name <- row.names(otu)
otu <- as.data.frame(sapply(otu, as.numeric))
row.names(otu) <- otu_name
net = blockwiseModules(otu, power = power, maxBlockSize = 6000,
                       TOMType = "unsigned", 
                       minModuleSize = 30, #每个基因模块最少的基因数目为
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       #saveTOMFileBase = "otu-TOM",
                       verbose = 3)
# 分步获得TOM矩阵
# adjacency <- adjacency(otu, power = power)
# tom_sim <- TOMsimilarity(adjacency)
# rownames(tom_sim) <- rownames(adjacency)
# colnames(tom_sim) <- colnames(adjacency)
# tom_sim[1:5,1:5] # 即为基因与基因间的tom矩阵，越接近1，共表达相似性越高
# #TOM 相异度 = 1 – TOM 相似度
# tom_dis  <- 1 - tom_sim
# #层次聚类树，使用中值的非权重成对组法的平均聚合聚类
# geneTree <- hclust(as.dist(tom_dis), method = 'average')
# plot(geneTree, xlab = '', sub = '', main = 'Gene clustering on TOM-based dissimilarity',
#      labels = FALSE, hang = 0.04)
#使用动态剪切树挖掘模块
# minModuleSize <- 30  #模块基因数目
# dynamicMods <- cutreeDynamic(dendro = geneTree, distM = tom_dis,
#                              deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
# 
# table(dynamicMods)
# #模块颜色指代
# dynamicColors <- labels2colors(dynamicMods)
# table(dynamicColors)
# plotDendroAndColors(geneTree, dynamicColors, 'Dynamic Tree Cut',
#                     dendroLabels = FALSE, addGuide = TRUE, hang = 0.03, guideHang = 0.05,
#                     main = 'Gene dendrogram and module colors')
# #基因表达聚类树和共表达拓扑热图
# plot_sim <- -(1-tom_sim)
# #plot_sim <- log(tom_sim)
# diag(plot_sim) <- NA
# TOMplot(plot_sim, geneTree, dynamicColors,
#         main = 'Network heatmap plot, selected genes')
# 模块特征基因并非某个具体的基因，而是一种“特征组合”


##绘画结果展示模块和特征聚类树
### open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
tom_tree <- net$dendrograms[[1]]
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# net对象中的net$dendrograms是基于bray距离的？
datExpr_tree<-hclust(dist(otu), method = "average")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)

## 如果这个时候样本是有性状，或者临床表型的，可以加进去看看是否聚类合理
#针对前面构造的样品矩阵添加对应颜色
# mapfile
sample_colors <- numbers2colors(as.numeric(factor(mapfile$Type)), 
                                colors = c("white","blue","red","green"),signed = FALSE)
sample_colors <- numbers2colors(mapfile ,signed = FALSE)
## 如果样品有多种分类情况，而且 datTraits 里面都是分类信息，那么可以直接用上面代码，当然，这样给的颜色不明显，意义不大。
#构造10个样品的系统聚类树及性状热图
par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(datExpr_tree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")

##结果保存###
# moduleLabels = net$colors
# moduleColors = labels2colors(net$colors)
# table(moduleColors)
# MEs = net$MEs;
# geneTree = net$dendrograms[[1]];
# save(MEs, moduleLabels, moduleColors, geneTree,
#      file = "otu.RData")


##表型与模块相关性##
## 这一步主要是针对于连续变量，如果是分类变量，需要转换成连续变量方可使用
moduleColorsWW = mergedColors
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(otu, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
# 在WGCNA中ME、
samples = model.matrix(~0+ mapfile$D1)
samples = mapfile
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

#par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.9,  yColorWidth=0.01, 
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
###导出网络到Cytoscape

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(otu, power = 13)
geneTree = net$dendrograms[[1]]
dissTOM = 1-TOM
nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(notus, size = nSelect)
selectTOM = dissTOM[select, select]
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

# Recalculate module eigengenes
# 模块特征基因
MEList = moduleEigengenes(otu, mergedColors)
MEs = MEList$eigengenes
# 通过模块特征基因进一步对共表达模块进行聚类







# 感兴趣性状的模块的具体基因分析





















