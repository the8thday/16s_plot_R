##网络格式转换

library(tidyverse)
library(igraph)
library(GGally)


# 合并相关矩阵和p值矩阵
#观测值的相关矩阵
cor_sparcc <- read.delim('/Users/congliu/prs/hbrm_Alveoli/cor_sparcc.out.txt', row.names = 1, sep = '\t', check.names = FALSE)
#伪 p 值矩阵
pvals <- read.delim('/Users/congliu/prs/hbrm_Alveoli/pvals.two_sided.txt', row.names = 1, sep = '\t', check.names = FALSE)
#保留 |相关性|≥0.8 且 p<0.01的值
cor_sparcc[abs(cor_sparcc) < 0.8] <- 0
pvals[pvals>=0.01] <- -1
pvals[pvals<0.01 & pvals>=0] <- 1
pvals[pvals==-1] <- 0

#筛选后的邻接矩阵
adj <- as.matrix(cor_sparcc) * as.matrix(pvals)
diag(adj) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
ggcorr(adj)
write.table(data.frame(adj, check.names = FALSE), '/Users/congliu/prs/hbrm_Alveoli/network.adj.txt', col.names = NA, sep = '\t', quote = FALSE)

# 后续进行R语言 igraph绘图
#输入数据，邻接矩阵
neetwork_adj <- read.delim('/Users/congliu/prs/hbrm_Alveoli/network.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)
head(neetwork_adj)[1:6]    #邻接矩阵类型的网络文件

#邻接矩阵 -> igraph 的邻接列表，获得含权的无向网络
g <- graph_from_adjacency_matrix(as.matrix(neetwork_adj), mode = 'undirected', weighted = TRUE, diag = FALSE)
g    #igraph 的邻接列表

#这种转换模式下，默认的边权重代表了 sparcc 计算的相关性（存在负值）
#由于边权重通常为正值，因此最好取个绝对值，相关性重新复制一列作为记录
E(g)$sparcc <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
is.weighted(g)
plot(g, layout=layout.circle)



# 其他类型网络文件 ----------------------------------------------------------------

#再转为其它类型的网络文件，例如
#再由 igraph 的邻接列表转换回邻接矩阵
adj_matrix <- as.matrix(get.adjacency(g, attr = 'sparcc'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(g, 'network.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(g, 'network.gml', format = 'gml')


# 边列表，节点属性列表 --------------------------------------------------------------

#边列表，也可以直接导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
edge <- data.frame(as_edgelist(g))

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  sparcc = E(g)$sparcc
)
head(edge_list)

write.table(edge_list, 'network.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#节点属性列表，对应边列表，记录节点属性，例如
node_list <- data.frame(
  nodes_id = V(g)$name,    #节点名称
  degree = degree(g)    #节点度
)
head(node_list)

write.table(node_list, 'network.node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)


