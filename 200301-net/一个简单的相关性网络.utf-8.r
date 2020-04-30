#############基于丰度相关性的微生物共发生网络
##计算微生物丰度间的相关系数
library(Hmisc)

#以属水平丰度为例，“genus_table.txt” 是一个属水平的微生物丰度表
genus <- read.delim('genus_table.txt', row.name = 1, check.names = FALSE)

#可选事先过滤一些低丰度或低频的类群
genus <- genus[which(rowSums(genus) >= 0.005), ]    #例如只保留相对丰度总和高于 0.005 的属

genus1 <- genus
genus1[genus1>0] <- 1
genus <- genus[which(rowSums(genus1) >= 5), ]    #例如只保留在 5 个及以上样本中出现的属

#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
genus_corr <- rcorr(t(genus), type = 'spearman')

#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
r <- genus_corr$r
r[abs(r) < 0.7] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(z, check.names = FALSE), 'genus_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##获得网络
library(igraph)

#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数 
g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
g

#自相关也可以通过该式去除
g <- simplify(g)

#孤立节点的删除（删除度为 0 的节点）
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

#为节点（微生物属）添加属性信息（界门纲目科属水平注释）
#“genus_taxonomy.txt” 记录了微生物的属性，读入该表后根据已知网络节点匹配对应的行
tax <- read.delim('genus_taxonomy.txt', row.name = 1, check.names = FALSE, stringsAsFactors = FALSE)
tax <- tax[as.character(V(g)$name), ]

V(g)$kingdom <- tax$kingdom
V(g)$phylum <- tax$phylum
V(g)$class <- tax$class
V(g)$order <- tax$order
V(g)$family <- tax$family
V(g)$genus <- tax$genus

#查看网络图
g
plot(g)

##网络文件输出，输出特定的网络文件类型，便于后续数据分析需求
#邻接矩阵，出了上述提到的在计算相关系数后，输出筛选后的相关系数矩阵外
#还可以由 igraph 的邻接列表转换
adj_matrix <- as.matrix(get.adjacency(g, attr = 'correlation'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

#边列表
edge <- data.frame(as_edgelist(g))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
    source = edge[[1]],
    target = edge[[2]],
    weight = E(g)$weight,
    correlation = E(g)$correlation
)
head(edge_list)

write.table(edge_list, 'network.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#节点属性列表
node_list <- data.frame(
    label = names(V(g)),
    kingdom = V(g)$kingdom,
    phylum = V(g)$phylum,
    class = V(g)$class,
    order = V(g)$order,
    family = V(g)$family,
    genus = V(g)$genus
)
head(node_list)

write.table(node_list, 'network.node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#边列表节点属性列表可以导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
#此外 igraph 也提供了可以被 gephi 或 cytoscape 等直接识别的格式
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(g, 'network.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(g, 'network.gml', format = 'gml')

#############微生物和环境、微生物和功能间相关性的探索
#以微生物和功能为例

##计算微生物类群丰度和功能基因丰度的相关系数
library(Hmisc)

#以门水平丰度为例，“phylum_table.txt” 是一个门水平的微生物丰度表
phylum <- read.delim('phylum_table.txt', row.name = 1, check.names = FALSE)

#“ARGs_table.txt”是一个抗生素抗性基因丰度表
ARGs <- read.delim('ARGs_table.txt', row.name = 1, check.names = FALSE)

#计算群落组成与功能的相关性，以 spearman 相关系数为例
phylum_ARGs_corr <- rcorr(as.matrix(phylum), as.matrix(ARGs), type = 'spearman')

#相关系数 r 值和显著性 p 值矩阵
r <- phylum_ARGs_corr$r
p <- phylum_ARGs_corr$P

#只保留微生物丰度-功能基因丰度的相关系数
#去除微生物-微生物、功能基因-功能基因之间的相关系数
r <- r[colnames(phylum),colnames(ARGs)]
p <- p[colnames(phylum),colnames(ARGs)]

#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
#该模式下，一定要注意负值的选择是否是合适的，因为很多情况下可能负相关无意义
r[abs(r) < 0.7] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p

#再转换为对称矩阵，igraph 只能识别这种样式的邻接矩阵类型
z1 <- phylum_ARGs_corr$r
z1[z1 != 0] <- 0
z1[rownames(z),colnames(z)] <- z
z1[colnames(z),rownames(z)] <- z

#write.table(data.frame(z1, check.names = FALSE), 'phylum_ARGs_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##获得网络
library(igraph)

#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物丰度和功能基因丰度间的 spearman 相关系数 
g <- graph.adjacency(z1, weighted = TRUE, mode = 'undirected')
g

#孤立节点的删除（删除度为 0 的节点）
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

#查看网络图
plot(g)

##网络文件输出，略，邻接矩阵、边列表等获取方式参考上文即可
#例如 gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(g, 'network.gml', format = 'gml')
