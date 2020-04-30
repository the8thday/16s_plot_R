#############���ڷ������Ե�΢���ﹲ��������
##����΢�����ȼ�����ϵ��
library(Hmisc)

#����ˮƽ���Ϊ������genus_table.txt�� ��һ����ˮƽ��΢�����ȱ�
genus <- read.delim('genus_table.txt', row.name = 1, check.names = FALSE)

#��ѡ���ȹ���һЩ�ͷ�Ȼ��Ƶ����Ⱥ
genus <- genus[which(rowSums(genus) >= 0.005), ]    #����ֻ������Է���ܺ͸��� 0.005 ����

genus1 <- genus
genus1[genus1>0] <- 1
genus <- genus[which(rowSums(genus1) >= 5), ]    #����ֻ������ 5 �������������г��ֵ���

#��������֮���Ƿ���ڷ�ȱ仯������ԣ��� spearman ���ϵ��Ϊ��
genus_corr <- rcorr(t(genus), type = 'spearman')

#��ֵɸѡ
#�� spearman ���ϵ������ 0.7 �Ĺ�ϵ�޳����� r>=0.7
r <- genus_corr$r
r[abs(r) < 0.7] <- 0

#ѡȡ������ p ֵС�� 0.05 �����ϵ������ p<0.05
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')    #��ѡ p ֵУ��������ʹ�� BH ��У�� p ֵ
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#��������ɸѡ�� r ֵ�� p ֵ��������
z <- r * p
diag(z) <- 0    #����ؾ����жԽ����е�ֵ������������أ�תΪ 0
head(z)[1:6,1:6]

#��˱�õ����ڽӾ����ʽ�������ļ���΢�����������ϵ������
write.table(data.frame(z, check.names = FALSE), 'genus_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##�������
library(igraph)

#���ڽӾ���ת��Ϊ igraph ������ڽ��б�
#������Ȩ���������磬Ȩ�ش�����΢���������ȵ� spearman ���ϵ�� 
g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
g

#�����Ҳ����ͨ����ʽȥ��
g <- simplify(g)

#�����ڵ��ɾ����ɾ����Ϊ 0 �Ľڵ㣩
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

#��ģʽ�£���Ȩ�ش��������ϵ��
#����Ȩ��ͨ��Ϊ��ֵ��������ȡ������ֵ�����ϵ�����¸���һ��
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

#Ϊ�ڵ㣨΢�����������������Ϣ�����Ÿ�Ŀ����ˮƽע�ͣ�
#��genus_taxonomy.txt�� ��¼��΢��������ԣ�����ñ�������֪����ڵ�ƥ���Ӧ����
tax <- read.delim('genus_taxonomy.txt', row.name = 1, check.names = FALSE, stringsAsFactors = FALSE)
tax <- tax[as.character(V(g)$name), ]

V(g)$kingdom <- tax$kingdom
V(g)$phylum <- tax$phylum
V(g)$class <- tax$class
V(g)$order <- tax$order
V(g)$family <- tax$family
V(g)$genus <- tax$genus

#�鿴����ͼ
g
plot(g)

##�����ļ����������ض��������ļ����ͣ����ں������ݷ�������
#�ڽӾ��󣬳��������ᵽ���ڼ������ϵ�������ɸѡ������ϵ��������
#�������� igraph ���ڽ��б�ת��
adj_matrix <- as.matrix(get.adjacency(g, attr = 'correlation'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

#���б�
edge <- data.frame(as_edgelist(g))    #igraph ���ڽ��б�תΪ���б�

edge_list <- data.frame(
    source = edge[[1]],
    target = edge[[2]],
    weight = E(g)$weight,
    correlation = E(g)$correlation
)
head(edge_list)

write.table(edge_list, 'network.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#�ڵ������б�
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

#���б�ڵ������б���Ե����� gephi �� cytoscape ��������ӻ�����н��б༭
#���� igraph Ҳ�ṩ�˿��Ա� gephi �� cytoscape ��ֱ��ʶ��ĸ�ʽ
#graphml ��ʽ����ʹ�� gephi ����򿪲����п��ӻ��༭
write.graph(g, 'network.graphml', format = 'graphml')

#gml ��ʽ����ʹ�� cytoscape ����򿪲����п��ӻ��༭
write.graph(g, 'network.gml', format = 'gml')

#############΢����ͻ�����΢����͹��ܼ�����Ե�̽��
#��΢����͹���Ϊ��

##����΢������Ⱥ��Ⱥ͹��ܻ����ȵ����ϵ��
library(Hmisc)

#����ˮƽ���Ϊ������phylum_table.txt�� ��һ����ˮƽ��΢�����ȱ�
phylum <- read.delim('phylum_table.txt', row.name = 1, check.names = FALSE)

#��ARGs_table.txt����һ�������ؿ��Ի����ȱ�
ARGs <- read.delim('ARGs_table.txt', row.name = 1, check.names = FALSE)

#����Ⱥ������빦�ܵ�����ԣ��� spearman ���ϵ��Ϊ��
phylum_ARGs_corr <- rcorr(as.matrix(phylum), as.matrix(ARGs), type = 'spearman')

#���ϵ�� r ֵ�������� p ֵ����
r <- phylum_ARGs_corr$r
p <- phylum_ARGs_corr$P

#ֻ����΢������-���ܻ����ȵ����ϵ��
#ȥ��΢����-΢������ܻ���-���ܻ���֮������ϵ��
r <- r[colnames(phylum),colnames(ARGs)]
p <- p[colnames(phylum),colnames(ARGs)]

#��ֵɸѡ
#�� spearman ���ϵ������ 0.7 �Ĺ�ϵ�޳����� r>=0.7
#��ģʽ�£�һ��Ҫע�⸺ֵ��ѡ���Ƿ��Ǻ��ʵģ���Ϊ�ܶ�����¿��ܸ����������
r[abs(r) < 0.7] <- 0

#ѡȡ������ p ֵС�� 0.05 �����ϵ������ p<0.05
p <- p.adjust(p, method = 'BH')    #��ѡ p ֵУ��������ʹ�� BH ��У�� p ֵ
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#��������ɸѡ�� r ֵ�� p ֵ��������
z <- r * p

#��ת��Ϊ�Գƾ���igraph ֻ��ʶ��������ʽ���ڽӾ�������
z1 <- phylum_ARGs_corr$r
z1[z1 != 0] <- 0
z1[rownames(z),colnames(z)] <- z
z1[colnames(z),rownames(z)] <- z

#write.table(data.frame(z1, check.names = FALSE), 'phylum_ARGs_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##�������
library(igraph)

#���ڽӾ���ת��Ϊ igraph ������ڽ��б�
#������Ȩ���������磬Ȩ�ش�����΢�����Ⱥ͹��ܻ����ȼ�� spearman ���ϵ�� 
g <- graph.adjacency(z1, weighted = TRUE, mode = 'undirected')
g

#�����ڵ��ɾ����ɾ����Ϊ 0 �Ľڵ㣩
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))

#��ģʽ�£���Ȩ�ش��������ϵ��
#����Ȩ��ͨ��Ϊ��ֵ��������ȡ������ֵ�����ϵ�����¸���һ��
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

#�鿴����ͼ
plot(g)

##�����ļ�������ԣ��ڽӾ��󡢱��б�Ȼ�ȡ��ʽ�ο����ļ���
#���� gml ��ʽ����ʹ�� cytoscape ����򿪲����п��ӻ��༭
write.graph(g, 'network.gml', format = 'gml')
