# /urs/bin/env Rscript

# RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types

# 利用RNA-SEQ表达数据，分析PBMC中免疫细胞类型的工具(t-sne)
if (!require("Rtsne")) {
  install.packages("Rtsne")
} else {
  require("Rtsne")
  print("Rtsne had been installed!")
}
library(RColorBrewer)


# 准备数据 --------------------------------------------------------------------
expr <- read.delim('/Users/congliu/Work/some_knowlege/GSE107011_Processed_data_TPM.txt.gz',
                   header = T, stringsAsFactors = F)
rownames(expr) <- expr[, 1]
expr <- subset(expr, select = -X)
# 只保留在3个及以上样本中有表达的基因
counts <- apply(expr, 1, function(x){sum(x != 0)})
expr.f <- expr[which(counts>3),]
expr.f <- log2(expr.f + 1)
expr.f <- t(expr.f)

# t-SNE -------------------------------------------------------------------
profvis::profvis({})
  #your code here
tsne_result <- Rtsne(
  expr.f,
  dims = 2,
  pca = FALSE,
  perplexity = 10, #小于(nrow-1)/3, 困惑度
  theta = 0.0, #速度与准确度间的平衡，0最精确
  max_iter = 1000
)
summary(tsne_result)
# 每个样本在两个维度的数据
tsne_result$Y
df <- as.data.frame(tsne_result$Y)
rownames(df) <- rownames(expr.f)
colnames(df) <- c("X","Y")
# 根据文献将行名拆分
names <- str_split(rownames(df),"_",simplify = T)
df$patient <- names[,1]
df$cell_types <- paste0(names[,2]," ", names[,3])
df$cell_types <- trimws(df$cell_types, which = "right")  # 去掉右边末端的空格
table(df$cell_types)
# 各panels的细胞类型
p1 <- c("PBMC","TFH","Treg","Th1","Th1.Th17","Th17","Th2","CD4 naive","CD4 TE")
p2 <- c("PBMC","VD2.","VD2..1","MAIT","CD8 naive","CD8 CM","CD8 EM","CD8 TE")
p3 <- c("PBMC","Progenitor","B naive","B NSM","B Ex","B SM","Plasmablasts")
p4 <- c("PBMC","Neutrophils","NK","C mono","I mono","NC mono","mDC","pDC","Basophils")

df_plot <- df
for (i in 1:4) {
  df_plot$types <- ifelse(df_plot$cell_types %in% get(paste0("p", i)), df_plot$cell_types, "Others")
  df_plot$types <- factor(df_plot$types, levels = c(get(paste0("p",i)),"others"), ordered = TRUE)
  assign(paste0("f",i),
         ggplot(df_plot, aes(x=X,y=Y,color=types))+
           geom_point()+
           scale_color_manual(values = c(brewer.pal(length(get(paste0("p",i))), "Set1"),"lightgrey")) + #指定颜色
           # 去掉网格线和x,y轴标签
           theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                 panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(),axis.line = element_line(colour = "black"))
  )
}

ggpubr::ggarrange(f1, f2, f3, f4, ncol = 2, nrow = 2, labels = paste0("Panel ",1:4))
