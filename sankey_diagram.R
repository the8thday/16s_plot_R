# /urs/bin/env Rscript
# sanky plot

# 参考 https://mp.weixin.qq.com/s/dgxgi_3PdjW5g-fOPlAe4Q?
  
#site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# install.packages("ggalluvial", repo=site)
# library(Seurat)
library(ggalluvial)
# library(riverplot)


# pbmc_small <- FindClusters(
#   object = pbmc_small,
#   resolution = c(seq(.4,1.6,.2))
# )
ggplot(data = pbmc_small@meta.data,
       aes(axis1 = RNA_snn_res.0.4, axis2 = RNA_snn_res.0.6,axis3 = RNA_snn_res.0.8,axis4 = RNA_snn_res.1,
           axis5 = RNA_snn_res.1.2,axis6 = RNA_snn_res.1.4,axis7 = RNA_snn_res.1.6)) +
  scale_x_discrete(limits = c(paste0("RNA_snn_res.",seq(.4,1.6,.2))), expand = c(.01, .05)) +
  geom_alluvium(aes(fill = RNA_snn_res.1.6)) +
  geom_stratum() + geom_text(stat = "stratum", infer.label = TRUE) +
  #coord_polar()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ggtitle("cell number in each cluster")


# 测试 ----------------------------------------------------------------------

titanic_wide <- data.frame(Titanic)

ggplot(data = titanic_wide,
       aes(axis1 = Class, axis2 = Sex, axis3 = Age,
           weight = Freq)) +
  scale_x_discrete(limits = c("Class", "Sex", "Age"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Survived)) +
  geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
  theme_minimal() +
  ggtitle("passengers on the maiden voyage of the Titanic",
          "stratified by demographics and survival")

titanic_long <- to_lodes(data.frame(Titanic),
                         key = "Demographic",
                         axes = 1:3)
head(titanic_long)
ggplot(data = titanic_long,
       aes(x = Demographic, stratum = stratum, alluvium = alluvium,
           weight = Freq, label = stratum)) +
  geom_alluvium(aes(fill = Survived)) +
  geom_stratum() + geom_text(stat = "stratum") +
  theme_minimal() +
  ggtitle("passengers on the maiden voyage of the Titanic",
          "stratified by demographics and survival")
# 产生和上图一样的图，只是数据源格式不同


ggplot(as.data.frame(Titanic),
       aes(weight = Freq,
           axis1 = Survived, axis2 = Sex, axis3 = Class)) +
  geom_alluvium(aes(fill = Class),
                width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", label.strata = TRUE, reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("Survived", "Sex", "Class")) +
  coord_flip() +
  ggtitle("Titanic survival by class and sex")


# -------------------------------------------------------------------------

# 实战1. 组间丰度变化 

# 编写测试数据
df=data.frame(
  Phylum=c("Ruminococcaceae","Bacteroidaceae","Eubacteriaceae","Lachnospiraceae","Porphyromonadaceae"),
  GroupA=c(37.7397,31.34317,222.08827,5.08956,3.7393),
  GroupB=c(113.2191,94.02951,66.26481,15.26868,11.2179),
  GroupC=c(123.2191,94.02951,46.26481,35.26868,1.2179),
  GroupD=c(37.7397,31.34317,222.08827,5.08956,3.7393)
)

# 数据转换长表格

melt_df = reshape2::melt(df)

# 绘制分组对应的分类学，有点像circos
ggplot(data = melt_df,
       aes(axis1 = Phylum, axis2 = variable,
           weight = value)) +
  scale_x_discrete(limits = c("Phylum", "variable"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Phylum)) +
  geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
  theme_minimal() +
  ggtitle("Phlyum abundance in each group")


# 组间各丰度变化 
ggplot(data = melt_df,
       
       aes(x = variable, weight = value, alluvium = Phylum)) +
  geom_alluvium(aes(fill = Phylum, colour = Phylum),
                alpha = .75, decreasing = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -30, hjust = 0)) +
  ggtitle("Phylum change among groups")








