#！/urs/bin/env Rscript

# 检验数据集是否为正态分布，方差齐性

infile <- '/Users/congliu/prs/ruijin_xirou/combined_otu_table_m2_std.txt'
map_file <- '/Users/congliu/prs/ruijin_xirou/map.txt'
outpath <- '/Users/congliu/prs/ruijin_xirou/'


otu <- read.delim(infile, row.names = 1, sep = '\t', 
                  stringsAsFactors = FALSE, check.names = FALSE) 
otumat <- subset(otu, taxonomy != 'Unassigned', select = -taxonomy)



# shapiro -----------------------------------------------------------------
# 频率直方图
# qq plot 是将两组分布的百分位树进行比对，比如对于一组数据和正态分布进行对比，或者
qqnorm(otumat$PRS029180041)
car::qqPlot(lm(chao1 ~ Type, data = input_data), simulate = TRUE, main = 'QQ Plot', labels = TRUE)
# 似乎对样本数目有要求, <5000,若结果中p值大于0.05，则接受原假设，数据分布符合正态性
shapiro.test(input_data$chao1)
# p值小于0.05 认为符合正态分布
apply(otumat, 2, shapiro.test)
shapiro.test(log(otumat$PRS029180041 + 1))
# 微生物组的数据一般不符合正态分布


# Kolmogorov - Smirnov ----------------------------------------------------
# Kolmogorov-Smirnov检验只能检验是否一个样本来自于一个已知样本，而Lilliefor检验可以检验是否来自未知总体
# 一种非参数检验的方法
ks.test(otumat$PRS029180041, 'pnorm')

# 方差齐性
car::leveneTest(otumat$PRS029180041, )
# 对于已经通过正态性检验的数据，推荐使用Bartlett检验来进行方差齐性检验
bartlett.test(chao1 ~ Type, data = input_data)
# 
fligner.test(chao1 ~ Type, data = input_data)


# 常见的统计方法 -----------------------------------------------------------------

# t-test
# anova
# 



# 相关性相关 -------------------------------------------------------------------

##Pearson、Spearman、Kendall 相关

#标准化不影响相关系数计算值，但可以让数据服从均值 0，标准差 1 的等方差结构
mtcars <- scale(mtcars)

#协方差计算，cov()
cov_pearson <- cov(mtcars, method = 'pearson')
cov_pearson

cov_spearman <- cov(mtcars, method = 'spearman')
cov_spearman

cov_kendall <- cov(mtcars, method = 'kendall')
cov_kendall

#相关系数计算，cor()
cor_pearson <- cor(mtcars, method = 'pearson')
cor_pearson

cor_spearman <- cor(mtcars, method = 'spearman')
cor_spearman

cor_kendall <- cor(mtcars, method = 'kendall')
cor_kendall

#相关图，例如
library(corrplot)

corrplot(cor_pearson, method = 'number', number.cex = 0.8, diag = FALSE, tl.cex = 0.8)
corrplot(cor_pearson, add = TRUE, type = 'upper', method = 'pie', diag = FALSE, tl.pos = 'n', cl.pos = 'n')


