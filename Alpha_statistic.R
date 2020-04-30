#！/urs/bin/env Rscript

## ---------------------------
##
## Script name: Alpha_statistic.R
##
## Purpose of script: alpha指数统计检验
##
## Author: liuc
##
## Date Created: 2020-04-29
##
## Copyright (c) 
## Email:
##
## ---------------------------
##
## Notes:
## Usage: Rscript Alpha_statistic.R m2.txt map.txt outpath   
##
## ---------------------------

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
infile <- args[1]
map_file <- args[2]
outpath <- args[3]

library(tidyverse)
library(ggrepel)
library(ggsci)


# infile <- '/Users/congliu/prs/WXRM/Reporter/combined_alpha_m1.txt'
# map_file <- '/Users/congliu/prs/WXRM/Reporter/map.txt'
# outpath <- '/Users/congliu/prs/WXRM/Reporter'
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# data tidy ---------------------------------------------------------------

alpha <- read_delim(infile, delim = '\t')
mapfile <- read_delim(map_file, delim = '\t') %>% 
  select(Type, Description)

input_data <- left_join(alpha, mapfile, by = c('Sample'='Description'))
fenzu <- unique(input_data$Type)
# fenzu <- c('metaplasia', 'Inflammation', 'cancer')

input_data <- input_data %>% filter(Type %in% fenzu)
# 可选择的level顺序
# input_data$Type <- ordered(input_data$Type, levels=c('Inflammation', 'metaplasia', 'cancer'))

# T分布检验和方差齐性检验 ------------------------------------------------------------

for (i in c('chao1', 'ace', 'shannon', 'simpson', 'goods_coverage')){
  shapiro_res <- tapply(input_data[[i]], input_data$Type, shapiro.test)
  shapiro_res
}

# shapiro_res_final <- data.frame(shapiro_res) %>% rownames_to_column() %>% 
#   tidyr::extract(rowname, into=c('Type', 'statistic'), regex="([[:alnum:]]+)\\.([[:alnum:]]+).*?")
shapiro_res_final <- data.frame(unlist(shapiro_res)) %>% rownames_to_column() %>% 
  separate(rowname, into=c('Type', 'statistic', 'X'), sep = '\\.') %>% 
  replace_na(list(X = 'M')) %>% mutate(zhi = unlist.shapiro_res.) %>% 
  select(-unlist.shapiro_res., -statistic) %>% filter(!zhi == 'X[[i]]') %>% 
  pivot_wider(names_from = X, values_from = zhi)

names(shapiro_res_final) <- c('Type', 'W', 'P', 'Method')

# 多分组levene检验中，Type各组间方差是否有差异
levene_all <- car::leveneTest(chao1 ~ Type, data = input_data)

levene_list <- NULL
for (i in seq(length(fenzu)-1)) {
  for (j in (i + 1):length(fenzu)) {
    input_data_ij <- subset(input_data, Type %in% c(fenzu[i], fenzu[j]))
    group_1 <- fenzu[i]
    group_2 <- fenzu[j]
    for (alpha in c('chao1', 'ace', 'shannon', 'simpson', 'goods_coverage')) {
      fml <- formula(paste(alpha,"~ Type"))
      levene_res <- car::leveneTest(fml, data = input_data_ij)
      levene_list <- rbind(levene_list,c(
        alpha,
        paste(group_1, group_2, sep = '/'),
        levene_res$Df[[1]],
        levene_res$`F value`[[1]],
        levene_res$`Pr(>F)`[[1]]
      )
      )
    }
  }
}
levene_list <- data.frame(levene_list)
names(levene_list) <- c('Type', 'Group', 'Df', 'F', 'P')
write.table(levene_list, file = file.path(outpath, 'Levene_test.txt'), sep = '\t',
            row.names = FALSE, quote = FALSE, na = '')

# Wilcoxon rank sum test or Mann-Whitney test -----------------------------

wilcoxon_list <- NULL
# pairwise.wilcox.test(input_data$chao1, input_data$Type, p.adjust.method = "BH")

for (i in seq(length(fenzu)-1)) {
  for (j in (i + 1):length(fenzu)) {
    cat(i, j, '\n')
    input_data_ij <- subset(input_data, Type %in% c(fenzu[i], fenzu[j]))
    group_1 <- fenzu[i]
    group_2 <- fenzu[j]
    for (alpha in c('chao1', 'ace', 'shannon', 'simpson', 'goods_coverage')) {
      fml <- formula(paste(alpha,"~ Type"))
      wilcox_res_ij <- wilcox.test(fml, data = input_data_ij, 
                                alternative = 'two.sided', paired = FALSE, exact=TRUE)
      wilcoxon_list <- rbind(wilcoxon_list, c(
        alpha,
        paste(group_1, group_2, sep = '/'),
        wilcox_res_ij$method,
        wilcox_res_ij$alternative,
        wilcox_res_ij$statistic,
        wilcox_res_ij$p.value
      ))
    }
    
  }
}
wilcoxon_list <- data.frame(wilcoxon_list)
names(wilcoxon_list) <- c('Type', 'Groups', 'Method', 'Alternative', 'W', 'P')
write.table(wilcoxon_list, file = file.path(outpath, 'Wilcoxon_rank_sum.txt'), sep = '\t',
            row.names = FALSE, quote = FALSE, na = '')


# Kruskal-Wallis Test -----------------------------------------------------

kruskal_list <- list()
kruskal_final <- NULL

for (i in c('chao1', 'ace', 'shannon', 'simpson', 'goods_coverage')){
  fml <- formula(paste(i,"~ Type"))
  kruskal_res <- kruskal.test(fml, data = input_data)
  kruskal_list[[i]] <- kruskal_res
  kruskal_final <- rbind(kruskal_final, c(
    i,
    kruskal_res$method,
    kruskal_res$statistic,
    kruskal_res$p.value
  ))
}

kruskal_final_2 <- data.frame(kruskal_final)
names(kruskal_final_2) <- c('Type', 'Method', 'chi_squared', 'P')
write.table(kruskal_final_2, file = file.path(outpath, 'Kruskal_Wallis_Test.txt'), sep = '\t',
            row.names = FALSE, quote = FALSE, na = '')

ggplot(data = input_data) +
  #geom_violin(aes(x = Type, y = chao1), trim = FALSE)+
  geom_boxplot(aes(x = Type, y = chao1, color = Type)) +
  scale_color_nejm() +
  theme_bw()
  
  # geom_segment(aes(x = Type, y = max(chao1), xend = Type, yend = ))
  # annotate('text')

# one-way analysis of variance  -------------------------------------------
# 没有去除离散点
one_way_anova <- NULL
tukey_res <- NULL

for (i in c('chao1', 'ace', 'shannon', 'simpson', 'goods_coverage')){
  fml <- formula(paste(i,"~ Type"))
  anova_res <- aov(fml, data = input_data)
  one_way_anova <- rbind(one_way_anova, c(
    i,
    summary(anova_res)[[1]][1,1],
    summary(anova_res)[[1]][1,4],
    summary(anova_res)[[1]][['Pr(>F)']][[1]]
  ))
  if (summary(anova_res)[[1]][['Pr(>F)']][[1]] < 0.05) {
    cat('You need do TukeyHSD')
  }
}
one_way_anova_2 <- data.frame(one_way_anova)
names(one_way_anova_2) <- c('Type', 'Df', 'F', 'P')
write.table(one_way_anova_2, file = file.path(outpath, 'One_way_ANOVA.txt'), sep = '\t',
            row.names = FALSE, quote = FALSE, na = '')

TukeyHSD(anova_res)
# plot(TukeyHSD(anova_res, conf.level = 0.95),las=1, col = "red")

# T-Test ------------------------------------------------------------------

ttest_list <- NULL

for (i in seq(length(fenzu)-1)) {
  for (j in (i + 1):length(fenzu)) {
    cat(i, j, '\n')
    input_data_ij <- subset(input_data, Type %in% c(fenzu[i], fenzu[j]))
    group_1 <- fenzu[i]
    group_2 <- fenzu[j]
    for (alpha in c('chao1', 'ace', 'shannon', 'simpson', 'goods_coverage')) {
      fml <- formula(paste(alpha,"~ Type"))
      tt_res_ij <- t.test(fml, data = input_data_ij, 
                                   alternative = 'two.sided', paired = FALSE)
      ttest_list <- rbind(ttest_list, c(
        alpha,
        paste(group_1, group_2, sep = '/'),
        tt_res_ij$conf.int[[1]],
        tt_res_ij$conf.int[[2]],
        tt_res_ij$estimate[[1]],
        tt_res_ij$estimate[[2]],
        tt_res_ij$p.value
      ))
    }
    
  }
}
ttest_list <- data.frame(ttest_list)
names(ttest_list) <- c('Type', 'Groups', 'conf_down', 'conf_up', 'estimate_group0', 'estimate_group1', 'P')
write.table(wilcoxon_list, file = file.path(outpath, 'Student_ttest.txt'), sep = '\t',
            row.names = FALSE, quote = FALSE, na = '')

# Friedman test -----------------------------------------------------------
# Friedman检验更适合非独立样本



# use ggpubr to make question simple --------------------------------------

for (i in seq(length(fenzu)-1)) {
  for (j in (i + 1):length(fenzu)) {
    input_data_ij <- subset(input_data, Type %in% c(fenzu[i], fenzu[j]))
    for (alpha in c('chao1')) {
      fml <- formula(paste(alpha,"~ Type"))
      compare_means(fml, data = input_data_ij, method = "wilcox.test", paired = FALSE)
    }
  }
}

compare_means(chao1~Type, data = input_data, method = "wilcox.test")

#  以上两个结果中p.value矫正后的结果不一致，可知后者是对各分组的holm矫正，但是分组很少，矫正并不需要





