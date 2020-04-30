#! /urs/bin/env Rscript
# 画alpha多样性差异的图片

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4)
infile <- args[1]
mapfile <- args[2]
outpath <- args[3]
compare_groups <- args[4]
sorted_group <- args[5]


library(ggpubr)
#library(reshape2)
library(patchwork)
library(tidyverse)
# infile <- '/Users/congliu/prs/many_reporter/SDFY02_map/combined_alpha_m1.txt'
# mapfile <- '/Users/congliu/prs/many_reporter/SDFY02_map/map.txt'
# outpath <- '/Users/congliu/prs/many_reporter/SDFY02_map/'
# compare_groups <- "list(c('adenoma','normal'),c('cancer','normal'),c('inflammation','normal'))"
# print(compare_groups)
# sorted_group <- "c('adenoma','cancer','inflammation','normal')"

alpha <- read_delim(infile, delim = '\t')
mapfile <- read_delim(mapfile, delim = '\t') %>% select(Type, Description)
dt <- right_join(mapfile, alpha, by = c('Description' = 'Sample')) %>% select(-Description)
##分组信息

#dt1 <- melt(dt,id.vars="Type",variable.name="OTUID",value.name="abundance")
dt1 <- pivot_longer(dt, cols = -Type, names_to = 'OTUID', values_to = 'abundance')
dt1 <- dt1 %>% filter(Type %in% eval(parse(text = (sorted_group))))
#chao1 boxplot
#df1 <-dt1[dt1$OTUID %in% c("chao1"),]
df1 <- dt1 %>% filter(OTUID %in% c("chao1"))
my_comparisons <- eval(parse(text = (compare_groups)))
print(my_comparisons)
df1$Type=factor(df1$Type, levels = eval(parse(text = (sorted_group))))
p1 <- ggboxplot(df1, x = "Type", y = "abundance",color = "Type", palette = "jco", add = c("mean","jitter"), 
          legend="right", error.plot = "linerange") +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise crossbar p-value
  stat_compare_means(label.y = 1700000)+     # Add global p-value
  xlab("")+
  ylab("chao1")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#ace boxplot
df2 <-dt1[dt1$OTUID %in% c("ace"),]
df2$Type=factor(df2$Type, levels = eval(parse(text = (sorted_group))))
p2<-ggboxplot(df2, x = "Type", y = "abundance",color = "Type", palette = "jco", add = c("mean","jitter"), 
              legend="right", error.plot = "linerange") +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise crossbar p-value
  stat_compare_means(label.y = 80000)+     # Add global p-value
  xlab("")+
  ylab("ace") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

#shannon boxplot
df3 <-dt1[dt1$OTUID %in% c("shannon"),]  #ɸѡotuid
df3$Type=factor(df3$Type, levels = eval(parse(text = (sorted_group))))
p3<-ggboxplot(df3, x = "Type", y = "abundance",color = "Type", palette = "jco", add = c("mean","jitter"), 
              legend="right", error.plot = "linerange") +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise crossbar p-value
  stat_compare_means(label.y = 20)+     # Add global p-value
   xlab("")+
  ylab("shannon") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#simpson	 boxplot
df4 <-dt1[dt1$OTUID %in% c("simpson"),]  #ɸѡotuid
df4$Type=factor(df4$Type, levels = eval(parse(text = (sorted_group))))
p4<-ggboxplot(df4, x = "Type", y = "abundance",color = "Type", palette = "jco", add = c("mean","jitter"), 
              legend="right", error.plot = "linerange") +
  stat_compare_means(comparisons = my_comparisons) + # Add pairwise crossbar p-value
  stat_compare_means(label.y = 1.015)+     # Add global p-value
   xlab("")+
  ylab("simpson") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p5<-(p1+p2)/(p3+p4)

p6 <- p5 + plot_annotation(title="Alpha diversity",theme = theme(plot.title = element_text(size = 18,hjust = 0.5)))
pdf(file.path(outpath, "Alpha_diversity.pdf"),width=15,height=15)
ggsave(file.path(outpath, "Alpha_diversity.png"),width=15,height=15, dpi = 800)
p6
dev.off()
