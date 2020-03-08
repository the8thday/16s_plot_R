#画附一图
library(tidyverse)
df <- read_delim('/Users/congliu/prs/fenxi/2019.11.19_苏附一徐俊项目_diff/High__vs__Low_meta_diff.txt', delim = '\t')
df_f <- df %>% filter(pvalues < 0.05) %>% arrange(desc(`+samples in group 1`, ))
df_20 <- select(df_f[1:20, ], c('OTU', '+samples in group 1'))

tiqu <- function(path){
  df <- read_delim(path, delim = '\t')
  df_f <- df %>% filter(pvalues < 0.05) %>% arrange(desc(`+samples in group 1`, ))
  df_20 <- select(df_f[1:20, ], c('OTU', '+samples in group 1'))
  write.table(df_20, 'all_filter.txt', append = T, row.names=F)
}

setwd('/Users/congliu/prs/fenxi/2019.11.19_苏附一徐俊项目_diff/')
paths <- list.files('/Users/congliu/prs/fenxi/2019.11.19_苏附一徐俊项目_diff/', pattern = '^[A-Z]')

for (i in paths) {
  print(i)
  tiqu(i)
}

foo <- read_delim('/Users/congliu/prs/fenxi/sfy/combined_otu_table_m2_std.txt', delim = '\t')

mapbiao <- read_delim('/Users/congliu/prs/fenxi/sfy/map_s.txt', delim = '\t')
h_s1 <- read_delim(paste('/Users/congliu/prs/fenxi/2019.11.19_苏附一徐俊项目_diff/', 'xag', sep = ''), delim = ' ')
bar <- t(as.data.frame(inner_join(foo, h_s1, by = c('OTU ID' = 'OTU'))))
colnames(bar) <- bar[1, ]
bar <- bar[-1, ]
bar <- as_tibble(bar, rownames = 'OTU')
mapbiao <- select(mapbiao, c('#Sample', 'Type5'))
#mapbiao <- select(mapbiao, c('#Sample', type))
hehe <- inner_join(bar, mapbiao, by = c('OTU' = '#Sample'))
#haha <- mutate(hehe, Type=Type5) %>% select(-c('OTU', 'Type5'))
haha <- hehe %>% filter(.[['Type5']] %in% c('MDA', 'MUC'))
haha <- hehe %>% filter(UQ(as.symbol('Type5')) %in% c('MDA', 'MUC'))
#此处有问题
#haha <- hehe[hehe['Type5'] == 'MDA']
#haha <- filter(hehe, Type5 %in% c('MDA','MUC')) %>% select(-c('OTU'))
p <- gather(haha, Sample, value, -Type)
p <- mutate(p, value = as.double(value))
#colnames(p) <- c('Type','Sample','value')
ggplot(data = p, aes(x=Sample, y=value, fill=Type)) + geom_bar(stat = "identity",position = "dodge")+
  ylab('Relative Abundance') + coord_flip() + theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + xlab('')

pppp <- function(type, x, a, b){
  name <- paste(a,b, sep = '_')
  mapbiao <- read_delim('/Users/congliu/prs/fenxi/sfy/map_s.txt', delim = '\t')
  h_s1 <- read_delim(paste('/Users/congliu/prs/fenxi/2019.11.19_苏附一徐俊项目_diff/', x, sep = ''), delim = ' ')
  bar <- t(as.data.frame(inner_join(foo, h_s1, by = c('OTU ID' = 'OTU'))))
  colnames(bar) <- bar[1, ]
  bar <- bar[-1, ]
  bar <- as_tibble(bar, rownames = 'OTU')
  #mapbiao <- select(mapbiao, c('#Sample', 'Type3'))
  mapbiao <- select(mapbiao, c('#Sample', type))
  hehe <- inner_join(bar, mapbiao, by = c('OTU' = '#Sample'))
  haha <- hehe %>% filter(UQ(as.symbol(type)) %in% c(a, b)) %>% select(-c('OTU'))
  #haha <- filter(hehe, Type5 %in% c('MDA','MUC')) %>% select(-c('OTU'))
  p <- gather(haha, Sample, value, -type)
  p <- mutate(p, value = as.double(value))
  colnames(p) <- c('Type','Sample','value')
  pp <- ggplot(data = p, aes(x=Sample, y=value, fill=Type)) + geom_bar(stat = "identity",position = "dodge")+
    ylab('Relative Abundance') + coord_flip() + theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + xlab('')
  ggsave(pp, filename = paste('/Users/congliu/prs/fenxi/2019.11.19_苏附一徐俊项目_diff/sfy/',name, '.pdf', sep = ''), width=8, height=8)
}
    
pppp('Type3', 'xaa', 'Healthy', 'Stage_I')
pppp('Type3', 'xab', 'Healthy', 'Stage_II')
pppp('Type3', 'xac', 'Healthy', 'Stage_III')
pppp('Type3', 'xad', 'Healthy', 'Stage_IV')

pppp('Type1', 'xae', 'Low', 'High')
pppp('Type2', 'xaf', 'Left', 'Right')

pppp('Type5', 'xag', 'MDA', 'MUC')
pppp('Type4', 'xah', 'Stable', 'Unstable')
pppp('Type3', 'xai', 'Stage_I', 'Stage_II')
pppp('Type3', 'xaj', 'Stage_I', 'Stage_III')
pppp('Type3', 'xak', 'Stage_I', 'Stage_IV')
pppp('Type3', 'xal', 'Stage_II', 'Stage_III')
pppp('Type3', 'xam', 'Stage_II', 'Stage_IV')
pppp('Type3', 'xan', 'Stage_III', 'Stage_IV')




