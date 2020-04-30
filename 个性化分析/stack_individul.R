 #! /urs/bin/env Rscript
##堆积图
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
infile_path <- args[1]
mapfile_path <- args[2]
outpath <- args[3]

if (!require(tidyverse)){
  install.packages('tidyverse')
} else {
  library(tidyverse)
}
#library(patchwork)

infile_path <- '/Users/congliu/prs/hbrm_Alveoli/combined_otu_table_m3_std.txt'
mapfile_path <- '/Users/congliu/prs/hbrm_Alveoli/map.txt'
outpath <- '/Users/congliu/prs/hbrm_Alveoli'
one_group <- c('lung_cancer', 'non_lung_cancer')
group1 <- one_group[1]
group2 <- one_group[2]

m3 <- read_delim(infile_path, delim = '\t')
mapfile <- read_delim(mapfile_path, delim = '\t') %>% select(Type, Description)
#mapfile_2 <- mapfile %>% separate(Type, c('Day', 'State'), sep = '_', remove = F)

taxnomy = c('s__', 'g__', 'f__', 'o__', 'c__', 'p__')
#Color = c("yellow", "darkgreen", "blue4", "lightblue", "maroon3", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
Color <- c('darkgreen', 'blue', 'green', 'darkorange', 'yellow', 'hotpink', 'grey', 
           'skyblue', 'violet', 'cadetblue', '#0072B2')
sample_num <- length(m3)-2
width <- round(0.1 * sample_num)

if (width < 8) {
  width <- 8
} else if(width > 30){
  width <- 30
} else {
  width <- width
}
cat("width: ", width)
n<- 10
duiji <- function(tax, n=10){
  m3_s <- m3 %>% filter(str_detect(`OTU ID`, 'p__'))
  m3_df <- m3_s %>% mutate(total = rowSums(.[2:(length(m3_s)-1)])) %>% arrange(desc(total)) %>%
    select(taxonomy, total, everything()) %>% select(-c(`OTU ID`))
  m3_top10 <- m3_df %>% slice(1:n)
  cc <- m3_df %>% slice(n+1:n())
  m3_final <- rbind(m3_top10, c("Others", colSums(cc[,2:length(cc)])))
  m3_final <- m3_final %>% mutate_at(vars(colnames(m3_final)[2:length(colnames(m3_final))]), as.double)
  write.table(m3_final, file.path(outpath, 'top10_tax.txt'), sep = '\t')
  
  m3_data <- m3_final %>% pivot_longer(cols = -taxonomy, names_to = "variable", values_to = "value")
  m3_data <- dplyr::left_join(m3_data, mapfile, by=c('variable'='Description')) %>% 
    arrange(Type)
  m3_data <- within(m3_data, variable <- factor(variable, levels = unique(pull(m3_data, variable))))
  m3_data <- within(m3_data, taxonomy <- factor(taxonomy, levels = pull(m3_final, taxonomy)))
  m3_data <- within(m3_data, Day <- factor(Day, levels = (c('D4', 'D7', 'D14', 'D21', 'D28'))))
  pp <- ggplot(data = m3_data, aes(x=Day, y=value, fill=taxonomy)) +
    geom_bar(stat = 'identity', position = 'fill') +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90,vjust=0.5,size = 10), 
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          plot.background = element_blank())
    #facet_wrap(Type~.,scales="free")
    #facet_grid(~Type, scales="free", space = 'free_x') + 
    #scale_fill_manual(values = Color)
  ggsave(pp, filename = file.path(outpath, paste('stack_day_', tax, '.png', sep = '')), height = 8, width = 8, dpi=800)
  ggsave(pp, filename = file.path(outpath, paste('stack_day_', tax, '.pdf', sep = '')), height = 8, width = 8)
}

for (tax in taxnomy) {
  print(tax)
  duiji(tax)
}
duiji('p__')

##加上聚类树
m3_s <- m3 %>% filter(str_detect(`OTU ID`, 's__'))
dis_bray <- vegan::vegdist(t(m3_s), method = 'bray')
tree <- hclust(dis_bray, method = 'average')
plot(tree)
