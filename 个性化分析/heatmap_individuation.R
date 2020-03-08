#! /urs/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
infile <- args[1]
mapfile <- args[2]
outpath <- args[3]

if (!require(tidyverse)){
  install.packages('tidyverse')
} else {
  library(tidyverse)
}

library(pheatmap)
require(RColorBrewer)
library(ComplexHeatmap)
library(ggsci)

infile_path <- '/Users/congliu/prs/nec_rice_analysis/BPD_noBPD/combined_otu_table_m3_std.txt'
mapfile_path <- '/Users/congliu/prs/nec_rice_analysis/BPD_noBPD/map.txt'
outpath <- '/Users/congliu/prs/nec_rice_analysis/BPD_noBPD'
one_group <- c('D4_BPD', 'D4_noBPD')
group1 <- one_group[1]
group2 <- one_group[2]

save_pheatmap_png <- function(x, filename, width=12, height=10, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf <- function(x, filename, width=12, height=12) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

taxnomy = c('s__', 'g__', 'f__', 'o__', 'c__', 'p__')
m3 <- read_delim(infile_path, delim = '\t')
map_file <- read.delim(mapfile_path, sep = '\t') %>% 
  select(c('Type', 'Description'))
mapfile_2 <- mapfile %>% separate(Type, c('Day', 'State'), sep = '_', remove = F)

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

heat <- function(tax){
  m3 <- m3 %>% filter(str_detect(`OTU ID`, tax)) %>% 
    select(-'OTU ID') %>%
    select(taxonomy, everything()) %>% 
    mutate(total=rowSums(.[2:(length(m3)-1)])) %>% 
    arrange(desc(total)) %>% 
    select(taxonomy, total, everything()) %>% 
    slice(1:20)
  
  m3_input <- as.data.frame(m3)
  rownames(m3_input) <- m3_input[,1]
  m3_input <- m3_input[-c(1,2)]
  m3_input <- log(m3_input + 1)
  m3_bray <- vegan::vegdist(t(m3_input), method = "bray")
  
  
  annotation_col <- data.frame(
    Type = map_file[,'Type']
    #Description = df[, 'Description']
  )
  rownames(annotation_col) <- map_file$Description
  
  #annotation_row <- data.frame()
  
  #m3_input <- t(m3_input)
  p <- pheatmap(m3_input,cluster_cols =T, cluster_rows = T, 
                clustering_distance_rows = vegan::vegdist(m3_input, method = "bray"),
                clustering_distance_cols = vegan::vegdist(t(m3_input), method = "bray"),
                 #color = colorRampPalette(colors = c("blue","yellow","red"))(100),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           #display_numbers = matrix(ifelse(m3_input > 5000, "*", ""), nrow(m3_input)), number_color = "black",
           #cellwidth = 15, cellheight = 15, main = "Example heatmap", 
           #scale = 'row', 
           clustering_method = 'average',
           distfun=dist_fun,
           #border = True, border_color = 'black',
           annotation_col = annotation_col, fontsize_col=4, fontsize_row=8,
           annotation_legend = TRUE
           #width = 8, height = 8
           )
  # dist_fun = function(x) dist(x, method="euclidean")
  dist_fun <- function(x) vegan::vegdist(x, method = 'bray')
  hclust_fun = function(x) hclust(x, method="average")
  pp <- gplots::heatmap.2( as.matrix(m3_input), trace="none", dendrogram="both", 
                   Rowv=TRUE, Colv=TRUE, distfun=dist_fun, hclustfun=hclust_fun, 
                   col=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
  
  row_dend = as.dendrogram(hclust(dist(mat)))
  row_dend = color_branches(row_dend, k = 2)
  ppp <- Heatmap(m3_input, cluster_rows = TRUE, cluster_columns = TRUE, name = '',
                 column_title = 'HeatMap', col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                 clustering_distance_rows = function(m) vegan::vegdist(m, method = 'bray'),
                 clustering_distance_columns = function(m) vegan::vegdist(m, method = 'bray'),
                 clustering_method_rows = 'average',
                 clustering_method_columns = 'average',
                 #cluster_rows = row_dend,
                 row_dend_reorder = TRUE,
                 row_names_gp = gpar(fontsize = 8),
                 #row_split =,
                 top_annotation = HeatmapAnnotation(df = annotation_col)
                 )
  png(file.path(outpath, paste('heatmap_', 's__', '.png', sep = '')), 
      width = 1000, height = 800)
  draw(ppp)
  dev.off()
  ggsave(p, filename = file.path(outpath, paste('heatmap_', 's__', '.png', sep = '')), 
         width = width, height = 10, dpi = 800)
  ggsave(p, filename = file.path(outpath, paste('heatmap_', 's__', '.pdf', sep = '')), width = width, height = 10)
  #save_pheatmap_pdf(pp, paste('/Users/congliu/prs/many_reporter/JSSRM02/', tax, '.pdf', sep = ''))
}

for (tax in taxnomy) {
  print(tax)
  heat(tax)
}
#heat('s__')

