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


taxnomy = c('s__', 'g__', 'f__', 'o__', 'c__', 'p__')
m3 <- read_delim(infile, delim = '\t')
map_file <- read.delim(mapfile, sep = '\t') %>% 
  select(c('Type', 'Description'))

sample_num <- length(m3)-2
width <- round(0.1 * sample_num)

if (width < 10) {
  width <- 10
} else if(width > 30){
  width <- 30
} else {
  width <- width
}
cat("width: ", width)

save_pheatmap_pdf <- function(x, filename, width=width, height=15) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png <- function(x, filename, width=width, height=10, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

heat <- function(tax){
  m3 <- m3 %>% filter(str_detect(`OTU ID`, tax)) %>% select(-'OTU ID') %>%
    select(taxonomy, everything()) %>% 
    mutate(total=rowSums(.[2:(length(m3)-1)])) %>% arrange(desc(total)) %>% 
    select(taxonomy, total, everything()) %>% 
    slice(1:50)
  
  m3_input <- as.data.frame(m3)
  rownames(m3_input) <- m3_input[,1]
  m3_input <- m3_input[-c(1,2)]
  m3_input <- log(m3_input + 1)
  
  annotation_col <- data.frame(
    Type = map_file[,'Type']
    #Description = df[, 'Description']
  )
  rownames(annotation_col) <- map_file$Description
  
  #annotation_row <- data.frame()
  
  #m3_input <- t(m3_input)
  pp <- pheatmap(m3_input,cluster_cols =T, cluster_rows = T, 
                 clustering_distance_rows = vegan::vegdist(m3_input, method = "bray"),
                 clustering_distance_cols = vegan::vegdist(t(m3_input), method = "bray"),
                 #color = colorRampPalette(colors = c("blue","yellow","red"))(100),
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                 #display_numbers = matrix(ifelse(m3_input > 5000, "*", ""), nrow(m3_input)), number_color = "black",
                 #cellwidth = 15, cellheight = 15, main = "Example heatmap", 
                 cellheight = 10,
                 #scale = 'row', 
                 clustering_method = 'average',
                 #border = True, border_color = 'black',
                 annotation_col = annotation_col, fontsize_col=4, fontsize_row=8,
                 annotation_legend = TRUE
                 #width = 8, height = 8
                 )
  
  #ggsave(pp, filename = file.path(outpath, paste('heatmap_', tax, '.png', sep = '')), width = width, height = 10, dpi = 800)
  ggsave(pp, filename = file.path(outpath, paste('heatmap_', tax, '.pdf', sep = '')), width = width, height = 10)
  #save_pheatmap_pdf(pp, file.path(outpath, paste('heatmap_', tax, '.pdf', sep = '')))
}

for (tax in taxnomy) {
  print(tax)
  heat(tax)
}
#heat('s__')


