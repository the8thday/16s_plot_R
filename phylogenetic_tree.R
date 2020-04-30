#! /urs/bin/env Rscript

library(phyloseq)
library(tidyverse)
library(microbiomeViz)
# tidy data ---------------------------------------------------------------
# 似乎biom格式并不是很一致，先注释
# rich_dense = import_biom('/Users/congliu/prs/ZJDXGW/otu_table.biom', parseFunction=parse_taxonomy_default)
# #rich_dense = import_biom(rich_dense_biom, parseFunction=parse_taxonomy_greengenes)
# tr = parsePhyloseq(rich_dense)
# p = tree.backbone(tr, size=1)
# p

# use m2 file -------------------------------------------------------------

infile <- '/Users/congliu/prs/ruijin_xirou/combined_otu_table_m2_std.txt'
map_file <- '/Users/congliu/prs/ruijin_xirou/map.txt'
outpath <- '/Users/congliu/prs/ruijin_xirou/'


otu <- read.delim(infile, row.names = 1, sep = '\t', 
                  stringsAsFactors = FALSE, check.names = FALSE) 
otumat <- subset(otu, taxonomy != 'Unassigned', select = -taxonomy)
mapfile <- read.delim(map_file, sep = '\t') %>% dplyr::select(Type, Description)

taxmat <- otu %>% dplyr::select(c(taxonomy)) %>% 
  filter(taxonomy != 'Unassigned') %>% 
  separate(col = taxonomy, 
           into = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'),
           sep = "; ")
rownames(taxmat) <- rownames(otumat)
taxmat <- as.matrix(taxmat)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)
physeq
GP = transform_sample_counts(physeq, function(otu) otu/sum(otu))
GP = filter_taxa(GP, function(x) max(x)>=0.01,TRUE)
GP = fix_duplicate_tax(GP)

tr = parsePhyloseq(GP)
p = tree.backbone(tr, size=1)
p
anno.data <- data.frame(node=c("g__Roseburia", "c__Clostridia", "s__Bacteroides_ovatus"),
                        color='red', stringsAsFactors = FALSE)
p <- clade.anno(p, anno.data)
p
