## learn to use phyloseq

library(tidyverse)
library(phyloseq)
# library(ggtree)
library(microbiomeViz)
# library(qiime2R)


# change data to phyloseq -------------------------------------------------

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
# physeq对象包含样本和


# R包microbiomeViz ---------------------------------------------------------
data("GlobalPatterns")
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

# qiime -------------------------------------------------------------------

qiimedata = import_qiime(otufile, mapfile, trefile, rs_file)

# ggtree ------------------------------------------------------------------

GP <- prune_taxa(taxa_sums(physeq) > 0, physeq)
GP.chl <- subset_taxa(GP, kingdom=="k__Bacteria")

ggtree(GP.chl) + geom_text2(aes(subset=!isTip, label=label), hjust=-.2, size=4) +
  geom_tiplab(aes(label=Genus), hjust=-.3) +
  geom_point(aes(x=x+hjust, color=SampleType, shape=Family, size=Abundance),na.rm=TRUE) +
  scale_size_continuous(trans=log_trans(5)) +
  theme(legend.position="right")

# plot lefse cladogram & LDA----------------------------------------------------

lefse_in <- read.table("/Users/congliu/prs/ZJDXGW/lefse_input_D_A.txt", head=TRUE, sep = '\t',
                       stringsAsFactors = FALSE)
# lefse_in <- read.table("/Users/congliu/prs_R/profiled_samples/merged_abundance_table.txt", head=TRUE, sep = '\t',
#                        stringsAsFactors = FALSE)
colnames(lefse_in) <- c('ID', LETTERS[1:length(lefse_in)-1])
#lefse_in <- lefse_in[,-ncol(lefse_in)]
## Use row means as a proxy for node size
dat <- data.frame(V1=lefse_in[,1], V2=rowMeans(lefse_in[,-1]), stringsAsFactors = FALSE)
#dat <- distinct(dat, V1, .keep_all = TRUE)
tr <- parseMetaphlanTSV(dat, node.size.offset=2, node.size.scale=0.8)
p <- tree.backbone(tr, size=0.5)
lefse_lists = data.frame(node=c('s__Haemophilus_parainfluenzae','p__Proteobacteria',
                                'f__Veillonellaceae','o__Selenomonadales',
                                'c__Negativicutes', 's__Streptococcus_parasanguinis',
                                'p__Firmicutes','f__Streptococcaceae',
                                'g__Streptococcus','o__Lactobacillales',
                                'c__Bacilli','s__Streptococcus_mitis'),
                         color=c(rep('darkgreen',6), rep('red','6')),
                         stringsAsFactors = FALSE
)
p <- clade.anno(p, lefse_lists, alpha=0.3)
p

# something about phyloseq ------------------------------------------------

# Unifrac距离
# 加权Unifrac距离
wei_unif_dis <- distance(physeq, method = 'wunifrac')
# 非加权 Unifrac 距离, only considers the presence/absence of taxa between sample pairs
# unwei_unif_dis <- distance(physeq, method = 'unifrac')
unwei_unifrac_dis <- UniFrac(physeq = physeq, weighted = FALSE)
# Bray-curtis 距离
bray_dis <- distance(physeq, method = 'bray')

















