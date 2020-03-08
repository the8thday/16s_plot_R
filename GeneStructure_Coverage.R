#!/usr/bin/env Rscript

library("getopt")
command=matrix(c("gtf",   "g", 1, "character",
                   "gtf file of gene structure annotation, [REQUIRED]",
                 "depth", "d", 1, "character",
                   "depth file from 'samtools depth [...]', [REQUIRED]",
                 "outdir","o", 1, "character",
                   "output dir of plot results, [REQUIRED]",
                 "genelst", "l", 2, "character",
                   "genelst file of wanted plot, optional argument",
                 "help",  "h", 0, "logical",
                   "print help infor and exit"), byrow=T, ncol=5)
args=getopt(command)
if (!is.null(args$help) || is.null(args$gtf) || is.null(args$depth) || is.null(args$outdir)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}
gtf <- args$gtf
depth <- args$depth
outdir <- args$outdir

library(ggplot2)
library(patchwork)

Annot <- read.table(gtf, header=F, sep='\t', quote="")
Annot$gene <- sapply(
  as.character(Annot[,9]),
  function(x) sub(".*gene_name\\s\"([^;]+)\";.*", "\\1", x)
)
Annot$transcript <- sapply(
  as.character(Annot[,9]),
  function(x) sub(".*transcript_name\\s\"([^;]+)\";.*", "\\1", x)
)
Annot <- Annot[,c(1, 3, 4, 5, 7, 10, 11)]
colnames(Annot) <- c('chro', 'feature', 'start', 'end', 'strand', 'gene', 'transcript')

Depth <- read.table(depth, header=F, sep='\t', quote="")
colnames(Depth) <- c('chro', 'posi', 'depth')

genelst <- unique(Annot$gene)
if (!is.null(args$genelst)){
  tmp <- unique(read.table(args$genelst, header=F, sep='\t', quote="")$V1)
  genelst <- tmp[tmp %in% genelst]
}
setwd(outdir)

for (geneid in genelst){
  tmp <- Annot[(Annot$gene == geneid),]
  Xchr <- as.character(unique(tmp$chro))
  Xstart <- min(tmp$start)
  Xend <- max(tmp$end)
  length_gene <- Xend-Xstart
  border <- if (length_gene < 10000) length_gene/35 else 300

  gene <- tmp[(tmp$feature) == 'transcript',]
  exon <- tmp[(tmp$feature) == 'exon',]
  cds <- tmp[(tmp$feature) == 'CDS',]
  p1 <- ggplot() +
    geom_segment(data=gene, aes(x=start, xend=end, y=transcript, yend=transcript),
                 size=0.7, color='#000000') +
    geom_segment(data=exon, aes(x=start, xend=end, y=transcript, yend=transcript),
                 size=3.5, color='#1800ff') +
    geom_segment(data=cds, aes(x=start, xend=end, y=transcript, yend=transcript),
                 size=3.5, color='#ff9800') +
    xlim(c(Xstart-border, Xend+border)) +
    ylab("") + xlab("") + theme_bw()
    # theme_classic() +
  # p1 <- ggplot() +
  #   geom_segment(data=gene, aes(x=start, xend=end, y=0, yend=0), size=1, color='#000000') +
  #   geom_rect(data=exon, aes(xmin=start, xmax=end, ymin=-.07, ymax=.07), fill='#1800ff') +
  #   geom_rect(data=cds, aes(xmin=start, xmax=end, ymin=-.07, ymax=.07), fill='#ff9800') +
  #   xlim(c(Xstart-border, Xend+border)) +
  #   ylim(c(-.5, .5))

  tmp2 <- Depth[(Depth$chro == Xchr) & (Depth$posi >= Xstart-border) & (Depth$posi <= Xend+border),]
  if (dim(tmp2)[1] == 0) {
    cat(paste('There is not reads target ot Gene: ', geneid, ', skipping...\n', sep=''))
    next
  }
  p2 <- ggplot(tmp2) +
    geom_bar(aes(x=posi, y=depth),
             width=1, stat="identity", color='#f44336') +
    # geom_histogram(aes(x=pos, y=depth), stat="identity", color='#4ec34a') +
    xlim(c(Xstart-border, Xend+border)) +
    ggtitle(geneid) +
    ylab("Depth Coverage") + xlab("") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

  number_transcript <- length(unique(tmp$transcript))
  hight <- if (number_transcript > 3) 1 else 2
  if (number_transcript > 5){
    hight <- 0.6
  }
  p <- p2 + p1 + plot_layout(ncol = 1, heights = c(hight, 1))
  ggsave(paste(geneid, '.GeneCoverage.png', sep=''), plot=p, height=4, width=6, type='cairo-png')
  # ggsave(paste(geneid, '.GeneCoverage.pdf', sep=''), plot=p, height=4, width=6)
}
