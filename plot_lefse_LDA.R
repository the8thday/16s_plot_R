#! /urs/bin/env Rscript
# plot lefse LDA 

#' Title
#'
#' @param lefse.result a .res file created by run_lefse.py
#' @param negate.class 
#' @param lda.threshold 
#' @param group 
#'
#' @return
#' @export
#'
#' @examples
lefse.lda.plot <- function(lefse.result, negate.class=NULL, lda.threshold=NULL, group=NULL){
  input <- read_tsv(lefse.result, col_names = c('feature','dummy','class','lda','pvalue'))
  
  tmp <- filter(input, !is.na(class)) %>% mutate(lda.scale=lda)
  if(!is.null(lda.threshold)){
    tmp <- filter(tmp, lda>lda.threshold)
  }
  if(!is.null(negate.class) && intersect(negate.class, tmp$class)==negate.class){
    tmp <- mutate(tmp, lda=ifelse(class %in% c(negate.class), lda*-1, lda)) %>%
      mutate(class=factor(class, levels=c(negate.class, setdiff(class, negate.class))))
  }
  tmp$group.origin <- ''
  if(!is.null(group)){
    colnames(group) <- c("X1","group")
    group$X1 <- as.character(group$X1) ## in case imported as factor
    tmp <- merge(tmp, group, by=1, all.x=TRUE)
    tmp$group[is.na(tmp$group)] <- 'Others'
    tmp$group.origin <- tmp$group
  }
  tmp <- mutate(tmp, group=paste0(class, group))
  
  plot.dat <- dplyr::arrange(tmp, class, lda) %>%
    mutate(group=factor(group, levels=rev(unique(group)), ordered = TRUE)) %>%
    mutate(feature=factor(feature, levels=feature, ordered = TRUE))
  
  tmp <- select(plot.dat, group, group.origin) %>% unique()
  anno <- tmp$group.origin
  names(anno) <- tmp$group
  
  p <- ggplot(plot.dat, aes(x=feature, y=lda, fill=class, label=feature)) +
    geom_bar(stat='identity') +
    geom_text(aes(y=0, x=feature), hjust=ifelse(plot.dat$lda<0, -0, 1), nudge_y = -sign(plot.dat$lda)*0.1) +
    coord_flip() +
    labs(x=NULL) +
    theme_classic() +
    facet_grid(group~., scale='free_y',  space = "free_y",
               labeller=labeller(group=anno)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'cm'), legend.position = 'left',
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          panel.spacing = unit(0, "lines"))
  
  if(!is.null(group)){
    p <- p + theme(strip.text.y = element_text(angle=0, margin = margin(0,3,0,3, "cm")))
  }
  p
}

lefse.lda.plot('/Users/congliu/prs/ZJDXGW/lefse_input_A_B.res', negate.class = 'A', lda.threshold = 2)