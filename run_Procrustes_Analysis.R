#! /urs/bin/env Rscript

library(tidyverse)
library(vegan)

#' Title
#'
#' @param otu otu table
#' @param env env table, same sample as otu
#' @param group sample group information
#' @param method method choose to perform for otu table, must be one of (PCA PCoA NMDS)
#'
#' @return a ggplot obj
#' @export
#'
#' @examples
#' Rscript run_Procrustes_Analysis.R
#' @import vegan
#' @importFrom tidyverse ggplot2
run_Procrustes <- function(otu, env, group, method="PCA"){
  if (method == "PCA"){
    env_pca <- rda(env, scale = TRUE)
    otu_hell <- decostand(otu, method = 'hellinger')
    otu_pca <- rda(otu_hell, scale = FALSE)
    site_env <- summary(env_pca, scaling = 1)$site
    site_otu <- summary(otu_pca, scaling = 1)$site
  } else if (method == "PCoA") {
    otu.bray <- vegdist(otu, method = 'bray')
    add <-  !(is.euclid(otu.bray))
    #env_pcoa <- cmdscale(env)
    # 对于环境因素一般采用pca，且其量纲不一致需要scale
    env_pca <- rda(env, scale = TRUE)
    otu_pcoa <- cmdscale(otu, k = nrow(otu)-1, eig = TRUE, add = add)
    site_env <- summary(env_pca, scaling = 1)$site
    site_otu <- data.frame({otu_pcoa$point})[1:2]
  } else if (method == "NMDS"){
    env_nmds <- metaMDS(env, distance = 'bray', k=2)
    otu_nmds <- metaMDS(otu, distance = 'bray', k=2)
    site_env <- data.frame(env_nmds$point)
    site_otu <- data.frame(otu_nmds$point)
  }
  proc <- procrustes(X=site_env, Y=site_otu, symmetric = TRUE)
  prot <- protest(X = env_pca, Y = otu_pca, permutations = how(nperm = 999))
  signif <- prot$signif
  m2 <- prot$ss
  
  Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
  X <- data.frame(proc$rotation)
  Y$sample <- rownames(Y)
  Y <- merge(Y, group, by.x = 'sample', by.y = 'samples')
  p <- ggplot(data = Y) +
    geom_point(aes(x = X1, y = X2, colour = groups), size = 1.5, shape = 16) +
    geom_point(aes(x = PC1, y = PC2, color = groups), size = 1.5, shape = 2) +
    geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.05, 'cm'))) +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    labs(x = 'Dimension 1', y = 'Dimension 2') +
    geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
    geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
    geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
    geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
    annotate('text', label = sprintf(paste('M^2 == ', signif)),
             x = -0.21, y = 0.42, size = 3, parse = TRUE) +
    annotate('text', label = paste('P < ', round(m2, 3)),
             x = -0.21, y = 0.38, size = 3, parse = TRUE)
}