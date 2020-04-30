#! /urs/bin/env Rscript
# a LDA analysis

library(MASS)
library(ggord)
library(yyplot)
ord <- lda(Species ~ ., iris, prior = rep(1, 3)/3)
ggplot(cbind(iris, predict(ord)$x), aes(LD1, LD2, color = Species)) +
  geom_point() +
  stat_ellipse(level = 0.95, show.legend = FALSE)


p + geom_ord_ellipse(ellipse_pro = .96, color='firebrick', size=1, lty=3) +
  geom_ord_ellipse(ellipse_pro = .99, lty=2) 


ord <- prcomp(iris[, 1:4])
ggord(ord, iris$Species) + geom_ord_ellipse(lty=2)


ord <- metaMDS(iris[, 1:4])
ggord(ord, iris$Species) + geom_ord_ellipse(lty=2)

# plot lefse LDA--------------------------------------------------------------

lefse.lda.plot('.res', )
