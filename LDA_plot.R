#! /urs/bin/env Rscript

library(MASS)
ord <- lda(Species ~ ., iris, prior = rep(1, 3)/3)
library(MASS)
ord <- lda(Species ~ ., iris, prior = rep(1, 3)/3)

library(yyplot)
p + geom_ord_ellipse(ellipse_pro = .96, color='firebrick', size=1, lty=3) +
  geom_ord_ellipse(ellipse_pro = .99, lty=2) 