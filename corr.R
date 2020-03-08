library(tidyverse)
library(Hmisc)
library(readxl)
library(ggcorr)


js <- read_excel('/Users/congliu/Downloads/cong.xlsx')

js %>% select(-`#OTU ID`) %>% cor()

js2 <- js %>% select(-`#OTU ID`)
cortest <- rcorr(as.matrix(js2), type = 'pearson')
corr.r <- as.data.frame(cortest$r)
corr.p <- as.data.frame(cortest$P)
corr.r[corr.p <= 0.05 & abs(corr.r) >= 0.6]


library(psych)
library(corrplot)

res <- corr.test(js2, use = "complete", method = "pearson", adjust = "fdr")
