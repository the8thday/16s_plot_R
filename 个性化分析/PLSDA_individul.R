#! /urs/bin/env Rscript


# 所需包 ---------------------------------------------------------------------

require(mixOmics)
require(tidyverse)
require(pROC)

# https://mixomicsteam.github.io/Bookdown/plsda.html#inputs-and-outputs
# tidy_data ---------------------------------------------------------------

infile <- '/Users/congliu/prs/lasso_data/crc2_4_45_46_47_48_49_mv34_qiime_silva_v2_m2_std.txt'
map_file <- '/Users/congliu/prs/lasso_data/map_20200108_add_ref_n_a_c.txt'
outpath <- '/Users/congliu/prs/nec_rice_analysis/brain_HL'
one_group <- c('c', 'n')


otu <- read.delim(infile, row.names = 1, sep = '\t', 
                  stringsAsFactors = FALSE, check.names = FALSE)
mapfile <- read.delim(map_file, sep = '\t') %>% dplyr::select(Type, Description)
otu <- subset(otu, taxonomy != 'Unassigned')
otu <- subset(otu, select = -taxonomy)
otu <- data.frame(t(otu))
#otu2 <- log(otu + 1)
otu2 <- otu
otu2$Description <- row.names(otu2)
otu3 <- merge(otu2, mapfile, by.x = 'Description', by.y = 'Description')
row.names(otu3) <- otu3$Description
otu3$Type <- as.character(otu3$Type)
otu3 <- subset(otu3, Type %in% one_group)
otu3$Type <- as.factor(otu3$Type)
otu4 <- subset(otu3, select = -c(Description, Type))

# PLSDA -------------------------------------------------------------------
# For PLS-DA, a loading >0.35 was chosen. feature selection, how could i choose

plsda_result <-plsda(otu4, otu3$Type, ncomp = 4)
#plotIndiv(plsda_result, ind.names = TRUE, style = 'ggplot2')
# plotVar(plsda_result, cutoff=0.5)
sample_site <- data.frame(plsda_result$variates)[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('plsda1', 'plsda2')
sample_site <- merge(sample_site, mapfile, by.x = 'names', by.y = 'Description', all.x = TRUE)
plsda_result_eig <- {plsda_result$explained_variance$X}[1:2]

p_plsda <- ggplot(sample_site, aes(plsda1, plsda2, color = Type, label = names)) +
  geom_point(size = 1.5, alpha = 0.6) + 
  stat_ellipse(show.legend = FALSE, linetype = 2) +
  scale_color_manual(values = c('#1D7ACC', '#F67433', '#00815F','#C673FF2E')) +
  scale_shape_manual(values = c(16, 17)) +
  theme(panel.grid = element_line(color = 'grey50'), 
        panel.background = element_rect(color = 'black', fill = 'transparent')) + 
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent')) +
  #labs(x = paste('PLS-DA axis1 ( explained variance ', round(100 * plsda_result_eig[1], 2), '% )', sep = ''), 
   #    y = paste('PLS-DA axis2 ( explained variance ', round(100 * plsda_result_eig[2], 2), '% )', sep = ''))
ggsave(p_plsda, filename = file.path(outpath, paste('plsda_', group1, '_', group2, '.pdf', sep = '')), height = 8, width = 8)

# go to further -----------------------------------------------------------

plotIndiv(plsda_result, ind.names = FALSE, legend = TRUE, 
          ellipse = TRUE, title = 'PLSDA')


splsda <- splsda(X=otu4, Y=otu3$Type, ncomp = 2, keepX=c(100,100))
selectVar(splsda, comp = 1)$value
plotLoadings(splsda, contrib = 'max', method = 'mean')

# tuning parameters ---------------------------------------从这里开始

plsda2 <- plsda(otu4, otu3$Type, ncomp = 10)
set.seed(42)
# 确定ncomp的个数
perf.plsda <- perf(plsda2, validation='Mfold', folds=10, 
                   progressBar=FALSE, nrepeat=50)
plot(perf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = 'horizontal')
perf.plsda$choice.ncomp
perf.plsda$error.rate

# 确定keepX的个数
list.keepX <- c(seq(10,300,20))
tune.splsda <- tune.splsda(otu4, otu3$Type, ncomp = 3,
                           validation = 'Mfold', 
                           folds = 10, 
                           progressBar = FALSE, 
                           measure = 'BER',
                           dist = 'max.dist',
                           test.keepX = list.keepX,
                           nrepeat = 50)
error <- tune.splsda$error.rate
ncomp <- tune.splsda$choice.ncomp$ncomp
select.keepX <- tune.splsda$choice.keepX[1:ncomp]
plot(tune.splsda, col = color.jet(3))

final.splsda <- splsda(otu4, otu3$Type, ncomp=ncomp, keepX=select.keepX)
plotIndiv(final.splsda, ind.names = FALSE, legend = TRUE, 
          comp = c(1,2),
          ellipse = TRUE, title = 'PLSDA comp 1 & 2')
auroc(final.splsda, roc.comp = 2) #仅是补充作用

select.var <- selectVar(final.splsda, comp = 1)$value
select.var <- select.var[order(select.var$value.var), colnames(select.var), drop=FALSE]
# select.var['OTU'] = rownames(select.var)
# select.var <- select.var %>% arrange(value.var)
plotLoadings(final.splsda, contrib = 'max', method = 'mean', comp = 1)
plotVar(final.splsda, cutoff = 0.5)
write.table(select.var, file = '/Users/congliu/prs/plsda_var.txt', sep = '\t')

set.seed(42)
perf.final <- perf(final.splsda, validation = "Mfold", folds = 5,
                   dist = 'max.dist', nrepeat = 10,
                   progressBar = FALSE)
perf.final$error.rate
plot(perf.final, col = color.mixo(5))
plot(perf.final$features$stable[[1]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 1', las =2)
plot(perf.final$features$stable[[2]], type = 'h', ylab = 'Stability', 
     xlab = 'Features', main = 'Comp 2', las =2)
# here we match the selected variables to the stable features
ind.match = match(selectVar(final.splsda, comp = 1)$name, 
                  names(perf.final$features$stable[[1]]))
#extract the frequency of selection of those selected variables
Freq = as.numeric(perf.final$features$stable[[1]][ind.match])

res <- data.frame(selectVar(final.splsda, comp = 1)$value, Freq)
write.table(res, file = '/Users/congliu/prs/plsda_res.txt', sep = '\t')


sample_site <- data.frame(final.splsda$variates)[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('plsda1', 'plsda2')
sample_site <- merge(sample_site, mapfile, by.x = 'names', by.y = 'Description', all.x = TRUE)
plsda_result_eig <- {final.splsda$explained_variance$X}[1:2]

final_p <- ggplot(sample_site, aes(plsda1, plsda2, color = Type, label = names)) +
  geom_point(size = 1.5, alpha = 0.6) + 
  stat_ellipse(show.legend = FALSE, linetype = 2) +
  scale_color_manual(values = c('#1D7ACC', '#F67433', '#00815F','#C673FF2E')) +
  scale_shape_manual(values = c(16, 17)) +
  theme(panel.grid = element_line(color = 'grey50'), 
        panel.background = element_rect(color = 'black', fill = 'transparent')) + 
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent')) +
  labs(x = paste('PLS-DA axis1 ( explained variance ', round(100 * plsda_result_eig[1], 2), '% )', sep = ''), 
      y = paste('PLS-DA axis2 ( explained variance ', round(100 * plsda_result_eig[2], 2), '% )', sep = ''))
final_p

# PLSDA模型 -----------------------------------------------------------------
set.seed(42)
DataSet <- subset(otu3, select = -Description)
train_test_index <- caret::createDataPartition(DataSet$Type, p=0.8, list = FALSE)
training_dataset <- DataSet[train_test_index,]
testing_dataset <- DataSet[-train_test_index,]
X_train <- subset(training_dataset, select = -c(Type))
y_train <- training_dataset$Type
X_test <- subset(testing_dataset, select = -c(Type))
y_test <- testing_dataset$Type


## For PLS-DA, train the model，如上文的final.splsda
plsda.train <- plsda(X_train, y_train, ncomp = 4)
# then predict
test.predict <- predict(plsda.train, X_test, dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,4] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = y_test, predicted = prediction)
get.BER(confusion.mat)

# The tuning should be performed on the training set only to avoid overfitting! 
splsda.train <- splsda(X_train, y_train, ncomp = 2, keepX = c(230, 50))
test.predict <- predict(splsda.train, X_test, dist = "all")
# store prediction for the 4th component,为什么那么多空值。。
prediction <- test.predict$class$max.dist[,2] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = y_test, predicted = prediction)
get.BER(confusion.mat)
#auroc(splsda.train, roc.comp = 1)
cbind(Y = as.character(y_test), prediction)










