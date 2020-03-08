#outlier

library(CORElearn)

train <- read.delim('/Users/congliu/train.txt', sep = '\t') %>% select(-Sample)
md <- CoreModel(Type ~ ., train, model="rf", rfNoTrees=100, 
                maxThreads=1)
outliers <- rfOutliers(md, train)
plot(abs(outliers))
#for a nicer display try 
plot(md, train, rfGraphType="outliers")

destroyModels(md) # clean up



dataset <- iris
md <- CoreModel(Species ~ ., dataset, model="rf", rfNoTrees=30, 
                maxThreads=1)
outliers <- rfOutliers(md, dataset)
plot(abs(outliers))
#for a nicer display try 
plot(md, dataset, rfGraphType="outliers")
