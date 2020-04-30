#！/urs/bin/env Rscript

#http://www.sthda.com/english/articles/36-classification-methods-essentials/150-stepwise-logistic-regression-essentials-in-r/
library(tidyverse)
library(MASS)
library(bestglm)

data("biopsy")

table(biopsy$class)
# 结果485 benign; 241 malignant，分类分布还可以

biopsy.omit <- na.omit(biopsy) %>% select(-ID)
# 在线性模型中，考虑分析共线性问题，即是变量之间相关性
corrplot::corrplot.mixed(cor(biopsy.omit[,1:9]))

# module -----------------------------------------------------------------

set.seed(42)

train_test_index <- caret::createDataPartition(biopsy.omit$class, p = 0.8, list = FALSE)
train_data <- biopsy.omit[train_test_index,]
test_data <- biopsy.omit[-train_test_index,]

lr_model <- glm(class~., family = binomial, data = train_data)
# >% stepAIC(trace = FALSE)
# 查看95%置信区间
confint(lr_model)
exp(coef(lr_model))
#car::vif(lr_model) 用以区分共线性问题，大于5为共线性
probabilities <- lr_model %>% predict(test_data, type = 'response')
predicted.classes <- ifelse(probabilities > 0.5, 'malignant', 'benign')
# Model accuracy
mean(predicted.classes==test_data$class)

# stepwise logistic models
step.model <- lr_model %>% stepAIC(trace = FALSE)
coef(step.model)


# 特征选择 --------------------------------------------------------------------

train_data <- train_data %>% mutate_at(vars(class), ~recode_factor(., !!!c(benign = 0, malignant=1)))
#train_data %>% mutate_at(.var = vars(class), .funs = forcats::fct_recode, 'benign' = '0', 'malignant'='1')
best_glm <- bestglm(train_data,
        IC = 'BIC',
        CVArgs = list(Method = "HTF",K=10,REP=1), #HTF表示方法为K折交叉验证，K指定了均分的份数量
        family = binomial)
# summary(best_glm)
best_glm
# 结果中V1 V2 V6 V7为较好的特征
BIC.fit <- glm(formula = class ~ V1 + V2 +V4+ V7, family = binomial, data = train_data)
test.BIC.pro <- predict(BIC.fit, newdata = test_data, type = "response")
# probabilities <- BIC.fit %>% predict(test_data, type = 'response')
# predicted.classes <- ifelse(probabilities > 0.05, 'malignant', 'benign')
# mean(predicted.classes==test_data$class)

test.BIC.pro <- ifelse(test.BIC.pro>= 0.05,1,0)#预测概率大于0.05默认为肿瘤发生，同样设置为1
test_data <- test_data %>% mutate_at(vars(class), ~dplyr::recode(., !!!c(benign = 0, malignant=1)))
pred.BIC <- ROCR::prediction(test.BIC.pro, test_data$class)
pred.BIC <- ROCR::performance(pred.BIC,"tpr","fpr")
ROCR::plot(pred.BIC,main = "ROC", add = TRUE, color=2)
# 混淆矩阵
caret::confusionMatrix()







