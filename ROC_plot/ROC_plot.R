library(pROC)
library(ggplot2)
setwd("/Users/congliu/prs_R/ROC_plot")
data <- read.table(file = "ROC_plot_input.txt", header = T, sep = "\t", row.names = 1)

#colnames(data)[2]
rocobj1 <- roc(data[,1], data[,2])
auc(rocobj1)
# plot(rocobj1)
# plot(rocobj1, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2), max.auc.polygon=TRUE,
#      auc.polygon.col="skyblue", print.thres=TRUE)


plot(rocobj1, print.auc=TRUE, auc.polygon=TRUE, partial.auc=c(1, 0.8),
     partial.auc.focus="sp", grid=c(0.1, 0.2), grid.col=c("green", "red"),
     max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)


mycoords <- coords(rocobj1, "all")
plot(mycoords["threshold",], mycoords["specificity",], type="l", 
     col="red", xlab="Cutoff", ylab="Performance")
lines(mycoords["threshold",], mycoords["sensitivity",], type="l", 
      col="blue")
legend(100, 0.4, c("Specificity", "Sensitivity"), 
       col=c("red", "blue"), lty=1)
best.coords <- coords(rocobj1, "best", best.method="youden", transpose = FALSE)
abline(v=best.coords["threshold"], lty=2, col="grey")
abline(h=best.coords["specificity"], lty=2, col="red")
abline(h=best.coords["sensitivity"], lty=2, col="blue")
