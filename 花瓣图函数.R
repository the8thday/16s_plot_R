#构建作图函数（参考自 https://www.cnblogs.com/xudongliang/p/7884667.html）
#花瓣图函数
flower_plot <- function(sample, otu_num, core_otu, start, a, b, r, ellipse_col, circle_col) {
  par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(1,1,1,1))
  plot(c(0,10),c(0,10),type='n')
  n   <- length(sample)
  deg <- 360 / n
  res <- lapply(1:n, function(t){
    draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180), 
                 y = 5 + sin((start + deg * (t - 1)) * pi / 180), 
                 col = ellipse_col[t],
                 border = ellipse_col[t],
                 a = a, b = b, angle = deg * (t - 1))
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         otu_num[t])
    
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) - start,
           adj = 1,
           cex = 1
      )
    } else {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) + start,
           adj = 0,
           cex = 1
      )
    }            
  })
  draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)
  text(x = 5, y = 5, paste('Core:', core_otu))
}


library(forestFloor)
library(randomForest)
library(rgl)
set.seed(1)

X<-data.frame(replicate(2,runif(2000)-.5))
y<--sqrt((X[,1])^4+(X[,2])^4)
Col<-fcol(X,1:2) #make colour pallete by x1 and x2
#insert outlier2 and colour it black
y2<-y;Col2<-Col
y2[1:100]<-rnorm(100,200,1);    #outliers
Col[1:100]="#000000FF" #black

#plot training set
plot3d(X[,1],X[,2],y,col=Col)
rf=randomForest(X,y)    #RF on clean data
rg=randomForest(X,y2)   #RF on contaminated data
vec.plot(rg,X,1:2,col=Col,grid.lines=200)
mean(abs(rf$predict[-c(1:100)]-y[-c(1:100)]))
mean(abs(rg$predict[-c(1:100)]-y2[-c(1:100)]))




