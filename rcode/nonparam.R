## 生成数据的函数
#n = 300
#U <- runif(n,min=0,max=2*pi)
#error <- rnorm(n,mean=0,sd=1)
#y <- 3*sin(U) + error
#plot(U, y)

nonparam <- function(U,y,U_hat,h,way = "kernal"){
   if (way == "kernal"){
     model <- ksmooth(U, y, "normal", bandwidth = 0.3, x.points = U_hat)
     macher <- data.frame(raw = U_hat,ind = 1:length(U_hat))
     ind <- macher[order(macher$raw),]$ind
     temp <- data.frame(x = model$x,y = model$y, ind = ind)
     result <- temp[order(temp$ind),]$y
     return (result)
   }
  else if(way == "loess"){
    model <- loess(y ~ U, span = 0.3, family="gaussian")
    result <- predict(model,U_hat)
    return (result)
  }
  else if(way == "spline"){
    model <- smooth.spline(U, y)
    result <- predict(model, U_hat)$y
    return (result)
  }
}

#plot(U,nonparam(U,y,U,0.3,"loess"))


#### NW-kernal
#library(KernSmooth)
#model <- ksmooth(U, y, "normal", bandwidth = 0.3, x.points = U)
#macher <- data.frame(raw = U,ind = 1:length(U))
#ind <- macher[order(macher$raw),]$ind
#test <- data.frame(x = model$x,y = model$y, ind = ind)
#result <- test[order(test$ind),]$y
#plot(U, result)

#### LOESS
#a <- loess(y ~ U, span = 0.3, family="gaussian")
#predict(a, U, se = FALSE)
#plot(U, predict(a, U))

#### 
#library(splines)
#a <- smooth.spline(U, y)
#predict(a, U)
#plot(U, predict(a, U)$y)

