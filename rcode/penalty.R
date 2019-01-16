library(glmnet)
library(ncvreg)

getbeta <- function(X, y, lambda=0, penalmode="lasso", family="gaussian"){
  if(penalmode == "lasso"){
    model <- glmnet(x = X, y = y, family=family,lambda = lambda, alpha = 1)
    beta <- as.data.frame(as.matrix(model$beta))$s0
  }
  else if(penalmode == "SCAD"){
    model <- ncvreg(X = X, y = y, family=family,penalty="SCAD",lambda = lambda, alpha = 1)
    beta <- as.data.frame(model$beta)[-1,1]
  }
  return (beta)
}


#a <- glmnet(x = X-g2, y = y-g1, family = "gaussian",lambda = 0.1, alpha = 1)
#plot(predict(a, X)+g1,y)
#beta <- as.data.frame(as.matrix(a$beta))$s0
#beta



#a <- ncvreg(X = X-g2, y = y-g1, family="gaussian",penalty="SCAD",lambda = 0.1)
#plot(predict(a, X)+g1,y)
#beta <- as.data.frame(a$beta)[-1,1]
#beta
