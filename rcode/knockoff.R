source('genedata.R')
source('penalty.R')

#p <- 5
#k <- 2
#n <- 3000
#rou <- 0.5
#data <- genedata(n,k,p,rou)
#X <- data$X
#U <- data$U
#y <- data$y


record <- c()
lambda_list <- c()
feature_index <- c()
W <- c()


dialectical <- function(X,y, lambda.high, lambda.low, precision, penalmode, family){

  beta.high <- getbeta(X = X, y = y, lambda=lambda.high, penalmode=penalmode, family=family)
  beta.low <- getbeta(X = X, y = y, lambda=lambda.low, penalmode=penalmode, family=family)
  element.new <- setdiff(which(beta.low != 0),which(beta.high != 0))
  
  if(length(element.new)==1){
    local.lambda.low <- lambda.low
    local.lambda.high <- lambda.high
    while(local.lambda.high - local.lambda.low > precision){
      lambda.mid <- (local.lambda.high + local.lambda.low)/2
      beta.mid <- getbeta(X = X, y = y, lambda=lambda.mid, penalmode=penalmode, family=family)
      find.mid <- which(beta.mid != 0)
      find.low <- which(beta.low != 0)
      if(length(setdiff(which(beta.low != 0),which(beta.mid != 0))) == 0){
        local.lambda.low <- lambda.mid
      }
      else{
        local.lambda.high <- lambda.mid
      }
    }
    record <<- c(record,element.new)
    feature_index <<- c(record,element.new)
    lambda_list <<- c(lambda_list,local.lambda.low)
  }
  
  else if(length(element.new)>1){
    lambda.mid <- (lambda.high + lambda.low)/2
    dialectical(X,y,
                lambda.high=lambda.high,
                lambda.low=lambda.mid,
                precision=precision,
                penalmode=penalmode,
                family=family)
    dialectical(X,y,
                lambda.high=lambda.mid,
                lambda.low=lambda.low,
                precision=precision,
                penalmode=penalmode,
                family=family)
    
  }

}


knockoff <- function(X,y,
                     precision=1e-5,
                     penalmode="lasso",
                     family="gaussian"){
  # clear the store memmory
  record <<- c()
  lambda_list <<- c()
  feature_index <<- c()
  W <<- c()
  
  lambda.low <- 0
  lambda.high <- lambda.low
  while (sum(getbeta(X = X, y = y, lambda=lambda.high, penalmode="SCAD", family="gaussian") != 0) != 0){
    lambda.high <- lambda.high+1
  }
  
  beta <- getbeta(X = X, y = y, lambda=0.2, penalmode="SCAD", family="gaussian")
  
  dialectical(X,y,
              lambda.high=lambda.high,
              lambda.low=lambda.low,
              precision=precision,
              penalmode=penalmode,
              family=family)

}

#knockoff(X,y)

#record
#lambda_list



#beta.high <- getbeta(X = X, y = y, lambda=0.2, penalmode="SCAD", family="gaussian")
#beta.low <- getbeta(X = X, y = y, lambda=0.02, penalmode="SCAD", family="gaussian")
