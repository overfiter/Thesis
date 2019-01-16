library(MASS)

genedata <- function(n,k,p,rou){
  sigma <- matrix(0, p, p)
  for (i in 1:p){
    for (j in 1:p){
      sigma[i,j] <- rou^abs(i-j)
    }
  }
  beta <- 5*sample(c(rep(1, k),rep(0, p-k)))
  X <- mvrnorm(n=n, mu=rep(0, p), sigma)
  U <- runif(n, min = 0, max = 2*pi)
  error <- rnorm(n,mean=0,sd=1)
  y <- X %*% beta + 10*sin(U) + error
  y <- unlist(y)
  
  return ( list(X = X,U = U,y = y) )
}


ortho <- function(X){
  n <- dim(X)[1]
  p <- dim(X)[2]
  if(n < 2*p){return (0)}
  
  U <- matrix(1,n,p)
  XU <- cbind(X,U)
  for(i in (1:p -1)){
    A <- XU[(n-p-i+1):n, 1:(p+i)]
    b <- -(t(XU[1:(n-p-i),1:(p+i)]) %*% XU[1:(n-p-i),(p+i+1)])
    XU[(n-p-i+1):n,(p+i+1)] <- solve(a=t(A),b=b)
  }
  return (XU[,((p+1):(2*p))])
}

#X <- genedata(10,3,3,0.5)$X
#cor(X)
#U <- ortho(X)
#t(X) %*% U
#write.csv(X, file = "X.csv", row.names = FALSE)


augment <- function(X){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  U <- ortho(X)
  X.norm <- t( t(X)/sqrt(colSums(X^2)) )
  U.norm <- t( t(U)/sqrt(colSums(U^2)) )
  sigma <- t(X.norm) %*% X.norm
  sigma.inv <- solve(sigma)
  
  s <- min(1, 2*min(eigen(sigma)$values)) * rep(1, p) - 0.0001
  C <- chol(2*diag(s) - diag(s)%*%sigma.inv%*%diag(s))
  X.tilde <- X.norm %*% (diag(rep(1,p)) - sigma.inv %*%diag(s)) + U.norm %*% C
  
  return (cbind(X.norm,X.tilde))
}
#X.tilde <- augment(X)

#cor(X)
#cor(X.tilde)
#cor(cbind(X,X.tilde))

#cbind(X,X.tilde)


