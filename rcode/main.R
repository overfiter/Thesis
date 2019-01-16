source('genedata.R')
source('nonparam.R')
source('knockoff.R')


p <- 5
k <- 2
n <- 3000
rou <- 0.5
data <- genedata(n,k,p,rou)
X <- data$X
XX <- augment(X)
U <- data$U
y <- data$y


g1 <- nonparam(U,y,U,0.3,"spline")
g2 <- 0*XX
for (i in 1:p){
  g2[,i] <- nonparam(U,XX[,i],U,0.3,"spline")
}

knockoff(XX-g2,y-g1)
lambda_list
feature_index

q <- 0.5
num.select <- max(1,floor(q*p))
selected <- feature_index[1:num.select]
X[,selected]

model <- lm((y-g1)~X[,selected])
y.hat <- predict(model,as.data.frame(X[,selected]))
plot(y,y.hat+g1)
