source("multiplecpfun.R")


####################################
#example, two thresholds.


# generate data
## covariance structure 1,2
dg.fun <- function(n = 300, p = 20, rho = 0.5) {
  
  sig <- rho^abs(row(diag(p))-col(diag(p)))
  X <- rmvnorm(n, sigma = sig)
#   X <- matrix(rnorm(n*p,0,1), n, p)
  
  # theta <- c(-0.25, sqrt(1-(-0.25)^2) )
  theta = c(0.75, -0.25, sqrt(1-(-0.75)^2-(-0.25)^2))
  subX <- X[,1:length(theta)]
  W <- subX%*%theta
  
  #Real threshods (qnorm(0.3),qnorm(0.6))=(-0.5244,0.2533)
  id1=which(W <= qnorm(0.3))
  id2=which(W <= qnorm(0.6) & W > qnorm(0.3))
  id3=which(W > qnorm(0.6))
  
  
  beta01=2
  beta11=c(rep(1,5), rep(0,p-5))
  
  beta02=1
  beta12=c(1,1,0,rep(0,2), rep(0,p-5))
  
  beta03=1
  beta13=c(0,2,0,rep(0,2), rep(0,p-5))
  #Real coefficeints:(beta01,beta11,beta02-beta01,beta12-beta11,beta03-beta02,beta13-beta12)
  #=(2,1,1,1,1,1,-1,0,0,-1,-1,-1,0,-1,1,0,0,0)
  
  X1=X[id1,]
  n1=length(id1)
  Y1=X1%*%beta11+beta01+rnorm(n1,0,sqrt(0.25))
  Z11=cbind(Y1,X1)
  
  X2=X[id2,]
  n2=length(id2)
  Y2=X2%*%beta12+beta02+rnorm(n2,0,sqrt(0.25))
  Z12=cbind(Y2,X2)
  
  X3=X[id3,]
  n3=length(id3)
  Y3=X3%*%beta13+beta03+rnorm(n3,0,sqrt(0.25))
  Z13=cbind(Y3,X3)
  
  ZZ=rbind(Z11,Z12,Z13)
  colnames(ZZ) <- c("y", paste("x", 1:p, sep = "_"))
  
  subg <- cumsum(c(n1, n2))
  
  return( list(ZZ = ZZ, tr.cp = subg) )
}

## covariance structure 3
dg.fun <- function(n = 300, p = 20, rho = 0.5) {
  
  sig <- matrix(rho, p,p) + diag(1-rho,p)
  X <- rmvnorm(n, sigma = sig)
  #   X <- matrix(rnorm(n*p,0,1), n, p)
  
  # theta <- c(-0.25, sqrt(1-(-0.25)^2) )
  theta = c(0.75, -0.25, sqrt(1-(-0.75)^2-(-0.25)^2))
  subX <- X[,1:length(theta)]
  W <- subX%*%theta
  
  #Real threshods (qnorm(0.3),qnorm(0.6))=(-0.5244,0.2533)
  id1=which(W <= qnorm(0.3))
  id2=which(W <= qnorm(0.6) & W > qnorm(0.3))
  id3=which(W > qnorm(0.6))
  
  
  beta01=2
  beta11=c(rep(1,5), rep(0,p-5))
  
  beta02=1
  beta12=c(1,1,0,rep(0,2), rep(0,p-5))
  
  beta03=1
  beta13=c(0,2,0,rep(0,2), rep(0,p-5))
  #Real coefficeints:(beta01,beta11,beta02-beta01,beta12-beta11,beta03-beta02,beta13-beta12)
  #=(2,1,1,1,1,1,-1,0,0,-1,-1,-1,0,-1,1,0,0,0)
  
  X1=X[id1,]
  n1=length(id1)
  Y1=X1%*%beta11+beta01+rnorm(n1,0,sqrt(0.25))
  Z11=cbind(Y1,X1)
  
  X2=X[id2,]
  n2=length(id2)
  Y2=X2%*%beta12+beta02+rnorm(n2,0,sqrt(0.25))
  Z12=cbind(Y2,X2)
  
  X3=X[id3,]
  n3=length(id3)
  Y3=X3%*%beta13+beta03+rnorm(n3,0,sqrt(0.25))
  Z13=cbind(Y3,X3)
  
  ZZ=rbind(Z11,Z12,Z13)
  colnames(ZZ) <- c("y", paste("x", 1:p, sep = "_"))
  
  subg <- cumsum(c(n1, n2))
  
  return( list(ZZ = ZZ, tr.cp = subg) )
}

## simulation function
sim.fun <- function(n = 300, p = 5, rho = 0) {
  
  data <- dg.fun(n = n, p = p, rho = rho)
  ZZ <- data$ZZ
  tr.cp <- data$tr.cp
  
  ############### main iter
  
  obj.fit <- iter.fun(ZZ, ini.theta = c(0.75, -0.25, sqrt(1-(-0.75)^2-(-0.25)^2)), tol = 1e-3, K = 10)
  
  fit <- obj.fit$sm.obj
  s.est <- obj.fit$s.est
  s.trace <- obj.fit$s.trace
  cpl1 <- rev(append(fit$nind, list(c(1:n)) )) 
  id.est <- rev(lapply(1:length(cpl1), function(i) setdiff(cpl1[[i]], unlist(cpl1[-(1:i)]))))
  
  return(list(obj = fit, s.est = s.est, s.trace = s.trace, id.est = id.est, tr.cp = tr.cp))
  
}

##
sim.res <- sim.fun(n = 300, p = 5, rho = 0.5)  

## 500 times run
library(doSNOW)
library(foreach)
cl<-makeCluster(5) #change the 6 to your number of CPU cores
registerDoSNOW(cl)

ptm <- proc.time()
sim.res<- foreach(i=1:500, .packages = c("plus", "grpreg", "BB", "ncvreg", "mvtnorm")) %dopar% {  
  #loop contents here
  sim.res <- sim.fun(n = 300, p = 20, rho = 0.5)  
  return(sim.res) 
  
} 
proc.time()-ptm
stopCluster(cl)




s.est <- sapply(sim.res, function(x) {x$s.est})
length(which(s.est == 2))/500
table(s.est)
theta.hat <- sapply(sim.res[which(s.est == 2)], function(x) {x$obj$theta.hat})
beta.hat <- sapply(sim.res[which(s.est == 2)], function(x) {x$obj$beta.hat})
delta.hat <- sapply(sim.res[which(s.est == 2)], function(x) {x$obj$delta.hat})
a.hat <- sapply(sim.res[which(s.est == 2)], function(x){x$obj$a.hat})


## plot
library(ggplot2)
qplot(c(a.hat), geom="histogram")

qplot(c(theta.hat), geom="histogram")

boxplot.matrix(t(rbind(beta.hat, delta.hat)))


##
n=300
p <- 5
beta <- c(2, rep(1,5), rep(0,p-5))
delta <- c( c(1,1,1,0,rep(0,2), rep(0,p-5))-beta, c(1,0,2,0,rep(0,2), rep(0,p-5)) - c(1,1,1,0,rep(0,2), rep(0,p-5)) )
a = c(qnorm(0.3), qnorm(0.6))
theta = c(0.75, -0.25, sqrt(1-(-0.75)^2-(-0.25)^2))

#### variable selection
coeff <- c(beta, delta)
coeff.hat <- rbind(beta.hat, delta.hat)  
select0 <- apply(coeff.hat, 2, function(x){which(x == 0)})
select1 <- apply(coeff.hat, 2, function(x){which(x != 0)})


## correct 0
zero.id <- which(coeff ==0)
mean(sapply(select0, function(x){ length(intersect(x, zero.id)) }))

## incorrect 0
nonzero.id <- which(coeff !=0)
mean(sapply(select0, function(x){ length(intersect(x, nonzero.id ) ) }))


## bias
bias.theta = rowMeans(apply(theta.hat, 2, function(x) {theta- x}))
bias.a = rowMeans(apply(a.hat, 2, function(x) {a- x}))
bias.beta = rowMeans(apply(beta.hat, 2, function(x) {beta- x}))
bias.delta = rowMeans(apply(delta.hat, 2, function(x) {delta- x}))

## sd
sd.theta = apply(theta.hat, 1, sd)
sd.a = apply(a.hat, 1, sd)
sd.beta = apply(beta.hat, 1, sd)
sd.delta = apply(delta.hat, 1, sd)

## rmse
mse.theta = sqrt(bias.theta^2 + sd.theta^2)
mse.a = sqrt(bias.a^2 + sd.a^2)
rmse.beta = sqrt(bias.beta^2 + sd.beta^2)
rmse.delta = sqrt(bias.delta^2 + sd.delta^2)



