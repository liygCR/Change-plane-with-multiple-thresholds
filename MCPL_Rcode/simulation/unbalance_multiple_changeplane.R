#example 4,two thresholds.
source("multiplecpfun.R")


iter.fun <- function(ZZ, ini.theta, tol= 1e-4, K= 10){
  
  theta.hat0 <- ini.theta
  lt <- length(theta.hat0) + 1
  theta.trace <- matrix(0, lt-1, K) ## fixed dimnesional theta
  s.trace <- matrix(0, 1, K)
  # ci <- list()
  k <- 1
  while (k <= K){
    
    ord=order(ZZ[,2:lt]%*%theta.hat0)
    ZZo=ZZ[ord,]
    
    n=dim(ZZo)[1]
    p=dim(ZZo)[2]-1
    Y=ZZo[,1]
    X=ZZo[,-1]
    #     Word = ZZo[,2:4]%*%c(1,theta.hat0)
    
    
    ## choose the best split c
    c <- seq(0.5,1.5,0.1)
    split.fit <- NULL
    for(j in 1:length(c)) {
      split.fit[[j]] <- split.fun(X, Y, c[j])
    }
    BIC <- sapply(split.fit, function(x){x$BIC}, simplify = T)
    mcp.cp1 <- sapply(split.fit, function(x){x$mcp.cp1}, simplify = F)
    mcp.cp1 <- mcp.cp1[[which.min(BIC)]]
    ci <- split.fit[[which.min(BIC)]]$ci
    s.est <- length(mcp.cp1)
    
    #     ## choose c = 1 in example 1
    #     split.fit <- split.fun(X, Y, c=1)
    #     mcp.cp1 <- split.fit$mcp.cp1
    #     ci <- split.fit$ci
    #     s.est <- length(mcp.cp1)
    
    
    if(s.est >= 1) {
      
      ############### step 2
      ## smoothing estimate with initial (theta, a, beta, delta)
      sm.obj <- smooth(ZZ = ZZo, s = s.est, ini.theta = theta.hat0, 
                       ini.a = c(-0.7, seq(0.7, 1, length.out = s.est-1)))
      theta.hat <- sm.obj$theta.hat
      
      theta.trace[,k] <- theta.hat
      s.trace[k] <- s.est
      k <- k + 1
      if ( all(theta.hat0 - theta.hat < tol) ) {
        print("Finish!")
        break
      } 
      theta.hat0 <- theta.hat
      
      
    } else {
      
      sm.obj <- lm(Y~X)
      break
      
    }  
    
    
  }
  
  
  return(list(s.est = s.est, sm.obj = sm.obj, s.trace = s.trace, theta.trace = theta.trace))
  
}


## generate data
## covariance structure 1,2
dg.fun <- function(n = 300, p = 20, rho = 0.5) {
  
  sig <- rho^abs(row(diag(p))-col(diag(p)))
  X <- rmvnorm(n, sigma = sig)
#   X <- matrix(rnorm(n*p,0,1), n, p)
  
  # theta <- c(-0.25, sqrt(1-(-0.25)^2) )
  theta = c(0.75, -0.25, sqrt(1-(-0.75)^2-(-0.25)^2))
  subX <- X[,1:length(theta)]
  W <- subX%*%theta
  
  #Real threshods (-1,1)
  id1=which(W <= -sqrt(2)/2)
  id2=which(W <= sqrt(2)/2 & W > -sqrt(2)/2)
  id3=which(W > sqrt(2)/2)
  
  
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
  
  #Real threshods (-1,1)
  id1=which(W <= -1)
  id2=which(W <= 0.5 & W > -1)
  id3=which(W > 0.5)
  
  
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

##
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
sim <- replicate(n=500, sim.fun(n = 300, p = 5, rho=0))

sim.res <- lapply(seq(dim(sim)[2]), function(x)sim[,x])

##	
n=300
p <- 5
beta <- c(2, rep(1,5), rep(0,p-5))
delta <- c( c(1,1,1,0,rep(0,2), rep(0,p-5))-beta, c(1,0,2,0,rep(0,2), rep(0,p-5)) - c(1,1,1,0,rep(0,2), rep(0,p-5)) )
a = c(-sqrt(2)/2, sqrt(2)/2)
theta = c(0.75, -0.25, sqrt(1-(-0.75)^2-(-0.25)^2))



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


