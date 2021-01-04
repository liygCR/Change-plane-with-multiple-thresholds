source("multiplecpfun.R")


####################################
#example 3,no thresholds.


dg.fun <- function(n = 300, p = 20, rho = 0.5) {
  
  sig <- rho^abs(row(diag(p))-col(diag(p)))
  X <- rmvnorm(n, sigma = sig)
  theta = c(0.75, -0.25, sqrt(1-(-0.75)^2-(-0.25)^2))
  subX <- X[,1:length(theta)]
  W <- subX%*%theta
  
  beta01=2
  beta11=c(rep(1,5), rep(0,p-5))
    
  Y=X%*%beta11+beta01+rnorm(n,0,sqrt(0.25))
  ZZ=cbind(Y,X)
  
  colnames(ZZ) <- c("y", paste("x", 1:p, sep = "_"))
  
  return( list(ZZ = ZZ) )
}


dg.fun <- function(n = 300, p = 20, rho = 0.5) {
  
  sig <- matrix(rho, p,p) + diag(1-rho,p)
  X <- rmvnorm(n, sigma = sig)

  theta = c(0.75, -0.25, sqrt(1-(-0.75)^2-(-0.25)^2))
  subX <- X[,1:length(theta)]
  W <- subX%*%theta
  
  beta01=2
  beta11=c(rep(1,5), rep(0,p-5))
    
  Y=X%*%beta11+beta01+rnorm(n,0,sqrt(0.25))
  ZZ=cbind(Y,X)
  
  colnames(ZZ) <- c("y", paste("x", 1:p, sep = "_"))
  
  return( list(ZZ = ZZ) )
}



sim.fun <- function(n = 300, p = 5, rho = 0) {
  
  data <- dg.fun(n = n, p = p, rho = rho)
  ZZ <- data$ZZ
  
  ############### main iter
  
  obj.fit <- iter.fun(ZZ, ini.theta = c(0.5, 0.5, sqrt(1-(0.5)^2-(-0.5)^2)), tol = 1e-3, K = 10)
  
  fit <- obj.fit$sm.obj
  s.est <- obj.fit$s.est

  return(list(obj = fit, s.est = s.est))
  
}


sim <- replicate(n=500, sim.fun(n = 300, p = 5, rho=0.5))

sim.res <- lapply(seq(dim(sim)[2]), function(x)sim[,x])

s.est <- sapply(sim.res, function(x) {x$s.est})
length(which(s.est == 2))/500
table(s.est)



