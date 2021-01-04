################### source functions ########################
library(grpreg)
library(plus)
library(BB)
library(ncvreg)
library(mvtnorm)


## 

obj.f <- function(a, theta, beta, delta, ZZ, h) {
  
  n <- nrow(ZZ)
  X <- cbind(1,ZZ[,-1])
  Y <- ZZ[,1]
  d <- ncol(X)
  s <- length(a)
  delta <- matrix(delta, d, s, byrow = F)
  lt <- length(theta)+1
  tmp <- NULL
  for (i in 1:s) {
    tmp[[i]] <- apply(X, 1, function(x, theta, a, delta, h) {
      (x%*%delta[,i]) * pnorm( drop( c(-1, x[2:lt]) %*% c(a[i], theta) ) / h)
    }, theta = theta, a = a, delta = delta, h = h)
  }
  
  sm <- Reduce("+", tmp)
  
  da <- cbind(sm, ZZ)
  Q <- apply(da, 1, function(x) {
    x[2] - c(1,x[-c(1:2)])%*%beta - x[1]
  })
  
  return(sum(Q^2)/n)
  
}


## given beta, theta, estimate a

obj.f_a <- function(a, theta, beta, delta, ZZ, h) {
  
  n <- nrow(ZZ)
  X <- cbind(1,ZZ[,-1])
  Y <- ZZ[,1]
  d <- ncol(X)
  s <- length(a)
  delta <- matrix(delta, d, s, byrow = F)
  lt <- length(theta)+1
  tmp <- NULL
  for (i in 1:s) {
    tmp[[i]] <- apply(X, 1, function(x, a, theta, delta, h) {
      (x%*%delta[,i]) * pnorm( drop( c(-1, x[2:lt]) %*% c(a[i], theta) ) / h)
    }, a = a, theta = theta, delta = delta, h = h)
  }
  
  sm <- Reduce("+", tmp)
  
  da <- cbind(sm, ZZ)
  Q <- apply(da, 1, function(x) {
    x[2] - c(1,x[-c(1:2)])%*%beta - x[1]
  })
  
  return(sum(Q^2)/n)
  
}

## given beta, a, estimate theta

obj.f_theta <- function(theta, a, beta, delta, ZZ, h) {
  
  n <- nrow(ZZ)
  X <- cbind(1,ZZ[,-1])
  Y <- ZZ[,1]
  d <- ncol(X)
  s <- length(a)
  lt <- length(theta)+1
  delta <- matrix(delta, d, s, byrow = F)
  tmp <- NULL
  for (i in 1:s) {
    tmp[[i]] <- apply(X, 1, function(x, theta, a, delta, h) {
      (x%*%delta[,i]) * pnorm( drop( c(-1, x[2:lt]) %*% c(a[i], theta) ) / h)
    }, theta = theta, a = a, delta = delta, h = h)
  }
  
  sm <- Reduce("+", tmp)
  
  da <- cbind(sm, ZZ)
  Q <- apply(da, 1, function(x) {
    x[2] - c(1,x[-c(1:2)])%*%beta - x[1]
  })
  
  return(sum(Q^2)/n)
  
}



## given theta, a, estimate beta
lm.f_beta <- function(theta, a, ZZ, h) {
  
  n <- nrow(ZZ)
  X <- cbind(1,ZZ[,-1,drop = FALSE])
  Y <- ZZ[,1, drop = FALSE]
  d <- ncol(X)
  s <- length(a)
  lt <- length(theta)+1
  
  XX <- X; ind <- NULL
  for(i in 1:s) {
    
    ind[[i]] <- which(drop(X[,2:lt]%*%theta) <= a[i])
    
    XX <- cbind( XX, t(apply(X, 1, function(x, theta, a, h){ x * pnorm(drop( c(-1, x[2:lt]) %*% 
                                                              c(a, theta) ) / h) }, theta= theta, a = a[i], h = h )) )
  }
  
  p.factor <- rep(1, ncol(XX[,-1]))
  for(i in 1:s){
    p.factor[((i+1)*(d-1)-2):((i+1)*(d-1))] <- 0
  }
   
  object <- ncvreg(XX[,-1], c(Y), penalty ="MCP", penalty.factor = p.factor)
  gcv <-  object$loss/(1 - colSums(object$beta!=0)/n)^2 
  lamb=object$lambda[which.min(gcv)]
  gamma.hat <- coef(object, lam = lamb)
  
  beta.hat <- gamma.hat[1:d]
  delta.hat <- gamma.hat[-c(1:d)]
  
  return(list(nind = ind, beta.hat = beta.hat, delta.hat = delta.hat))
  
}


split.fun <- function(X, Y, c) {
  
  n=length(Y)
  p=dim(X)[2]
  m=ceiling(sqrt(n*c))
  q=floor(n/m)
  
  ###########transform###################
  K_temp=matrix(0, nrow = q, ncol=q, byrow=TRUE)
  
  X_temp=cbind(1,X)
  
  Y_temp=c(Y)
  for(i in 1:q)
    K_temp[i,1:i]=rep(1,i)
  
  
  x=NULL
  y=NULL
  x[[1]]=as.matrix(X_temp[1:((n-(q-1)*m)),])
  y[[1]]=Y_temp[1:((n-(q-1)*m))]
  
  for(i in 2:q) {
    x[[i]] = as.matrix(X_temp[(n - (q - i + 1) * m + 1):((n - (q - i) * m)),])
    y[[i]] = Y_temp[(n - (q - i + 1) * m + 1):((n - (q - i) * m))]
  }
  X_temp1 <- lapply(1:length(x), function(j, mat, list)
    kronecker(K_temp[j, , drop = FALSE], x[[j]]), mat = K_temp, list = x)
  Xn=do.call("rbind",X_temp1)
  
  
  Yn=NULL
  for(i in 1:q) {
    Yn=c(Yn,y[[i]])
  }
  
  group <- kronecker(c(1:q),rep(1,p+1))
  colnames(Xn) <- group
  fit <- grpreg::grpreg(Xn[,-1], Yn, group[-1], penalty = "grMCP")
  #  mcp.coef <- grpreg::select(fit,"BIC" )$beta
  #  mcp.lamb <- grpreg::select(fit,"BIC" )$lambda
  
  gcv <-  fit$loss/(1 - colSums(fit$beta!=0)/n)^2 
  mcp.lamb <- fit$lambda[which.min(gcv)]
  mcp.coef <- fit$beta[,which.min(gcv)]
  
  
  #     v = 1.2
  #     BIC.lam <- fit$loss + v*fit$df*log(1/fit$lambda)
  #     mcp.coef <- fit$beta[,which.min(BIC.lam)]
  #     mcp.lamb <- fit$lambda[which.min(BIC.lam)]
  
  
  mcp.coef.s <- sum(abs(mcp.coef))
  
  mcp.coef.v.m <- abs(matrix(c(mcp.coef), q, (p + 1), byrow = TRUE))
  mcp.coef.m <- c(apply(mcp.coef.v.m, 1, max))
  
  mcp.cp <- which(mcp.coef.m!=0)
  
  if (length(mcp.cp) > 1) {
    for (i in 2:length(mcp.cp))
    {
      if (mcp.cp[i] - mcp.cp[i - 1] == 1)
        mcp.cp[i] = 0
    }
  }
  
  mcp.cp1 <- mcp.cp[mcp.cp > 1 & mcp.cp < q]
  
  ## consistent change location interval
  s.est <- length(mcp.cp1)
  if(s.est > 1){
    
    ci <- NULL
    ci[[1]] <- c( (n-(q - mcp.cp1[1]+2)*m) : (n-(q - mcp.cp1[1])*m) )
    numb <- n-(q-1)*m + (mcp.cp1[1] - 2)*m
    for(i in 2:s.est) {
      ci[[i]] = c( (n-(q - mcp.cp1[i]+2)*m) : (n-(q - mcp.cp1[i])*m) )
      numb <- c( numb, n-(q-1)*m+(mcp.cp1[i]-2)*m )
    }
    
    numb <- c(0, numb, n)
    #     rss <- c()
    #     for(i in 1: (length(numb)-1)) {
    #       rss[i] <- lmres( X[(numb[i]+1):numb[i+1],], Y[(numb[i]+1):numb[i+1]])$rss
    #     }
    #     rss <- sum(rss)
    #     BIC <- log(n)*((s.est+1)*(p+1))+n*log(rss/n)
    
    xx_temp <- NULL
    xx_temp[[1]] <- X_temp
    for(i in 1:s.est) {
      xx_temp[[i+1]] <- Xn[,which(colnames(Xn)==mcp.cp1[i])]
    }
    xx <- do.call(cbind, xx_temp)
    rss <- sum(lm(Yn~ xx-1 )$res^2)    
    BIC <- log(n)*((s.est+1)*(p+1))+n*log(rss/n)
    
    
  } else if (s.est == 1) {
    
    ci <- c( (n-(q - mcp.cp1+2)*m) : (n-(q - mcp.cp1)*m) )
    numb <- c(0, n-(q-1)*m + (mcp.cp1 - 2)*m, n)
    #     rss <- c()
    #     for(i in 1: (length(numb)-1)) {
    #       rss[i] <- lmres( X[(numb[i]+1):numb[i+1],], Y[(numb[i]+1):numb[i+1]])$rss
    #     }
    #     rss <- sum(rss)
    #     BIC <- log(n)*((s.est+1)*(p+1)) + n*log(rss/n)
    
    
    xx_temp <- NULL
    xx_temp[[1]] <- X_temp
    for(i in 1:s.est) {
      xx_temp[[i+1]] <- Xn[,which(colnames(Xn)==mcp.cp1[i])]
    }
    xx <- do.call(cbind, xx_temp)
    rss <- sum(lm(Yn~ xx-1 )$res^2)
    BIC <- log(n)*((s.est+1)*(p+1))+n*log(rss/n)
    
    
  }  else {
    
    rss <- sum(lm(Y~ X_temp -1 )$res^2)
    BIC <- log(n)*(p+1)+n*log(rss/n)
    #     BIC <- Inf
    mcp.cp1 = mcp.cp1
    ci = 1:n
    #     rss = Inf
    numb = c(0,n)
  }
  
  
  #     rss <- sum((Yn - Xn%*%mcp.coef)^2)
  #   BIC <- log(n)*((s.est+1)*(p+1))+n*log(rss/n)
  
  return(list(BIC = BIC, mcp.cp1 = mcp.cp1, ci = ci, rss = rss, numb = numb, lamb = mcp.lamb) )
  
  
}



smooth <- function(ZZ, s.est, ini.a = c(-0.3, seq(0.2, 0.8, length.out = s.est-1)), 
                   ini.theta = c(0.75, -0.25, sqrt(1-(-0.75)^2-(-0.25)^2)),
                   tol= 1e-4, K = 20) {
  
  p <- ncol(ZZ)-1
  n <- nrow(ZZ)
  k <- 1; d <- p+1
  #   tol = 1e-3
  ## initial
  a.hat0 <- ini.a
  theta.hat0 <- ini.theta
  #   beta.hat0 <- ini.beta
  #   delta.hat0 <- ini.delta
  h <- n^(-1)*log(n)
  lt <- length(theta.hat0) + 1
  beta.trace <- matrix(0, d, K)
  delta.trace <- matrix(0, s.est*d, K)
  theta.trace <- matrix(0, lt-1, K)
  a.trace <- matrix(0, s.est, K)
  while (k <= K) {
    
    ## given theta, a, estimate beta
    ord=order(ZZ[,2:lt]%*%theta.hat0)
    ZZ=ZZ[ord,]
    fit <- lm.f_beta(theta = theta.hat0, a = a.hat0, ZZ = ZZ, h = h)
    beta.hat <- fit$beta.hat
    delta.hat <- fit$delta.hat
    # nind <- fit$nind
    
    ## given beta, theta, estimate a
    fit.a <- BBoptim(par = a.hat0, fn = obj.f_a, lower = rep(-1,s.est), upper = rep(1,s.est), 
                     theta = theta.hat0, beta = beta.hat, delta = delta.hat, ZZ = ZZ, h = h )
    a.hat <- fit.a$par 
    
    ## given beta, a, estimate theta
    fit.theta <- BBoptim(par = theta.hat0, fn = obj.f_theta, lower = rep(-1,lt-1), upper = rep(1,lt-1),
                         a = a.hat, beta = beta.hat, delta = delta.hat, ZZ = ZZ, h = h )
    theta.hat <- fit.theta$par/norm(fit.theta$par, "2")
    
    beta.trace[,k] <- beta.hat
    delta.trace[,k] <- delta.hat
    theta.trace[,k] <- theta.hat
    a.trace[,k] <- a.hat
    k <- k + 1
    if ( all( abs(a.hat0 - a.hat) < tol, abs(theta.hat0 - theta.hat) < tol)  ) {
      print("ok!")
      break
    } 
    theta.hat0 <- theta.hat
    a.hat0 <- a.hat
    
    # print(k)
  }
  
  ## given theta, a, estimate beta
  ord=order(ZZ[,2:lt]%*%theta.hat)
  ZZo=ZZ[ord,]
  fit <- lm.f_beta(theta = theta.hat, a = a.hat, ZZ = ZZo, h = h)
  beta.hat <- fit$beta.hat
  delta.hat <- fit$delta.hat
  nind <- fit$nind
  
  deviance <- n*obj.f(theta = theta.hat, a = a.hat, beta = beta.hat, delta = delta.hat,
                      ZZ= ZZo, h = h)
  
  
  return(obj <- list(beta.hat = beta.hat, delta.hat = delta.hat, theta.hat = theta.hat, 
                     a.hat = a.hat, nind = nind, deviance = deviance))
  
}


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
    
#     ## choose c = 1
#     split.fit <- split.fun(X, Y, c)
#     mcp.cp1 <- split.fit$mcp.cp1
#     ci <- split.fit$ci
#     s.est <- length(mcp.cp1)
    
    
    if(s.est >= 1) {
      
      ############### step 2
      ## smoothing estimate with initial (theta, a, beta, delta)
      sm.obj <- smooth(ZZ = ZZo, s = s.est, ini.theta = theta.hat0)
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

