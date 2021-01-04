library(grpreg)
library(plus)
library(BB)
library(ncvreg)
library(mvtnorm)

################### source functions ########################

## step I: given estimation gamma, estimate theta by nls
obj.f <- function(theta, ZZ, beta, delta, h) {
  
  n <- nrow(ZZ)
  X <- cbind(1,ZZ[,-1])
  Y <- ZZ[,1]
  d <- ncol(X)
  lt <- length(theta)
  tmp <- apply(X, 1, function(x, theta, a, delta, h) {
    (x%*%delta) * pnorm( drop(x[1:lt] %*% theta) / h)
  }, theta = theta, delta = delta, h = h)
  
  da <- cbind(tmp, ZZ)
  Q <- apply(da, 1, function(x) {
    x[2] - c(1,x[-c(1:2)])%*%beta - x[1]
  })
  
  return(sum(Q^2)/n)
  
}


## step2: given estimation theta, then estimate (beta, delta) by lm

lm.f_beta <- function(theta, ZZ, h) {
  
  n <- nrow(ZZ)
  X <- cbind(1,ZZ[,-1,drop = FALSE])
  Y <- ZZ[,1, drop = FALSE]
  d <- ncol(X)
  lt <- length(theta)
  
  ind <- which(drop(X[,1:lt]%*%theta) < 0)
  
  XX <- t(apply(X, 1, function(x, theta, h){cbind(x, x*pnorm(x[1:lt]%*%theta/h))}, 
                theta= theta, h = h ))
  
  
  #   object <- plus(XX, c(Y), method="mcp", gamma = 3, intercept = F, normalize = F,eps = 1e-30)
  #   bic=dim(XX)[1]*log(as.vector((1-object$r.square)*sum(Y^2))/length(Y)) + log(dim(XX)[1])*object$dim                                                                                 
  #   lamb=object$lam[which.min(bic)]  
  #   gamma.hat <- coef(object, lam = lamb)
  
  object <- ncvreg(XX[,-1], c(Y), penalty ="MCP")
  bic <-  n*log(object$loss/n)  + colSums(object$beta!=0)*(log(n) + 2 * 1 * log(dim(XX)[2])) ## EBIC
  lamb=object$lambda[which.min(bic)]   
  gamma.hat <- coef(object, lam = lamb)
  
  beta.hat <- gamma.hat[1:d]
  delta.hat <- gamma.hat[-c(1:d)]
  # delta.hat <- matrix(gamma.hat[-c(1:d)], d, s)
  
  return(list(nind = ind, beta.hat = beta.hat, delta.hat = delta.hat))
  
}

## smooth iteration
smooth0 <- function(ZZ, ini.theta = c(-0.15,0.3,sqrt(1-(-0.15)^2-(0.3)^2)), tol= 1e-4, K = 20){
  
  p <- ncol(ZZ)-1
  n <- nrow(ZZ)
  k <- 1; d <- p+1
  #   tol = 1e-3
  ## initial
  theta.hat0 <- ini.theta
  h <- n^(-1)*log(n)
  lt <- length(theta.hat0)
  theta.trace <- matrix(0, lt, K)
  
  while (k <= K) {
    
    ## given theta, a, estimate beta
    ord=order(cbind(1, ZZ[,2:lt])%*%theta.hat0)
    ZZ=ZZ[ord,]
    fit <- lm.f_beta(theta = theta.hat0, ZZ = ZZ, h = h)
    beta.hat <- fit$beta.hat
    delta.hat <- fit$delta.hat
    
    ## given beta, a, estimate theta
    fit.theta <- BBoptim(par = theta.hat0, fn = obj.f, lower = rep(-1,lt), upper = rep(1,lt),
                         beta = beta.hat, delta = delta.hat, ZZ = ZZ, h = h )
    theta.hat <- fit.theta$par/norm(fit.theta$par, "2")
    
    theta.trace[,k] <- theta.hat
    k <- k + 1
    if ( all(abs(theta.hat0 - theta.hat) < tol)  ) {
      print("ok!")
      break
    } 
    theta.hat0 <- theta.hat
    # print(k)
    
  }
  
  ## given theta, a, estimate beta
  ord=order(cbind(1, ZZ[,2:lt])%*%theta.hat)
  ZZo=ZZ[ord,]
  fit <- lm.f_beta(theta = theta.hat, ZZ = ZZo, h = h)
  beta.hat <- fit$beta.hat
  delta.hat <- fit$delta.hat
  nind <- fit$nind
  
  deviance <- n*obj.f(theta = theta.hat, beta = beta.hat, delta = delta.hat,
                      ZZ= ZZo, h = h)
  
  return(list(beta.hat = beta.hat, delta.hat = delta.hat, theta.hat = theta.hat, 
              nind = nind, deviance = deviance))
  
}

## simulation data generate
dgp <- function(n = 300, p = 20, rho = 0 ) {
  
  # generate data
  
  sig <- rho^abs(row(diag(p))-col(diag(p)))
  X <- rmvnorm(n, sigma = sig)
  theta = c(-0.15,0.3,sqrt(1-(-0.15)^2-(0.3)^2))
  
  
  #Real threshods qnorm(0.5) = 0
  subX <- cbind(1,X[,1:(length(theta)-1)])
  W <- subX%*%theta
  id1 <- which(W < qnorm(0.5))
  id2 <- which(W >= qnorm(0.5))
  
  beta01=2
  beta11=c(rep(1,5), rep(0,p-5))
  
  beta02=1
  beta12=c(1,1,0,rep(0,2), rep(0,p-5))
  
  X1=X[id1,]
  n1=length(id1)
  Y1=X1%*%beta11+beta01+rnorm(n1,0,sqrt(0.25))
  Z11=cbind(Y1,X1)
  
  X2=X[id2,]
  n2=length(id2)
  Y2=X2%*%beta12+beta02+rnorm(n2,0,sqrt(0.25))
  Z12=cbind(Y2,X2)
  
  ZZ <- rbind(Z11,Z12)
  colnames(ZZ) <- c("y", paste("x", 1:p, sep = "_"))
  
  
  return(list(ZZ = ZZ, cl = n1))
}

dgp <- function(n = 300, p = 20, rho = 0 ) {
  
  # generate data
  
  sig <- matrix(rho, p,p) + diag(1-rho,p)
  X <- rmvnorm(n, sigma = sig)
  theta = c(-0.15,0.3,sqrt(1-(-0.15)^2-(0.3)^2))
  
  
  #Real threshods qnorm(0.5) = 0
  subX <- cbind(1,X[,1:(length(theta)-1)])
  W <- subX%*%theta
  id1 <- which(W < qnorm(0.5))
  id2 <- which(W >= qnorm(0.5))
  
  beta01=2
  beta11=c(rep(1,5), rep(0,p-5))
  
  beta02=1
  beta12=c(1,1,0,rep(0,2), rep(0,p-5))
  
  X1=X[id1,]
  n1=length(id1)
  Y1=X1%*%beta11+beta01+rnorm(n1,0,sqrt(0.25))
  Z11=cbind(Y1,X1)
  
  X2=X[id2,]
  n2=length(id2)
  Y2=X2%*%beta12+beta02+rnorm(n2,0,sqrt(0.25))
  Z12=cbind(Y2,X2)
  
  ZZ <- rbind(Z11,Z12)
  colnames(ZZ) <- c("y", paste("x", 1:p, sep = "_"))
  
  
  return(list(ZZ = ZZ, cl = n1))
}





## simulation function
sim.fun <- function(n = 300, p = 20, rho = 0) {
  
  dat <- dgp(n = n, p = p, rho = rho)
  ZZ <- dat$ZZ
  Y <- ZZ[,1]
  X <- ZZ[,-1]
  p <- ncol(ZZ)-1
  tr.cl <- dat$cl
  
  ## initial estimation by equally split data
#   ini.theta <- c(-0.15,0.3,sqrt(1-(-0.15)^2-(0.3)^2))
  ini.theta <- c(0.5,0.5,sqrt(1-(0.5)^2-(0.5)^2))

  
  ## main
  
  fit <- smooth0(ZZ = ZZ, ini.theta = ini.theta)
  
  return(list(fit = fit, tr.cl= tr.cl ))
  
}

ptm <- proc.time()
res <- replicate(500, sim.fun(n = 300, p = 5, rho = 0) )
proc.time() - ptm

###
theta.hat <- apply(res, 2, function(x) {x$fit$theta.hat})
beta.hat <- apply(res, 2, function(x) {x$fit$beta.hat})
delta.hat <- apply(res, 2, function(x) {x$fit$delta.hat})

coeff <- rbind(beta.hat, delta.hat)  
select0 <- apply(coeff, 2, function(x){which(x == 0)})
select1 <- apply(coeff, 2, function(x){which(x != 0)})


## correct 0
mean(sapply(select0, function(x){ length(intersect(x, c(8,9))) }))

## incorrect 0
mean(sapply(select0, function(x){ length(intersect(x, setdiff(1:12, 8:9) ) ) }))


####################bias, sd and mse for theta.hat, beta.hat, delta.hat
beta = c(2,1,1,1,1,1)
delta = c(1,1,1,0,0,0)-beta
theta = c(-0.15,0.3,sqrt(1-(-0.15)^2-(0.3)^2))


## bias
bias.theta = rowMeans(apply(theta.hat, 2, function(x) {theta- x}))
bias.beta = rowMeans(apply(beta.hat, 2, function(x) {beta- x}))
bias.delta = rowMeans(apply(delta.hat, 2, function(x) {delta- x}))


## sd
sd.theta =  apply(theta.hat, 1, sd)
sd.beta = apply(beta.hat, 1, sd)
sd.delta = apply(delta.hat, 1, sd)


## rmse
rmse.theta = sqrt(bias.theta^2+ sd.theta^2)
rmse.beta = sqrt(bias.beta^2 + sd.beta^2)
rmse.delta = sqrt(bias.delta^2 + sd.delta^2)


paste0(round(rbind(bias.beta, sd.beta,rmse.beta),3),"&", collapse = "" )

paste0(round(rbind(bias.theta,sd.theta,rmse.theta), 3), "&", collapse = " " )


# load("single_sim_n300p5d3e3")
################# misclassification rate
n=300
id.fun <- function(x) {
  id <- rep(1,n)
  id[x] <- 0
  id
}

##
cl.hat <- apply(res, 2, function(x){ x$fit$nind } )
cl.tr <- apply(res, 2, function(x){ c(1:x$tr.cl)})

group.hat <- sapply(cl.hat, id.fun, simplify = F)
group.tr <- sapply(cl.tr, id.fun, simplify = F)
mean(mapply(function(x,y){mean(abs(x-y))}, group.hat, group.tr ))
sd(mapply(function(x,y){mean(abs(x-y))}, group.hat, group.tr ))




### NMI
NMI.fun <- function(clusters, true_labels) {
  
  tbl = table(clusters, true_labels)
  conv_df = as.data.frame.matrix(tbl)
  mutual_information = 0.0
  for (i in 1:nrow(conv_df)) {
    
    for (j in 1:ncol(conv_df)) {
      
      if (conv_df[i,j] > 0.0) {
        
        mutual_information = mutual_information + ((conv_df[i,j] / sum(tbl)) * log((sum(tbl) * 
                                  conv_df[i,j]) / (sum(conv_df[i,]) * sum(conv_df[,j]))))
      }
    }
  }
  entr_cluster = sum(apply(conv_df, 1, function(x) -(sum(x) / sum(tbl)) * log(sum(x) / sum(tbl))))
  
  entr_class = sum(apply(conv_df, 2, function(x) -(sum(x) / sum(tbl)) * log(sum(x) / sum(tbl))))
  
  NMI = (mutual_information / ((entr_cluster + entr_class) / 2.0))
  
  return(NMI)
  
}

boxplot(mapply(NMI.fun, group.hat, group.tr))

dat = mapply(NMI.fun, group.hat, group.tr)
dat <- stack(as.data.frame(dat))

library(ggplot2)
xlab11 = expression(paste("Example 1", " (", Sigma[1], ", ", "n=150", ", ", "p=5", ")" ))
ex1n150p5 <- ggplot(dat) + ggtitle(xlab11)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

xlab12 = expression(paste("Example 1", " (", Sigma[1], ", ", "n=300", ", ", "p=5", ")" ))
ex1n300p5 <- ggplot(dat) + ggtitle(xlab12)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )


xlab11s2 = expression(paste("Example 1", " (", Sigma[2], ", ", "n=150", ", ", "p=5", ")" ))
ex1n150p5s2 <- ggplot(dat) + ggtitle(xlab11s2)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

xlab12s2 = expression(paste("Example 1", " (", Sigma[2], ", ", "n=300", ", ", "p=5", ")" ))
ex1n300p5s2 <- ggplot(dat) + ggtitle(xlab12s2)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

xlab11s3 = expression(paste("Example 1", " (", Sigma[3], ", ", "n=150", ", ", "p=5", ")" ))
ex1n150p5s3 <- ggplot(dat) + ggtitle(xlab11s3)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

xlab12s3 = expression(paste("Example 1", " (", Sigma[3], ", ", "n=300", ", ", "p=5", ")" ))
ex1n300p5s3 <- ggplot(dat) + ggtitle(xlab12s3)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )


