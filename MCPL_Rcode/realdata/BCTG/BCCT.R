source("BCCT_fun.R")


#### 
data <- read.csv("mrss.csv")
data <- na.omit(data)
dat <- cbind(mrss=data$mrss, data[, -c(1:2,14)])

X = model.matrix(~ haq+pga+dlcop+fvcp+over+pain+fev1p+ethnic+sex+durdis+age, data = dat)


ZZ <- cbind(scale(dat[,1]), apply(X[, -c(1,9,10)], 2, scale), X[,9:10] )
colnames(ZZ)[1] = "mrss"

ZZ <- as.data.frame(ZZ)



## our
theta0 = c(0.75, -0.25, sqrt(1-(0.8)^2-(-0.25)^2))
obj.fit <- iter.fun( as.matrix(ZZ), ini.theta = theta0, tol = 1e-3, K = 10)

fit <- obj.fit$sm.obj
s.est <- obj.fit$s.est
s.trace <- obj.fit$s.trace
cpl <- rev(append(fit$nind, list(c(1:nrow(ZZ)))))
names(cpl[[1]]) <- rownames(ZZ) 
id.est <- rev(lapply(1:length(cpl), function(i) setdiff(cpl[[i]], unlist(cpl[-(1:i)]))))
sapply(id.est, length) ## locations

beta.hat <- obj.fit$sm.obj$beta.hat
delta.hat <- obj.fit$sm.obj$delta.hat
a.hat <- obj.fit$sm.obj$a.hat
theta.hat <- obj.fit$sm.obj$theta.hat
subX <- as.matrix(ZZ)[,2:(length(theta.hat)+1)]



########### MSE 
training <- ZZ
##
xx <- cbind(1, training[,-1])
X.train <- xx
for(i in 1:s.est){
  X.train <- cbind(X.train, t(apply(xx, 1, function(x, theta, a, h){
    x * (drop( x[2:(length(theta)+1)]%*%theta > a[i]) )
  }, theta = theta.hat, a = a.hat )))
}

y.train <- training[,1]
y.fit <- as.matrix(X.train) %*% c(beta.hat, delta.hat) 
mse.is <- mean((y.train - y.fit)^2)


##### sd and p-value for our
coef <- c(beta.hat, delta.hat)
coef.id <- which(coef!=0)


XX <- as.matrix(X.train)
sigma <- sqrt(nrow(ZZ)*mse.is/(nrow(ZZ)- ncol(XX)))
se <- sigma * sqrt(diag(solve( t(XX[,coef.id])%*%XX[,coef.id] )))

n <- nrow(XX)
t.value <- coef[coef.id]/se
p.value <- 2*pt(-abs(t.value), n-length(coef.id))

round(cbind(coef=coef[coef.id], se, t.value, p.value),3)



## lm
lm.fit <- lm(mrss ~ ., ZZ)
summary(lm.fit)

y.fit.lm <- predict(lm.fit, training)
mse.is.lm <- mean((y.train - y.fit.lm)^2)



