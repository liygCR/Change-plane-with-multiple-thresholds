source("ACTG_fun.R")

####
library(speff2trial)

data <- ACTG175

dat <- data.frame(cd420 = data$cd420, age = data$age, wtkg = data$wtkg, 
                  karnof = data$karnof, cd40 =data$cd40, cd80 = data$cd80, 
                  hemo = as.factor(data$hemo), homo = as.factor(data$homo), 
                  drugs = as.factor(data$drugs), race = as.factor(data$race),
                  gender = as.factor(data$gender), str2 = as.factor(data$str2), 
                  symptom = as.factor(data$symptom), treatment = as.factor(data$arms))
rownames(dat) <- data$pidnum

tmp <- apply(dat[,c(1:6)], 2, scale)
dd <- cbind(tmp, dat[, -c(1:6)])

X <- model.matrix(~ hemo + gender + cd40+ str2 + age + wtkg+ karnof+ cd80 +  homo +   
                    drugs + race +  symptom + treatment , dd)
ZZ <- data.frame( cd420 = dd$cd420, X[,-1] )


n <- nrow(ZZ)
theta0 = c(-0.25, -0.2, 0.85, -0.25, 0.2 ) ## initial estimate theta
obj.fit <- iter.fun( as.matrix(ZZ), ini.theta = theta0, tol = 1e-3, K = 10)


fit <- obj.fit$sm.obj
s.est <- obj.fit$s.est
s.trace <- obj.fit$s.trace
cpl <- rev(append(fit$nind, list(c(1:nrow(ZZ)))))
names(cpl[[1]]) <- rownames(ZZ) 
id.est <- rev(lapply(1:length(cpl), function(i) setdiff(cpl[[i]], unlist(cpl[-(1:i)]))))
sapply(id.est, length)

beta.hat <- obj.fit$sm.obj$beta.hat
delta.hat <- obj.fit$sm.obj$delta.hat
a.hat <- obj.fit$sm.obj$a.hat
theta.hat <- obj.fit$sm.obj$theta.hat



######################### in-sample mse ##############
xx <- cbind(1, ZZ[,-1])
XX <- xx
for(i in 1:s.est){
  XX <- cbind(XX, t(apply(xx, 1, function(x, theta, a, h){
    x * (drop( x[2:(length(theta)+1)]%*%theta > a[i]) )
  }, theta = theta.hat, a = a.hat )))
}

gamma.hat<- c(beta.hat, delta.hat)
y <- ZZ[,1]
y.fit <- as.matrix(XX) %*% gamma.hat 
mse <- mean((y - y.fit)^2)

##### sd and p-value for our
n <- nrow(ZZ)
coef.id <- which(gamma.hat!=0)
sigma <- sqrt(mse)

XXX <- as.matrix(XX)
sd <- sigma * sqrt(diag(solve( t(XXX[,coef.id])%*%XXX[,coef.id] )))

t.value <- gamma.hat[coef.id]/sd
p.value <- 2*pt(-abs(t.value), n-length(coef.id))

round(cbind(coef=gamma.hat[coef.id], sd, t.value, p.value),3)


########### lm mse ###########################
lm.fit <- lm(cd420 ~ ., ZZ)
summary(lm.fit)

## in-sample
training <- ZZ
y.train <- training[,1]
y.fit.lm <- predict(lm.fit, training)
mse.is.lm <- mean((y.train - y.fit.lm)^2)
