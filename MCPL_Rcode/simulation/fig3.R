# load("unbalance_sim_n500p5d3_e3")

n=500
p <- 5
beta <- c(2, rep(1,5), rep(0,p-5))
delta <- c( c(1,1,1,0,rep(0,2), rep(0,p-5))-beta, c(1,0,2,0,rep(0,2), rep(0,p-5)) - c(1,1,1,0,rep(0,2), rep(0,p-5)) )
a = c(qnorm(0.3), qnorm(0.6))
theta = c(0.75, -0.25, sqrt(1-(-0.75)^2-(-0.25)^2))
gamma0 <- c(beta, delta)


s.est <- sapply(sim.res, function(x) {x$s.est})
length(which(s.est == 2))/500
table(s.est)
theta.hat <- sapply(sim.res[which(s.est == 2)], function(x) {x$obj$theta.hat})
beta.hat <- sapply(sim.res[which(s.est == 2)], function(x) {x$obj$beta.hat})
delta.hat <- sapply(sim.res[which(s.est == 2)], function(x) {x$obj$delta.hat})
a.hat <- sapply(sim.res[which(s.est == 2)], function(x){x$obj$a.hat})


dd <- as.data.frame(t(rbind(beta.hat, delta.hat)))
colnames(dd) <- 1:18
dat <- stack(dd)
colnames(dat) <- c("Coefficients", "Index")

### plot ##
xlab21 = expression(paste("Example 2", " (", Sigma[1], ", ", "n=150", ")" ))
ex2n150p5 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab21)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5))


xlab22 = expression(paste("Example 2", " (", Sigma[1], ", ", "n=300", ")" ))
ex2n300p5 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab22)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5))

xlab23 = expression(paste("Example 2", " (", Sigma[1], ", ", "n=500", ")" ))
ex2n500p5 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab23)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) ) 

xlab21s2 = expression(paste("Example 2", " (", Sigma[2], ", ", "n=150", ")" ))
ex2n150p5s2 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab21s2)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) ) 

xlab22s2 = expression(paste("Example 2", " (", Sigma[2], ", ", "n=300", ")" ))
ex2n300p5s2 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab22s2)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5)) 


xlab23s2 = expression(paste("Example 2", " (", Sigma[2], ", ", "n=500", ")" ))
ex2n500p5s2 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab23s2)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) )


xlab21s3 = expression(paste("Example 2", " (", Sigma[3], ", ", "n=150", ")" ))
ex2n150p5s3 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab21s3)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) )


xlab22s3 = expression(paste("Example 2", " (", Sigma[3], ", ", "n=300", ")" ))
ex2n300p5s3 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab22s3)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) )

xlab23s3 = expression(paste("Example 2", " (", Sigma[3], ", ", "n=500", ")" ))
ex2n500p5s3 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab23s3)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) )


## 
xlab41 = expression(paste("Example 4", " (", Sigma[1], ", ", "n=150", ")" ))
ex4n150p5 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab41)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) )

xlab42 = expression(paste("Example 4", " (", Sigma[1], ", ", "n=300", ")" ))
ex4n300p5 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab42)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) )

xlab43 = expression(paste("Example 4", " (", Sigma[1], ", ", "n=500", ")" ))
ex4n500p5 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab43)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) )


xlab41s2 = expression(paste("Example 4", " (", Sigma[2], ", ", "n=150", ")" ))
ex4n150p5s2 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab41s2)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) )

xlab42s2 = expression(paste("Example 4", " (", Sigma[2], ", ", "n=300", ")" ))
ex4n300p5s2 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab42s2)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) )

xlab43s2 = expression(paste("Example 4", " (", Sigma[2], ", ", "n=500", ")" ))
ex4n500p5s2 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab43s2)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) )


xlab41s3 = expression(paste("Example 4", " (", Sigma[3], ", ", "n=150", ")" ))
ex4n150p5s3 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab41s3)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) )

xlab42s3 = expression(paste("Example 4", " (", Sigma[3], ", ", "n=300", ")" ))
ex4n300p5s3 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab42s3)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) )

xlab43s3 = expression(paste("Example 4", " (", Sigma[3], ", ", "n=500", ")" ))
ex4n500p5s3 <- ggplot(dat,aes(x = Index, y = Coefficients) ) + ggtitle(xlab43s3)+ xlab("")+
  geom_boxplot()+ geom_point(data = data.frame(x= factor(1:18), y= gamma0),
                             aes(x=x, y=y), color = 'red', shape = 8) +
  theme(plot.title = element_text(hjust = 0.5) )



multiplot( ex2n150p5, ex2n300p5, ex2n500p5 , 
           ex2n150p5s2, ex2n300p5s2, ex2n500p5s2,
           ex2n150p5s3, ex2n300p5s3, ex2n500p5s3, cols=3)

multiplot( ex4n150p5, ex4n300p5, ex4n500p5 , 
           ex4n150p5s2, ex4n300p5s2, ex4n500p5s2,
           ex4n150p5s3, ex4n300p5s3, ex4n500p5s3, cols=3)
