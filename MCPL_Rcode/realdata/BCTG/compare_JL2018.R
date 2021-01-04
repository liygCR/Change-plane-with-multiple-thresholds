#### 
data <- read.csv("mrss.csv")
data <- na.omit(data)
dat <- cbind(mrss=data$mrss, data[, -c(1:2,14)])

X = model.matrix(~ haq+pga+dlcop+fvcp+over+pain+fev1p+ethnic+sex+durdis+age, data = dat)


ZZ <- cbind(scale(dat[,1]), apply(X[, -c(1,9,10)], 2, scale), X[,9:10] )
colnames(ZZ)[1] = "mrss"

ZZ <- as.data.frame(ZZ)


##################### compare with other methods ##############
# ZZ <- as.data.frame(ZZ)
########## Jin and Li 2017 annals of statastics ###############
source("cpfunctionscad.R")

ord=order(ZZ[,2]) ## 2-4
ZZ=ZZ[ord,]

##estimate thresholds
Y=ZZ[,1]
X=ZZ[,-1]
n = nrow(ZZ)
delta=rep(1,n)
p = ncol(X)

n = nrow(ZZ)
c=seq(0.1,1.5,0.1)
m=ceiling(c*sqrt(n))
c=c[which(m>p+1)]




bicy=c
tsmc=NULL
for(i in 1:length(c)){
  
  tsm=TSMCP(Y,X,delta,c[i])
  bicy[i]=log((length(tsm[[1]])+1)*(p+1))+n*log(tsm[[3]])
  #   bicy[i]=log(log(n))* log((length(tsm[[1]])+1)*(p+1))+n*log(tsm[[3]])
  
  tsmc[[i]]=tsm
  print(i)
}

tsmcp=tsmc[[which(bicy==min(bicy))[1]]]
sapply(tsmc, function(x){x$cp})
#choose the optimal results by BIC
tsmcp[[1]] # change points
tsmcp[[3]] #variance of error. 

## haq
z1 = X[,1]
s.tsmcp1 = length(tsmcp[[1]]) # change points
a.tsmcp1 = X[tsmcp[[1]],1] #thresholds.   
coef.tsmcp1 = tsmcp[[2]] # coefficients. 

## pga
z2 = X[,2]
s.tsmcp2 = length(tsmcp[[1]]) # change points
a.tsmcp2 = X[tsmcp[[1]],2] #thresholds.   
coef.tsmcp2 = tsmcp[[2]] # coefficients.

## dlcop
z3 = X[,3]
s.tsmcp3 = length(tsmcp[[1]]) # change points
a.tsmcp3 = X[tsmcp[[1]],3] #thresholds.   
coef.tsmcp3 = tsmcp[[2]] # coefficients.


#########################
#####equal weight
######################

delta=rep(1,n)

z1=rowMeans(ZZ[,c(2:4)])
o.z=order(z1)
x.tsmcd=X[o.z,]
y.tsmcd=Y[o.z]
p=ncol(x.tsmcd)
c=seq(0.5,1.5,0.1)
m=ceiling(c*sqrt(n))
c=c[which(m>p+1)]

bicy=c
tsmc=NULL
for(i in 1:length(c)){
  tsm = TSMCP(y.tsmcd,x.tsmcd,delta,c[i])
  #   df = sum(tsm[[2]]!=0)
  #   bicy[i] = log(n)*(df)*log(df)/3+n*log(tsm[[3]])  
  bicy[i]=log((length(tsm[[1]])+1)*(p+1))+n*log(tsm[[3]])
  tsmc[[i]] = tsm
  print(i)
}
tsmcp.ew = tsmc[[which(bicy==min(bicy))[1]]]
sapply(tsmc, function(x){x$cp})
tsmcp.ew[[1]] # change points
z1[o.z][tsmcp.ew[[1]]] #thresholds. 
tsmcp.ew[[2]] # coefficients.
tsmcp.ew[[3]] #variance of error. 

z.ew = z1[o.z]
s.tsmcp.a = length(tsmcp.ew[[1]]) # change points
a.tsmcp.a = z1[o.z][tsmcp.ew[[1]]]  #thresholds.   
coef.tsmcp.a = tsmcp.ew[[2]] # coefficients.
