source("cpfunctionscad.R")

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


########### COMPARE WITH OTHER METHODS #########################

########### Jin and Li 2017 annals of statastics ###############

ord=order(ZZ[,6]) ## 4,6
training1=ZZ[ord,]

##estimate thresholds
Y=training1[,1]
X=training1[,-1]
n = nrow(training1)
p = ncol(X)

delta=rep(1,n)
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
#   sapply(tsmc, function(x){x$cp})
# choose the optimal results by BIC
length(tsmcp[[1]]) # number of change points
X[tsmcp[[1]],3] #thresholds.   



## cd40
z1 = X[,3]
s.tsmcp1 = length(tsmcp[[1]]) # change points
a.tsmcp1 = X[tsmcp[[1]],3] #thresholds.   
coef.tsmcp1 = tsmcp[[2]] # coefficients. 
## in-sample
mse.is.tsmcp = tsmcp[[3]] #variance of error. 

## age
z2 = X[,5]
s.tsmcp2 = length(tsmcp[[1]]) # change points
a.tsmcp2 = X[tsmcp[[1]],5] #thresholds.   
coef.tsmcp2 = tsmcp[[2]] # coefficients.




#   #########################
#   #####equal weight
#   ######################

Y=ZZ[,1]
X=ZZ[,-1]
n = nrow(ZZ)
p = ncol(X)

z1=rowMeans(X[,1:5])
o.z=order(z1)
x.tsmcd=X[o.z,]
y.tsmcd=Y[o.z]

c=seq(0.1,1.5,0.1)
m=ceiling(c*sqrt(n))
c=c[which(m>p+1)]
delta=rep(1,n)

bicy=c
tsmc=NULL
for(i in 1:length(c)){
  tsm = TSMCP(y.tsmcd,x.tsmcd,delta,c[i])
  df = sum(tsm[[2]]!=0)
  # bicy[i] = log(n)*(df)*log(df)/3+n*log(tsm[[3]])  
  bicy[i] = log((length(tsm[[1]])+1)*(p))+n*log(tsm[[3]])  
  tsmc[[i]] = tsm
  print(i)
}

tsmcp.a = tsmc[[which(bicy==min(bicy))[1]]]
#   sapply(tsmc, function(x){x$cp})
s.tsmcp.a = length(tsmcp.a[[1]]) # change points
a.tsmcp.a = z1[o.z][tsmcp.a[[1]]] #thresholds. Real threshods (qnorm(0.3),qnorm(0.6))=(-0.5244,0.2533)
coef.tsmcp.a = tsmcp.a[[2]] # coefficients. Real value:(2,1,1,1,1,1,-1,0,0,-1,-1,-1,0,-1,1,0,0,0)

