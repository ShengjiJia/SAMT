library("KernSmooth")
library("np")

###############useful functions
localMax <-function(y, span = 10){
  if (length(y) < span * 2 + 1)  return(NULL)
  n  = length(y)
  index = NULL
  for (i in (span + 1) : (n - span) ) {
    if ( y[i] == max(y[(i - span) : (i + span)]) ) index = c(index, i)
  }
  return (index)
}

screening<-function(y, x, h){
  n=length(x)
  rkernal<-function(t){
    0.75*(1-((x-t)/h)^2)*(((x-t)/h)>=0)*(((x-t)/h)<=1)
  }
  lkernal<-function(t){
    0.75*(1-((x-t)/h)^2)*(((x-t)/h)>=-1)*(((x-t)/h)<0)
  }
  rweight<-function(t){
    rkernal(t)*(sum(rkernal(t)*(x-t)^2)-(x-t)*sum(rkernal(t)*(x-t)))
  }
  lweight<-function(t){
    lkernal(t)*(sum(lkernal(t)*(x-t)^2)-(x-t)*sum(lkernal(t)*(x-t)))
  }
  rlimit<-function(t){
    sum(rweight(t)*y)/sum(rweight(t))
  }
  llimit<-function(t){
    sum(lweight(t)*y)/sum(lweight(t))
  }
  right=rep(0,n)
  for (i in 1:(n-1)){
    right[i]=rlimit(x[i])   
  }
  left=rep(0,n)
  for (i in 2:n){
    left[i]=llimit(x[i])
  }
  L=right-left
  for(i in 1:floor(h*n)){
    L[i]=0
  }
  for(i in (n-floor(h*n)):n){
    L[i]=0
  }
  return(L)
}

localtest<-function (y, x, t0, h){
  W=function(t,x,h){diag(1/h*3/4*(1-((x-t)/h)^2)*((1-((x-t)/h)^2)>0))}
  D=function(t,x,h){matrix(c(rep(1,times=length(x)),(x-t)/h),nrow=length(x),ncol=2)}
  s=function(t,x,h){t(c(1,0))%*%solve(t(D(t,x,h))%*%W(t,x,h)%*%D(t,x,h))%*%t(D(t,x,h))%*%W(t,x,h)}
  S=NULL                                                            #smooth matrix S
  for(i in 1:length(x)){S=rbind(S,s(x[i],x,h))}
  z=1*(x>t0)
  Z=(diag(1,length(x))-S)%*%z                             #tildeZ
  Y=(diag(1,length(x))-S)%*%y                             #tildeY
  r=lm(Y~0+Z)$residuals
  fit=npreg(r^2~x,regtype="ll",ckertype="epanechnikov")
  sigma=diag(fit$mean)
  Wald=(t(Z)%*%Y)^2/sum(fit$mean*(Z^2))
  Wald
}

###############Simulation in supplementary material
set.seed(23)
n=100
beta=c(0,0.1,0.2,0.3,0.4,0.5)
test=matrix(0,nr=3,nc=500)
power=matrix(0,nr=3,nc=6)
bandwidth=c(0.09,0.12,0.15)
x=(1:100)/100
for (k in 1:3){
  b=bandwidth[k]
  reject=rep(0,times=500)          #reject or not to calculate power
  for (i in 1:500){
    y=exp(-x)+rnorm(n,mean=0,sd=0.1*(1+x)) 
    test[k,i]=localtest(y,x,0.6,b)                       #under null
    if(localtest(y,x,0.6,b)>qchisq(0.95,df=1)) {reject[i]=1}
  }
  power[k,1]=mean(reject)
  for (j in 2:6){                                             #under alternative
    reject=rep(0,times=500)                 
    for (i in 1:500){
      x=(1:100)/100
      y=exp(-x)+beta[j]*(x>0.6)+rnorm(n,mean=0,sd=0.1*(1+x))
      if(localtest(y,x,0.6,b)>qchisq(0.95,df=1)) {reject[i]=1}
    }
    power[k,j]=mean(reject)
  }
}

par(mfrow=c(2,3))
plot(qchisq((1:99)/100,df=1),quantile(test[1,],probs=(1:99)/100), xlab="Theoritical Quantiles", ylab="Sample Quantiles",main="h=0.09")
lines(qchisq((1:99)/100,df=1),qchisq((1:99)/100,df=1),lwd=2)
plot(qchisq((1:99)/100,df=1),quantile(test[2,],probs=(1:99)/100), xlab="Theoritical Quantiles", ylab="Sample Quantiles",main="h=0.12")
lines(qchisq((1:99)/100,df=1),qchisq((1:99)/100,df=1),lwd=2)
plot(qchisq((1:99)/100,df=1),quantile(test[3,],probs=(1:99)/100), xlab="Theoritical Quantiles", ylab="Sample Quantiles",main="h=0.15")
lines(qchisq((1:99)/100,df=1),qchisq((1:99)/100,df=1),lwd=2)
plot(x=beta,y=power[1,], xlab="beta",ylab="power",ylim=c(0,1),main="h=0.09", type="b")
plot(x=beta,y=power[2,], xlab="beta",ylab="power",ylim=c(0,1),main="h=0.12", type="b")
plot(x=beta,y=power[3,], xlab="beta",ylab="power",ylim=c(0,1),main="h=0.15", type="b")

set.seed(18)
bandwidth=c(0.09,0.12,0.15)
estimate=matrix(0,nr=3,nc=500)
test=matrix(0,nr=3,nc=500)
for (k in 1:3){
  b=bandwidth[k]
  for (l in 1:500){
    y=exp(-x)+beta[6]*(x>0.6)+rnorm(n,mean=0,sd=0.1*(1+x))
    L=abs(screening(y,x,b))
    estimate[k,l]=which.max(L)/n
  }
  for (l in 1:500){
    y=exp(-x)+rnorm(n,mean=0,sd=0.1*(1+x))
    test[k,l]=localtest(y,x,estimate[k,l],b)                       #under null
  }
}

par(mfrow=c(2,3))
hist(estimate[1,], breaks=20, xlab="location", ylab="frequency", main="h=0.09")
hist(estimate[2,], breaks=20, xlab="location", ylab="frequency", main="h=0.12")
hist(estimate[3,], breaks=20, xlab="location", ylab="frequency", main="h=0.15")
plot(qchisq((1:99)/100,df=1),quantile(test[1,],probs=(1:99)/100), xlab="Theoritical Quantiles", ylab="Sample Quantiles",main="h=0.09")
lines(qchisq((1:99)/100,df=1),qchisq((1:99)/100,df=1),lwd=2)
plot(qchisq((1:99)/100,df=1),quantile(test[2,],probs=(1:99)/100), xlab="Theoritical Quantiles", ylab="Sample Quantiles",main="h=0.12")
lines(qchisq((1:99)/100,df=1),qchisq((1:99)/100,df=1),lwd=2)
plot(qchisq((1:99)/100,df=1),quantile(test[3,],probs=(1:99)/100), xlab="Theoritical Quantiles", ylab="Sample Quantiles",main="h=0.15")
lines(qchisq((1:99)/100,df=1),qchisq((1:99)/100,df=1),lwd=2)
