library("KernSmooth")
library("np")

######################################functions
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

#################################Simulation1 
set.seed(18)
n=2000
J=20
x=(1:n)/n
tau=sort(sample(1:49, size=J, replace = FALSE))/50
fi=runif(1, min=0, max=2*pi)
thi=runif(1, min=0, max=2*pi)
beta=sample(c(-1,-0.5,0.5,1), size=J, replace = TRUE)
signal=0.1*(sin(20*pi*x+fi)+2*sin(8*pi*x+thi)) 
for(i in 1:J){
  signal<-signal+beta[i]*(x>tau[i])
}
h1=c(0.01*0.85^2,0.01*0.85,0.01)
lambda=c(0.24,0.22,0.20,0.18)
num1=matrix(rep(0,12),nr=3,ncol=4)       #average number of candidates with theta=0
NUM1=matrix(rep(0,12),nr=3,ncol=4)      #average number of final estimates with theta=0
NUM1SD=matrix(rep(0,12),nr=3,ncol=4)   #sd of number of final estimates with theta=0
FDP1=matrix(rep(0,12),nr=3,ncol=4)        #average false discovery proportion with theta=0
FDP1SD=matrix(rep(0,12),nr=3,ncol=4)        #sd of false discovery proportion with theta=0
num2=matrix(rep(0,12),nr=3,ncol=4)       # theta=0.5  
NUM2=matrix(rep(0,12),nr=3,ncol=4)    
NUM2SD=matrix(rep(0,12),nr=3,ncol=4)    
FDP2=matrix(rep(0,12),nr=3,ncol=4)   
FDP2SD=matrix(rep(0,12),nr=3,ncol=4)            
################################theta=0
for(k in 1:3){
  h=h1[k]
  for(l in 1:4){
    fp=1:100                 #false positives
    fdp=1:100               #false discovery proportion
    n1=1:100                #number of candidates 
    n2=1:100                #number of final change points 
    for(j in 1:100){
      y=signal+rnorm(n,mean=0,sd=0.1)  
      L=abs(screening(y,x,h))
      x1=localMax(L,floor(h*n))/n
      candidate=x1[which(L[localMax(L,floor(h*n))]>lambda[l])]
      x2=c(0,(candidate[1:(length(candidate)-1)]+candidate[2:length(candidate)])/2,1.01)
      pvalue=1:length(candidate)
      for(i in 1:length(candidate)){
        subx=x[which((x>x2[i])&(x<x2[i+1]))]
        suby=y[which((x>x2[i])&(x<x2[i+1]))]
        test=localtest(suby,subx,candidate[i],max(10,0.1*length(subx))/n)
        pvalue[i]=pchisq(test,df=1,lower.tail=FALSE)
      }
      estimate=candidate[which(p.adjust(pvalue, "BH")<0.1)]
      fp[j]=0
      for(i in 1:length(estimate)){
        if(min(abs(estimate[i]-tau))>(2/n))
          fp[j]=fp[j]+1
      }
      n1[j]=length(candidate)
      n2[j]=length(estimate)
      fdp[j]=fp[j]/n2[j]
    }
    FDP1[k,l]=mean(fdp)
    FDP1SD[k,l]=sd(fdp)
    num1[k,l]=mean(n1)
    NUM1[k,l]=mean(n2)
    NUM1SD[k,l]=sd(n2)
  }
}

######################################theta=0.5
for(k in 1:3){
  h=h1[k]
  for(l in 1:4){
    fp=1:100                 #flase positives
    fdp=1:100               #false discovery proportion
    n1=1:100                #number of candidates 
    n2=1:100                #number of final change points 
    for(j in 1:100){
      y=signal+rnorm(n,mean=0,sd=0.1*(1+0.5*sin(2*pi*x)))  
      L=abs(screening(y,x,h))
      x1=localMax(L,floor(h*n))/n
      candidate=x1[which(L[localMax(L,floor(h*n))]>lambda[l])]
      x2=c(0,(candidate[1:(length(candidate)-1)]+candidate[2:length(candidate)])/2,1.01)
      pvalue=1:length(candidate)
      for(i in 1:length(candidate)){
        subx=x[which((x>x2[i])&(x<x2[i+1]))]
        suby=y[which((x>x2[i])&(x<x2[i+1]))]
        test=localtest(suby,subx,candidate[i],max(10,0.1*length(subx))/n)
        pvalue[i]=pchisq(test,df=1,lower.tail=FALSE)
      }
      estimate=candidate[which(p.adjust(pvalue, "BH")<0.1)]
      fp[j]=0
      for(i in 1:length(estimate)){
        if(min(abs(estimate[i]-tau))>(2/n))
          fp[j]=fp[j]+1
      }
      n1[j]=length(candidate)
      n2[j]=length(estimate)
      fdp[j]=fp[j]/n2[j]
    }
    FDP2[k,l]=mean(fdp)
    FDP2SD[k,l]=sd(fdp)
    num2[k,l]=mean(n1)
    NUM2[k,l]=mean(n2)
    NUM2SD[k,l]=sd(n2)
  }
}

#####################################Figure
par(mfrow=c(2,3))
plot(x=lambda, y=num1[1, ], pch=2, col=3, type="b", ylim=c(10,45), xlab = "lambda",ylab ="number of change points",main="h=0.01*0.85^2")
lines(x=lambda, y=NUM1[1, ], col=2, type="b")
legend("topright", legend=c("screening","multiple testing"), pch=c(2,1), col=c(3,2), bty="n")
plot(x=lambda, y=num1[2, ], pch=2, col=3, type="b", ylim=c(10,45), xlab = "lambda",ylab ="number of change points",main="h=0.01*0.85")
lines(x=lambda, y=NUM1[2, ], col=2, type="b")
legend("topright", legend=c("screening","multiple testing"), pch=c(2,1), col=c(3,2), bty="n")
plot(x=lambda, y=num1[3, ], pch=2, col=3, type="b", ylim=c(10,45), xlab = "lambda",ylab ="number of change points",main="h=0.01")
lines(x=lambda, y=NUM1[3, ], col=2, type="b")
legend("topright", legend=c("screening","multiple testing"), pch=c(2,1), col=c(3,2), bty="n")

plot(x=lambda, y=num2[1, ], pch=2, col=3, type="b", ylim=c(10,45), xlab = "lambda",ylab ="number of change points",main="h=0.01*0.85^2")
lines(x=lambda, y=NUM2[1, ], col=2, type="b")
legend("topright", legend=c("screening","multiple testing"), pch=c(2,1), col=c(3,2), bty="n")
plot(x=lambda, y=num2[2, ], pch=2, col=3, type="b", ylim=c(10,45), xlab = "lambda",ylab ="number of change points",main="h=0.01*0.85")
lines(x=lambda, y=NUM2[2, ], col=2, type="b")
legend("topright", legend=c("screening","multiple testing"), pch=c(2,1), col=c(3,2), bty="n")
plot(x=lambda, y=num2[3, ], pch=2, col=3, type="b", ylim=c(10,45), xlab = "lambda",ylab ="number of change points",main="h=0.01")
lines(x=lambda, y=NUM2[3, ], col=2, type="b")
legend("topright", legend=c("screening","multiple testing"), pch=c(2,1), col=c(3,2), bty="n")

######################################Table
rbind(t(rbind(NUM1[1,],NUM1SD[1,],FDP1[1,],FDP1SD[1,],NUM1[2,],NUM1SD[2,],FDP1[2,],FDP1SD[2,],NUM1[3,],NUM1SD[3,],FDP1[3,],FDP1SD[3,])),
      t(rbind(NUM2[1,],NUM2SD[1,],FDP2[1,],FDP2SD[1,],NUM2[2,],NUM2SD[2,],FDP2[2,],FDP2SD[2,],NUM2[3,],NUM2SD[3,],FDP2[3,],FDP2SD[3,])))
