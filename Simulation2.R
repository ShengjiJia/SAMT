library(KernSmooth)
library(np)
library(wavethresh)
library(cumSeg)
library(DNAcopy)

######################################useful functions
#Estimate the standard deviation of the intensities
estimateSigma<-function (Y, h = 10) {  
  n = length(Y)
  YBar = rep(0, n)
  for (i in 1:n) {
    a = min(n, i + h)
    b = max(1, i - h)
    YBar[i] = mean(Y[b:a])
  }
  return(sqrt(var(Y - YBar) * (2 * h + 1)/(2 * h)))
}

#Calculate the value for local diagnostic function
localDiagnostic<-function (y, h) { 
  yy = c(rep(0, h - 1), y, rep(0, h))
  n = length(y)
  z = rep(0, n)
  for(i in 1:n){
    z[i]=sum(yy[i:(h+i-1)])/h-sum(yy[(h+i):(2*h-1+i)])/h
  }
  return(z)
}

#Get the local maximizers of local diagnostic function
localMax<-function (y, span = 5) {  
  if (length(y) < span * 2 + 1) 
    return(NULL)
  n = length(y)
  index = NULL
  for (i in (span + 1):(n - span)) {
    if (y[i] == max(y[(i - span):(i + span)])) 
      index = c(index, i)
  }
  return(index)
}

clean<-function (LocalM, h) 
{
  len <- length(LocalM)
  rm.list <- NULL
  for (i in 1:(len - 1)) {
    if (LocalM[i] >= LocalM[i + 1] - h) {
      rm.list <- c(rm.list, i)
    }
  }
  if (length(rm.list) > 0) {
    LocalM <- LocalM[-as.vector(rm.list)]
  }
  return(LocalM = LocalM)
}

SARAp<-function (Y, h, hh = 2 * h, sigma = NULL) { 
  n = length(Y)
  LDF = localDiagnostic(Y, h)
  LDF.pos = LDF
  LDF.neg = -LDF
  if (is.null(sigma)) 
    sigma = estimateSigma(Y, h = max(3, 2 * floor(log(n))))
  pV.pos = 1 - 2 * pnorm(LDF.pos/(sqrt(2/h) * sigma))
  LocalMax = localMax(LDF.pos, span = hh)
  LocalMax = clean(LocalMax, h)
  LocalMaxValue = pV.pos[LocalMax]
  pV.neg = 1 - 2 * pnorm(LDF.neg/(sqrt(2/h) * sigma))
  LocalMin = localMax(LDF.neg, span = hh)
  LocalMin = clean(LocalMin, h)
  LocalMinValue = pV.neg[LocalMin]
  LocalExt <- c(LocalMax, LocalMin)
  LocalExtValue <- c(LocalMaxValue, LocalMinValue)
  LocalExtValue <- LocalExtValue[order(LocalExt)]
  LocalExt <- sort(LocalExt)
  return(list(index = LocalExt, pV = LocalExtValue))
}

#Get the inverse cumulative distribution function of local min p-values
fInverse<-function (n = 10000, h = 10, hh = 2 * h, precise = 10000, simT = 100) { 
  empirical = NULL
  for (i in 1:simT) {
    Y = rnorm(n)
    LDF = localDiagnostic(Y, h)
    LDF.pos = LDF
    LDF.neg = -LDF
    sigma = 1
    index.pos = localMax(y = LDF.pos, span = hh)
    pV.pos = 1 - 2 * pnorm(LDF.pos[index.pos]/(sqrt(2/h) * 
                                                 sigma))
    index.neg = localMax(y = LDF.neg, span = hh)
    pV.neg = 1 - 2 * pnorm(LDF.neg[index.neg]/(sqrt(2/h) * 
                                                 sigma))
    index <- c(index.pos, index.neg)
    pv <- c(pV.pos, pV.neg)
    pv <- pv[order(index)]
    index <- sort(index)
    len <- length(index)
    rm.list <- NULL
    for (j in 1:(len - 1)) {
      if (index[j] >= index[j + 1] - h) {
        rm.list <- c(rm.list, j)
      }
    }
    if (length(rm.list) > 0) {
      pv <- pv[-rm.list]
    }
    empirical <- c(empirical, pv)
    if (length(empirical) > 10 * precise) 
      break
  }
  return(quantile(empirical, probs = c(0:precise)/precise))
}

SARA<-function (Y, h = 10, hh = 2 * h, FINV = NULL, sigma = NULL, precise = 10000) {
  object = SARAp(Y = Y, h = h, hh = hh, sigma = sigma)
  index = object$index
  pV = object$pV
  if (is.null(FINV)) 
    FINV = fInverse(n = length(Y), h = h, hh = hh, precise = precise, 
                    simT = 100)
  pVcorr = pV
  for (i in 1:length(pV)) {
    pVcorr[i] = (length(which(FINV < pV[i])))/(precise + 
                                                 1)
  }
  return(list(index = index, pV = pVcorr))
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

#################################Simulation2
set.seed(9)
n=2048
J=20
x=(1:n)/n
h=0.01
lambda=0.2
################################case1
method1=matrix(0,nr=10,ncol=100)          
method2=matrix(0,nr=10,ncol=100)        
method3=matrix(0,nr=10,ncol=100)  
method4=matrix(0,nr=10,ncol=100) 
method5=matrix(0,nr=10,ncol=100)        
method6=matrix(0,nr=10,ncol=100)  
method7=matrix(0,nr=10,ncol=100) 
n1=1:100
n2=1:100
n3=1:100
n4=1:100
n5=1:100
n6=1:100
n7=1:100
for(j in 1:100){
  tau=(((1:J)-1)*100+sample(1:100,size=J,replace=TRUE))/n
  beta1=sample(c(-0.5,0.5), size=J, replace = TRUE)
  beta2=sample(c(-0.5,0.5), size=J, replace = TRUE)
  beta3=sample(c(-0.5,0.5), size=J, replace = TRUE)
  beta4=sample(c(-0.5,0.5), size=J, replace = TRUE)
  signal1=0.1*(sin(20*pi*x+runif(1,min=0,max=2*pi))+2*sin(50*pi*x+runif(1,min=0,max=2*pi))) 
  for(i in 1:J) {signal1<-signal1+beta1[i]*(x>tau[i])}
  signal2=0.1*(sin(20*pi*x+runif(1,min=0,max=2*pi))+2*sin(50*pi*x+runif(1,min=0,max=2*pi))) 
  for(i in 1:J) {signal2<-signal2+beta2[i]*(x>tau[i])}
  signal3=0.1*(sin(20*pi*x+runif(1,min=0,max=2*pi))+2*sin(50*pi*x+runif(1,min=0,max=2*pi))) 
  for(i in 1:J) {signal3<-signal3+beta3[i]*(x>tau[i])}
  signal4=0.1*(sin(20*pi*x+runif(1,min=0,max=2*pi))+2*sin(50*pi*x+runif(1,min=0,max=2*pi))) 
  for(i in 1:J) {signal4<-signal4+beta4[i]*(x>tau[i])}
  y1=signal1+rnorm(n,mean=0,sd=0.1*(1+runif(1,min=0,max=0.5)*sin(2*pi*x)))  
  y2=signal2+rnorm(n,mean=0,sd=0.1*(1+runif(1,min=0,max=0.5)*sin(2*pi*x))) 
  y3=signal3+rnorm(n,mean=0,sd=0.1*(1+runif(1,min=0,max=0.5)*sin(2*pi*x))) 
  y4=signal4+rnorm(n,mean=0,sd=0.1*(1+runif(1,min=0,max=0.5)*sin(2*pi*x))) 
  L1=abs(screening(y1,x,h))
  L2=abs(screening(y2,x,h))
  L3=abs(screening(y3,x,h))
  L4=abs(screening(y4,x,h))
  ##############################method1
  CBS=DNAcopy::segment(CNA(y1, rep(1,n), 1:n))
  estimate1=CBS$output[,4]
  estimate1=estimate1[-length(CBS$output[,4])]/n
  for(i in 1:10){
    if(min(abs(tau[2*i]-estimate1))>(2/n))
      method1[i,j]=0
    else
      method1[i,j]=1
  }
  n1[j]=length(estimate1)
  ##############################method2
  estimate2=jumpoints(y1,k=60,output="2")$psi/n
  for(i in 1:10){
    if(min(abs(tau[2*i]-estimate2))>(2/n))
      method2[i,j]=0
    else
      method2[i,j]=1
  }
  n2[j]=length(estimate2)
  ##############################method3
  sara=SARA(y1,h=20)
  estimate3=sara$index[which(p.adjust(sara$pV, "BH")<0.1)]/n
  for(i in 1:10){
    if(min(abs(tau[2*i]-estimate3))>(2/n))
      method3[i,j]=0
    else
      method3[i,j]=1
  }
  n3[j]=length(estimate3)
  ##############################method4
  wavelet=wd(y1, filter.number=1, family="DaubExPhase")
  threshold=threshold(wavelet, levels =1:(nlevelsWT(wavelet)-1), type = "hard", policy = "universal",dev = madmad)
  estimate4=c(which(accessD(threshold,level=7)!=0)/2^7-2^(-8), which(accessD(threshold,level=8)!=0)/2^8-2^(-9),
              which(accessD(threshold,level=9)!=0)/2^9-2^(-10))
  for(i in 1:10){
    if(min(abs(tau[2*i]-estimate4))>(2/n))
      method4[i,j]=0
    else
      method4[i,j]=1
  }
  n4[j]=length(estimate4)
  ##############################method5
  L=L1
  x1=localMax(L,floor(h*n))/n
  candidate=x1[which(L[localMax(L,floor(h*n))]>lambda)]
  x2=c(0,(candidate[1:(length(candidate)-1)]+candidate[2:length(candidate)])/2,1.01)
  pvalue=1:length(candidate)
  for(i in 1:length(candidate)){
    subx=x[which((x>x2[i])&(x<x2[i+1]))]
    suby=y1[which((x>x2[i])&(x<x2[i+1]))]
    test=localtest(suby,subx,candidate[i],max(10,0.1*length(subx))/n)
    pvalue[i]=pchisq(test,df=1,lower.tail=FALSE)
  }
  estimate5=candidate[which(p.adjust(pvalue, "BH")<0.1)]
  for(i in 1:10){
    if(min(abs(tau[2*i]-estimate5))>(2/n))
      method5[i,j]=0
    else
      method5[i,j]=1
  }
  n5[j]=length(estimate5)
  ##########################method6&7
  L=L1^2+L2^2+L3^2+L4^2
  x1=localMax(L,floor(h*n))/n
  candidate=x1[which(L[localMax(L,floor(h*n))]>lambda)]
  x2=c(0,(candidate[1:(length(candidate)-1)]+candidate[2:length(candidate)])/2,1.01)
  pvalue=1:length(candidate)
  pSIM=1:length(candidate)
  for(i in 1:length(candidate)){
    subx=x[which((x>x2[i])&(x<x2[i+1]))]
    suby1=y1[which((x>x2[i])&(x<x2[i+1]))]
    suby2=y2[which((x>x2[i])&(x<x2[i+1]))]
    suby3=y3[which((x>x2[i])&(x<x2[i+1]))]
    suby4=y4[which((x>x2[i])&(x<x2[i+1]))]
    test1=localtest(suby1,subx,candidate[i],max(10,0.1*length(subx))/n)
    test2=localtest(suby2,subx,candidate[i],max(10,0.1*length(subx))/n)
    test3=localtest(suby3,subx,candidate[i],max(10,0.1*length(subx))/n)
    test4=localtest(suby4,subx,candidate[i],max(10,0.1*length(subx))/n)
    test=test1+test2+test3+test4
    pvalue[i]=pchisq(test,df=4,lower.tail=FALSE)
    p1=pchisq(test1,df=1,lower.tail=FALSE)
    p2=pchisq(test2,df=1,lower.tail=FALSE)
    p3=pchisq(test3,df=1,lower.tail=FALSE)
    p4=pchisq(test4,df=1,lower.tail=FALSE)
    pSIM[i]=pnorm(0.5*qnorm(p1)+0.5*qnorm(p2)+0.5*qnorm(p3)+0.5*qnorm(p4))
  }
  estimate6=candidate[which(p.adjust(pvalue, "BH")<0.1)]
  estimate7=candidate[which(p.adjust(pSIM, "BH")<0.1)]
  for(i in 1:10){
    if(min(abs(tau[2*i]-estimate6))>(2/n))
      method6[i,j]=0
    else
      method6[i,j]=1
  }
  for(i in 1:10){
    if(min(abs(tau[2*i]-estimate7))>(2/n))
      method7[i,j]=0
    else
      method7[i,j]=1
  }
  n6[j]=length(estimate6)
  n7[j]=length(estimate7)
}

################################case2
Method1=matrix(0,nr=10,ncol=100)          
Method2=matrix(0,nr=10,ncol=100)        
Method3=matrix(0,nr=10,ncol=100)
Method4=matrix(0,nr=10,ncol=100)  
Method5=matrix(0,nr=10,ncol=100)        
Method6=matrix(0,nr=10,ncol=100)
Method7=matrix(0,nr=10,ncol=100)     
N1=1:100
N2=1:100
N3=1:100
N4=1:100
N5=1:100
N6=1:100
N7=1:100
for(j in 1:100){
  tau=(((1:J)-1)*100+sample(1:100,size=J,replace=TRUE))/n
  beta1=sample(c(-1,-0.5,0,0.5,1), size=J, replace = TRUE)
  beta2=sample(c(-1,-0.5,0,0.5,1), size=J, replace = TRUE)
  beta3=sample(c(-1,-0.5,0,0.5,1), size=J, replace = TRUE)
  beta4=sample(c(-1,-0.5,0,0.5,1), size=J, replace = TRUE)
  signal1=0.1*(sin(20*pi*x+runif(1,min=0,max=2*pi))+2*sin(50*pi*x+runif(1,min=0,max=2*pi))) 
  for(i in 1:J) {signal1<-signal1+beta1[i]*(x>tau[i])}
  signal2=0.1*(sin(20*pi*x+runif(1,min=0,max=2*pi))+2*sin(50*pi*x+runif(1,min=0,max=2*pi))) 
  for(i in 1:J) {signal2<-signal2+beta2[i]*(x>tau[i])}
  signal3=0.1*(sin(20*pi*x+runif(1,min=0,max=2*pi))+2*sin(50*pi*x+runif(1,min=0,max=2*pi))) 
  for(i in 1:J) {signal3<-signal3+beta3[i]*(x>tau[i])}
  signal4=0.1*(sin(20*pi*x+runif(1,min=0,max=2*pi))+2*sin(50*pi*x+runif(1,min=0,max=2*pi))) 
  for(i in 1:J) {signal4<-signal4+beta4[i]*(x>tau[i])}
  y1=signal1+rnorm(n,mean=0,sd=0.1*(1+runif(1,min=0,max=0.5)*sin(2*pi*x)))  
  y2=signal2+rnorm(n,mean=0,sd=0.1*(1+runif(1,min=0,max=0.5)*sin(2*pi*x))) 
  y3=signal3+rnorm(n,mean=0,sd=0.1*(1+runif(1,min=0,max=0.5)*sin(2*pi*x))) 
  y4=signal4+rnorm(n,mean=0,sd=0.1*(1+runif(1,min=0,max=0.5)*sin(2*pi*x))) 
  L1=abs(screening(y1,x,h))
  L2=abs(screening(y2,x,h))
  L3=abs(screening(y3,x,h))
  L4=abs(screening(y4,x,h))
  ##############################method1
  CBS=DNAcopy::segment(CNA(y1, rep(1,n), 1:n))
  estimate1=CBS$output[,4]
  estimate1=estimate1[-length(CBS$output[,4])]/n
  for(i in 1:10){
    if(min(abs(tau[2*i]-estimate1))>(2/n))
      Method1[i,j]=0
    else
      Method1[i,j]=1
  }
  N1[j]=length(estimate1)
  ##############################method2
  estimate2=jumpoints(y1,k=60,output="2")$psi/n
  for(i in 1:10){
    if(min(abs(tau[2*i]-estimate2))>(2/n))
      Method2[i,j]=0
    else
      Method2[i,j]=1
  }
  N2[j]=length(estimate2)
  ##############################method3
  sara=SARA(y1,h=20)
  estimate3=sara$index[which(p.adjust(sara$pV, "BH")<0.1)]/n
  for(i in 1:10){
    if(min(abs(tau[2*i]-estimate3))>(2/n))
      Method3[i,j]=0
    else
      Method3[i,j]=1
  }
  N3[j]=length(estimate3)
  ##############################method4
  wavelet=wd(y1, filter.number=1, family="DaubExPhase")
  threshold=threshold(wavelet, levels =1:(nlevelsWT(wavelet)-1), type = "hard", policy = "universal",dev = madmad)
  estimate4=c(which(accessD(threshold,level=7)!=0)/2^7-2^(-8), which(accessD(threshold,level=8)!=0)/2^8-2^(-9),
              which(accessD(threshold,level=9)!=0)/2^9-2^(-10))
  for(i in 1:10){
    if(min(abs(tau[2*i]-estimate4))>(2/n))
      Method4[i,j]=0
    else
      Method4[i,j]=1
  }
  N4[j]=length(estimate4)
  ##############################method5
  L=L1
  x1=localMax(L,floor(h*n))/n
  candidate=x1[which(L[localMax(L,floor(h*n))]>lambda)]
  x2=c(0,(candidate[1:(length(candidate)-1)]+candidate[2:length(candidate)])/2,1.01)
  pvalue=1:length(candidate)
  for(i in 1:length(candidate)){
    subx=x[which((x>x2[i])&(x<x2[i+1]))]
    suby=y1[which((x>x2[i])&(x<x2[i+1]))]
    test=localtest(suby,subx,candidate[i],max(10,0.1*length(subx))/n)
    pvalue[i]=pchisq(test,df=1,lower.tail=FALSE)
  }
  estimate5=candidate[which(p.adjust(pvalue, "BH")<0.1)]
  for(i in 1:10){
    if(min(abs(tau[2*i]-estimate5))>(2/n))
      Method5[i,j]=0
    else
      Method5[i,j]=1
  }
  N5[j]=length(estimate5)
  ##########################method6&7
  L=L1^2+L2^2+L3^2+L4^2
  x1=localMax(L,floor(h*n))/n
  candidate=x1[which(L[localMax(L,floor(h*n))]>lambda)]
  x2=c(0,(candidate[1:(length(candidate)-1)]+candidate[2:length(candidate)])/2,1.01)
  pvalue=1:length(candidate)
  pSIM=1:length(candidate)
  for(i in 1:length(candidate)){
    subx=x[which((x>x2[i])&(x<x2[i+1]))]
    suby1=y1[which((x>x2[i])&(x<x2[i+1]))]
    suby2=y2[which((x>x2[i])&(x<x2[i+1]))]
    suby3=y3[which((x>x2[i])&(x<x2[i+1]))]
    suby4=y4[which((x>x2[i])&(x<x2[i+1]))]
    test1=localtest(suby1,subx,candidate[i],max(10,0.1*length(subx))/n)
    test2=localtest(suby2,subx,candidate[i],max(10,0.1*length(subx))/n)
    test3=localtest(suby3,subx,candidate[i],max(10,0.1*length(subx))/n)
    test4=localtest(suby4,subx,candidate[i],max(10,0.1*length(subx))/n)
    test=test1+test2+test3+test4
    pvalue[i]=pchisq(test,df=4,lower.tail=FALSE)
    p1=pchisq(test1,df=1,lower.tail=FALSE)
    p2=pchisq(test2,df=1,lower.tail=FALSE)
    p3=pchisq(test3,df=1,lower.tail=FALSE)
    p4=pchisq(test4,df=1,lower.tail=FALSE)
    pSIM[i]=pnorm(0.5*qnorm(p1)+0.5*qnorm(p2)+0.5*qnorm(p3)+0.5*qnorm(p4))
  }
  estimate6=candidate[which(p.adjust(pvalue, "BH")<0.1)]
  estimate7=candidate[which(p.adjust(pSIM, "BH")<0.1)]
  for(i in 1:10){
    if(min(abs(tau[2*i]-estimate6))>(2/n))
      Method6[i,j]=0
    else
      Method6[i,j]=1
  }
  for(i in 1:10){
    if(min(abs(tau[2*i]-estimate7))>(2/n))
      Method7[i,j]=0
    else
      Method7[i,j]=1
  }
  N6[j]=length(estimate6)
  N7[j]=length(estimate7)
}

#####show table
rbind(c(mean(n1),apply(method1,1,mean)),
      c(mean(n2),apply(method2,1,mean)),
      c(mean(n3),apply(method3,1,mean)),
      c(mean(n4),apply(method4,1,mean)),
      c(mean(n5),apply(method5,1,mean)),
      c(mean(n6),apply(method6,1,mean)),
      c(mean(n7),apply(method7,1,mean)),
      c(mean(N1),apply(Method1,1,mean)),
      c(mean(N2),apply(Method2,1,mean)),
      c(mean(N3),apply(Method3,1,mean)),
      c(mean(N4),apply(Method4,1,mean)),
      c(mean(N5),apply(Method5,1,mean)),
      c(mean(N6),apply(Method6,1,mean)),
      c(mean(N7),apply(Method7,1,mean)))
