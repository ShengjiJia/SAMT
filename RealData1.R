library(imputeTS)
library(cumSeg)
library(DNAcopy)
library(KernSmooth)
library(np)
library(wavethresh)

###############################################real data 
CGHdata <- read.csv("C:/Users/Acer/Desktop/My Document/Documents/Research/Projects/change points(SaMT)/CGHdataset.csv", header=TRUE)
index=c(6,19,20,21)
data=CGHdata[2:2301,(1+3*index)]               
n=nrow(data)
d=ncol(data)
x=(1:n)/n
for (i in 1:d){
  data[,i]=na_ma(data[,i],k=5,weighting="linear")        #imputation
}

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
  h3=dpill(x, r^2)
  fit=npreg(r^2~x,bws=h3,ckertype="epanechnikov")
  sigma=diag(fit$mean)
  Wald=(t(Z)%*%Y)^2/sum(fit$mean*(Z^2))
  Wald
}

optbandwidth<-function (y, x, h1, h2){
  nn=length(x)
  n1=length(h1)
  n2=length(h2)
  x_odd=x[2*(1:floor(nn/2))-1]
  x_even=x[2*(1:floor(nn/2))]
  y_odd=y[2*(1:floor(nn/2))-1]
  y_even=y[2*(1:floor(nn/2))]
  RSS=seq(from=0,to=0,length.out=n1*n2)
  hh1=as.vector(matrix(rep(h1,n2),nr=n2, byrow=TRUE))
  hh2=rep(h2,n1)
  for(j in 1:(n1*n2)){
    h=hh1[j]
    L=abs(screening(y_even,x_even,h))
    lambda=4*mad(L)
    n=length(x_even)
    x1=localMax(L,floor(h*n))/n
    candidate=x1[which(L[localMax(L,floor(h*n))]>lambda)]
    x2=c(0,(candidate[1:(length(candidate)-1)]+candidate[2:length(candidate)])/2,1.01)
    pvalue=1:length(candidate)
    for(i in 1:length(candidate)){
      subx=x[which((x>x2[i])&(x<x2[i+1]))]
      suby=y1[which((x>x2[i])&(x<x2[i+1]))]
      test=localtest(suby,subx,candidate[i],max(10,hh2[j]*length(subx))/n)
      pvalue[i]=pchisq(test,df=1,lower.tail=FALSE)
    }
    estimate1=candidate[which(p.adjust(pvalue, "BH")<0.1)]
    jump=0 
    for(i in 1:length(estimate1)){
      jump<-jump+L[(estimate1[i]*n)]*(x_even>estimate1[i])
    }
    h3=dpill(x_even,(y_even-jump))
    fitted=locpoly(x_even,(y_even-jump),bandwidth=h3,gridsize=nn+1,range.x=c(0,1))$y[2*(1:floor(nn/2))]+jump
    RSS[j]=sum(y_odd-fitted)^2
  }
  return(c(hh1[which.min(RSS)],hh2[which.min(RSS)]))
}

##############################################1st sequence  
y1=data[,1]
num1=1:5
#####CBS
set.seed(2)
CBS=DNAcopy::segment(CNA(y1, rep(1,n), 1:n))                
num1[1]=length(CBS$output[,4])-1
######cumSeg
num1[2]=jumpoints(y1,k=60,output="2")$n.psi
######SaRa                 
num1[3]=sum(p.adjust(SARA(y1,h=15)$pV, "BH")<0.1)
######Wavelet
yy=c(y1,y1[2300-(1:1796)])
wavelet=wd(yy, filter.number=1, family="DaubExPhase")
threshold=threshold(wavelet, levels =1:(nlevelsWT(wavelet)-1), type = "hard", policy = "universal",dev = madmad)
a=4096*c(which(accessD(threshold,level=8)!=0)/2^8-2^(-9), which(accessD(threshold,level=9)!=0)/2^9-2^(-10), which(accessD(threshold,level=10)!=0)/2^10-2^(-11))
num1[4]=sum(diff(sort(a*(a<2048)+(4096-a)*(a>2048)))>16)+1
#######SaMT
bandwidth1=optbandwidth(y1,x,h1=c(10,15,20,25,30)/n,h2=c(0.05,0.1,0.15))
h=bandwidth1[1]
L=abs(screening(y1,x,h))
lambda=4*mad(L)
x1=localMax(L,floor(h*n))/n
candidate=x1[which(L[localMax(L,floor(h*n))]>lambda)]
x2=c(0,(candidate[1:(length(candidate)-1)]+candidate[2:length(candidate)])/2,1.01)
pvalue=1:length(candidate)
for(i in 1:length(candidate)){
  subx=x[which((x>x2[i])&(x<x2[i+1]))]
  suby=y1[which((x>x2[i])&(x<x2[i+1]))]
  test=localtest(suby,subx,candidate[i],max(10,bandwidth1[2]*length(subx))/n)
  pvalue[i]=pchisq(test,df=1,lower.tail=FALSE)
}
estimate1=candidate[which(p.adjust(pvalue, "BH")<0.1)]
num1[5]=length(estimate1)
######adaptive Neyman test 
index=CBS$output[,4][1:num1[1]]
Index=c(0,index,n)
fit=NULL
for(i in 1:(length(index)+1)){
  m=mean(y1[(Index[i]+1):Index[i+1]])
  fit=c(fit,rep(m, Index[i+1]-Index[i]))
}
r1=y1-fit
r1=r1-mean(r1)
fft1=NULL
for(j in 1:n/2){
  a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r1)
  b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r1)
  fft1=c(fft1, a, b)
}
v1=var(fft1[floor(n/4):n])        #sigma^2
T1=NULL
for(m in 1:n){
  t=(2*m*v1^2)^(-0.5)*sum((fft1[1:m])^2-v1)
  T1=c(T1,t)
}
T1n=sqrt(2*log(log(n)))*max(T1)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi) 

##############################################2nd sequence  
y2=data[,2]
num2=1:5
#####CBS
set.seed(2)
CBS=DNAcopy::segment(CNA(y2, rep(1,n), 1:n))                
num2[1]=length(CBS$output[,4])-1
######cumSeg
num2[2]=jumpoints(y2,k=60,output="2")$n.psi
######SaRa                 
num2[3]=sum(p.adjust(SARA(y2,h=15)$pV, "BH")<0.1)
######Wavelet
yy=c(y2,y2[2300-(1:1796)])
wavelet=wd(yy, filter.number=1, family="DaubExPhase")
threshold=threshold(wavelet, levels =1:(nlevelsWT(wavelet)-1), type = "hard", policy = "universal",dev = madmad)
a=4096*c(which(accessD(threshold,level=8)!=0)/2^8-2^(-9), which(accessD(threshold,level=9)!=0)/2^9-2^(-10), which(accessD(threshold,level=10)!=0)/2^10-2^(-11))
num2[4]=sum(diff(sort(a*(a<2048)+(4096-a)*(a>2048)))>16)+1
#######SaMT
bandwidth2=optbandwidth(y2,x,h1=c(10,15,20,25,30)/n,h2=c(0.05,0.1,0.15))
h=bandwidth2[1]
L=abs(screening(y2,x,h))
lambda=4*mad(L)
x1=localMax(L,floor(h*n))/n
candidate=x1[which(L[localMax(L,floor(h*n))]>lambda)]
x2=c(0,(candidate[1:(length(candidate)-1)]+candidate[2:length(candidate)])/2,1.01)
pvalue=1:length(candidate)
for(i in 1:length(candidate)){
  subx=x[which((x>x2[i])&(x<x2[i+1]))]
  suby=y2[which((x>x2[i])&(x<x2[i+1]))]
  test=localtest(suby,subx,candidate[i],max(10,bandwidth2[2]*length(subx))/n)
  pvalue[i]=pchisq(test,df=1,lower.tail=FALSE)
}
estimate2=candidate[which(p.adjust(pvalue, "BH")<0.1)]
num2[5]=length(estimate2)
######adaptive Neyman test 
index=CBS$output[,4][1:num2[1]]
Index=c(0,index,n)
fit=NULL
for(i in 1:(length(index)+1)){
  m=mean(y2[(Index[i]+1):Index[i+1]])
  fit=c(fit,rep(m, Index[i+1]-Index[i]))
}
r1=y2-fit
r1=r1-mean(r1)
fft1=NULL
for(j in 1:n/2){
  a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r1)
  b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r1)
  fft1=c(fft1, a, b)
}
v1=var(fft1[floor(n/4):n])        #sigma^2
T1=NULL
for(m in 1:n){
  t=(2*m*v1^2)^(-0.5)*sum((fft1[1:m])^2-v1)
  T1=c(T1,t)
}
T2n=sqrt(2*log(log(n)))*max(T1)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi) 

##############################################3rd sequence  
y3=data[,3]
num3=1:5
#####CBS
set.seed(2)
CBS=DNAcopy::segment(CNA(y3, rep(1,n), 1:n))                
num3[1]=length(CBS$output[,4])-1
######cumSeg
num3[2]=jumpoints(y3,k=60,output="2")$n.psi
######SaRa                 
num3[3]=sum(p.adjust(SARA(y3,h=15)$pV, "BH")<0.1)
######Wavelet
yy=c(y3,y3[2300-(1:1796)])
wavelet=wd(yy, filter.number=1, family="DaubExPhase")
threshold=threshold(wavelet, levels =1:(nlevelsWT(wavelet)-1), type = "hard", policy = "universal",dev = madmad)
a=4096*c(which(accessD(threshold,level=8)!=0)/2^8-2^(-9), which(accessD(threshold,level=9)!=0)/2^9-2^(-10), which(accessD(threshold,level=10)!=0)/2^10-2^(-11))
num3[4]=sum(diff(sort(a*(a<2048)+(4096-a)*(a>2048)))>16)+1
#######SaMT
bandwidth3=optbandwidth(y3,x,h1=c(10,15,20,25,30)/n,h2=c(0.05,0.1,0.15))
h=bandwidth3[1]
L=abs(screening(y3,x,h))
lambda=4*mad(L)
x1=localMax(L,floor(h*n))/n
candidate=x1[which(L[localMax(L,floor(h*n))]>lambda)]
x2=c(0,(candidate[1:(length(candidate)-1)]+candidate[2:length(candidate)])/2,1.01)
pvalue=1:length(candidate)
for(i in 1:length(candidate)){
  subx=x[which((x>x2[i])&(x<x2[i+1]))]
  suby=y3[which((x>x2[i])&(x<x2[i+1]))]
  test=localtest(suby,subx,candidate[i],max(10,bandwidth3[2]*length(subx))/n)
  pvalue[i]=pchisq(test,df=1,lower.tail=FALSE)
}
estimate3=candidate[which(p.adjust(pvalue, "BH")<0.1)]
num3[5]=length(estimate3)
######adaptive Neyman test 
index=CBS$output[,4][1:num3[1]]
Index=c(0,index,n)
fit=NULL
for(i in 1:(length(index)+1)){
  m=mean(y3[(Index[i]+1):Index[i+1]])
  fit=c(fit,rep(m, Index[i+1]-Index[i]))
}
r1=y3-fit
r1=r1-mean(r1)
fft1=NULL
for(j in 1:n/2){
  a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r1)
  b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r1)
  fft1=c(fft1, a, b)
}
v1=var(fft1[floor(n/4):n])        #sigma^2
T1=NULL
for(m in 1:n){
  t=(2*m*v1^2)^(-0.5)*sum((fft1[1:m])^2-v1)
  T1=c(T1,t)
}
T3n=sqrt(2*log(log(n)))*max(T1)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi) 

##############################################4th sequence  
y4=data[,4]
num4=1:5
#####CBS
set.seed(2)
CBS=DNAcopy::segment(CNA(y4, rep(1,n), 1:n))                
num4[1]=length(CBS$output[,4])-1
######cumSeg
num4[2]=jumpoints(y4,k=60,output="2")$n.psi
######SaRa                 
num4[3]=sum(p.adjust(SARA(y4,h=15)$pV, "BH")<0.1)
######Wavelet
yy=c(y4,y4[2300-(1:1796)])
wavelet=wd(yy, filter.number=1, family="DaubExPhase")
threshold=threshold(wavelet, levels =1:(nlevelsWT(wavelet)-1), type = "hard", policy = "universal",dev = madmad)
a=4096*c(which(accessD(threshold,level=8)!=0)/2^8-2^(-9), which(accessD(threshold,level=9)!=0)/2^9-2^(-10), which(accessD(threshold,level=10)!=0)/2^10-2^(-11))
num4[4]=sum(diff(sort(a*(a<2048)+(4096-a)*(a>2048)))>16)+1
#######SaMT
bandwidth4=optbandwidth(y4,x,h1=c(10,15,20,25,30)/n,h2=c(0.05,0.1,0.15))
h=bandwidth4[1]
L=abs(screening(y4,x,h))
lambda=4*mad(L)
x1=localMax(L,floor(h*n))/n
candidate=x1[which(L[localMax(L,floor(h*n))]>lambda)]
x2=c(0,(candidate[1:(length(candidate)-1)]+candidate[2:length(candidate)])/2,1.01)
pvalue=1:length(candidate)
for(i in 1:length(candidate)){
  subx=x[which((x>x2[i])&(x<x2[i+1]))]
  suby=y4[which((x>x2[i])&(x<x2[i+1]))]
  test=localtest(suby,subx,candidate[i],max(10,bandwidth4[2]*length(subx))/n)
  pvalue[i]=pchisq(test,df=1,lower.tail=FALSE)
}
estimate4=candidate[which(p.adjust(pvalue, "BH")<0.1)]
num4[5]=length(estimate4)
######adaptive Neyman test 
index=CBS$output[,4][1:num4[1]]
Index=c(0,index,n)
fit=NULL
for(i in 1:(length(index)+1)){
  m=mean(y4[(Index[i]+1):Index[i+1]])
  fit=c(fit,rep(m, Index[i+1]-Index[i]))
}
r1=y4-fit
r1=r1-mean(r1)
fft1=NULL
for(j in 1:n/2){
  a=sqrt(2/n)*sum(cos(2*pi*j*(1:n)/n)*r1)
  b=sqrt(2/n)*sum(sin(2*pi*j*(1:n)/n)*r1)
  fft1=c(fft1, a, b)
}
v1=var(fft1[floor(n/4):n])        #sigma^2
T1=NULL
for(m in 1:n){
  t=(2*m*v1^2)^(-0.5)*sum((fft1[1:m])^2-v1)
  T1=c(T1,t)
}
T4n=sqrt(2*log(log(n)))*max(T1)-2*log(log(n))-0.5*log(log(log(n)))+0.5*log(4*pi) 

###################################################common change points
h=min(c(bandwidth1[1],bandwidth2[1],bandwidth3[1],bandwidth4[1]))
L1=abs(screening(y1,x,h))
L2=abs(screening(y2,x,h))
L3=abs(screening(y3,x,h))
L4=abs(screening(y4,x,h))
L=L1^2+L2^2+L3^2+L4^2
x1=localMax(L,floor(h*n))/n
lambda=(4*mad(L1))^2+(4*mad(L2))^2+(4*mad(L3))^2+(4*mad(L4))^2
candidate=x1[which(L[localMax(L,floor(h*n))]>lambda)]
x2=c(0,(candidate[1:(length(candidate)-1)]+candidate[2:length(candidate)])/2,1.01)
pvalue=1:length(candidate)
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
}
common=candidate[which(p.adjust(pvalue, "BH")<0.1)]

###################################################show Table
rbind(c(num1, T1n), c(num2, T2n), c(num3, T3n), c(num4, T4n))

###################################################show Figure
par(mfrow=c(2,2))
###sequence1
jump=0 
LL=screening(y1,x,h)
for(i in 1:length(estimate1)){
  jump<-jump+LL[(estimate1[i]*n)]*(x>estimate1[i])
}
h3=dpill(x,(y1-jump))
fitted=npreg((y1-jump)~x,bws=h3,ckertype="epanechnikov")$mean+jump
plot(x*n, y1, xlim=c(1801,2300), ylim=c(-0.5,0.4), xlab="locations", ylab="Log 2 ratio", main="X1333-4", pch=20, col=8)
lines(x*n ,fitted, xlim=c(1801,2300), col=2, lwd=2) 
abline(v=common*n, xlim=c(1801,2300), lty=2)
###sequence2
jump=0 
LL=screening(y2,x,h)
for(i in 1:length(estimate2)){
  jump<-jump+LL[(estimate2[i]*n)]*(x>estimate2[i])
}
h3=dpill(x,(y2-jump))
fitted=npreg((y2-jump)~x,bws=h3,ckertype="epanechnikov")$mean+jump
plot(x*n, y2, xlim=c(1801,2300), ylim=c(-1,1), xlab="locations", ylab="Log 2 ratio", main="X1533-1", pch=20, col=8)
lines(x*n ,fitted, xlim=c(1801,2300), col=2, lwd=2) 
abline(v=common*n, xlim=c(1801,2300), lty=2)
###sequence3
jump=0 
LL=screening(y3,x,h)
for(i in 1:length(estimate3)){
  jump<-jump+LL[(estimate3[i]*n)]*(x>estimate3[i])
}
h3=dpill(x,(y3-jump))
fitted=npreg((y3-jump)~x,bws=h3,ckertype="epanechnikov")$mean+jump
plot(x*n, y3, xlim=c(1801,2300), ylim=c(-1,1), xlab="locations", ylab="Log 2 ratio", main="X1533-10", pch=20, col=8)
lines(x*n ,fitted, xlim=c(1801,2300), col=2, lwd=2) 
abline(v=common*n, xlim=c(1801,2300), lty=2)
###sequence4
jump=0 
LL=screening(y4,x,h)
for(i in 1:length(estimate4)){
  jump<-jump+LL[(estimate4[i]*n)]*(x>estimate4[i])
}
h3=dpill(x,(y4-jump))
fitted=npreg((y4-jump)~x,bws=h3,ckertype="epanechnikov")$mean+jump
plot(x*n, y4, xlim=c(1801,2300), ylim=c(-1,1), xlab="locations", ylab="Log 2 ratio", main="X1533-13", pch=20, col=8)
lines(x*n ,fitted, xlim=c(1801,2300), col=2, lwd=2) 
abline(v=common*n, xlim=c(1801,2300), lty=2)