set.seed(2000)
rm(list=ls())
gc()

mytime <- Sys.time()
library(sn)
library(stats)
library(ggplot2)
library(LaplacesDemon)
library(actuar)
library(msm)  #trunc_norm



a=1:100
m=4*a+1
K=1000000

B1=B2=B3=B4=B5=B6=matrix(0,nrow=length(a),ncol=3)

#######################3
mytsfun<-function(x){
q1=min(x)
q2=quantile(x,0.25)
q3=quantile(x,0.5)
q4=quantile(x,0.75)
q5=max(x)


S3=max(2.65*log(0.6*n)/sqrt(n)*abs(q1+q5-2*q3)/(q5-q1),abs(q4+q2-2*q3)/(q4-q2))

ZD=(q5-q4)/(q2-q1)
ZD=log(ZD)

TB3=kn*((q1-2*q3+q5)*wa+(q2-2*q3+q4)*(1-wa))/((q5-q1)*wb+(q4-q2)*(1-wb))

result=list(S3=S3,ZD=ZD, TB3=TB3)
return(result) 
}


######################3

for (i in a){
n=m[i]


sc3=3/sqrt(n)-40/n^3

######################3 
c2_1=3.23/log(n^1.25-3.1)+8.45/(n^1+1.5)
######
cb=1/((15*0.05+0.6)^0.4/100*n+(26*0.05+0.8)^0.4)



if (n<=100) {wa=1} else {wa=1-0.7*(log10(n)-2)}  
if (n<=200) {wb=0.1} else {wb=1}
if (n<=200) {kn=0.55} else {kn=1.65}



##########################
A1=A2=A3=A4=A5=A6=matrix(0,nrow=K,ncol=3)

for (k in 1:K){

#####################

x=rlnorm(n, meanlog = 0, sdlog = 1)

result=mytsfun(x)

S3=result$S3
ZD=result$ZD
TB3=result$TB3


if (abs(ZD)>c2_1){A1[k,1]=1}
if (S3>sc3){A1[k,2]=1}
if (abs(TB3)>cb){A1[k,3]=1}
##################
x=rweibull(n, shape=1, scale = 3) 

result=mytsfun(x)

S3=result$S3
ZD=result$ZD
TB3=result$TB3


if (abs(ZD)>c2_1){A2[k,1]=1}
if (S3>sc3){A2[k,2]=1}
if (abs(TB3)>cb){A2[k,3]=1}

################
##################
x1 <- rnorm(n,-2,1 )
x2 <- rnorm(n,2,1)
p<- runif(n,0,1)
pa=0.1   
y1<- p <pa
x<- x1 *y1+x2*(1-y1 )



result=mytsfun(x)

S3=result$S3
ZD=result$ZD
TB3=result$TB3

if (abs(ZD)>c2_1){A3[k,1]=1}
if (S3>sc3){A3[k,2]=1}
if (abs(TB3)>cb){A3[k,3]=1}
################
##################
x=rtnorm(n,0,1,lower=0,upper=Inf)  #trunc_norm

result=mytsfun(x)

S3=result$S3
ZD=result$ZD
TB3=result$TB3

if (abs(ZD)>c2_1){A4[k,1]=1}
if (S3>sc3){A4[k,2]=1}
if (abs(TB3)>cb){A4[k,3]=1}

################
##################

x=rsn(n, 0, 1, -10)  #-10 

result=mytsfun(x)

S3=result$S3
ZD=result$ZD
TB3=result$TB3


if (abs(ZD)>c2_1){A5[k,1]=1}
if (S3>sc3){A5[k,2]=1}
if (abs(TB3)>cb){A5[k,3]=1}

#####################

x=rbeta(n,shape1=20,shape2=2) 

result=mytsfun(x)

S3=result$S3
ZD=result$ZD
TB3=result$TB3

if (abs(ZD)>c2_1){A6[k,1]=1}
if (S3>sc3){A6[k,2]=1}
if (abs(TB3)>cb){A6[k,3]=1}



################
}
B1[i,]=colMeans(A1)
B2[i,]=colMeans(A2)
B3[i,]=colMeans(A3)
B4[i,]=colMeans(A4)
B5[i,]=colMeans(A5)
B6[i,]=colMeans(A6)
}
B1
B2
B3
B4
B5
B6
##############################3

Sys.time()-mytime

############################3

#save(B1,B2, B3,B4,B5, m, file="power_ratioBala_h.RData")

