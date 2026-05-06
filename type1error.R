set.seed(2000)
rm(list=ls())
gc()

mytime <- Sys.time()
library(sn)
library(stats)
library(ggplot2)
library(actuar)



a=1:100
m=4*a+1
K=3000000

B1=matrix(0,nrow=length(a),ncol=3)

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

for (i in 1:length(a)){
n=m[i]

sc3=3/sqrt(n)-40/n^3
sc2=2.65/sqrt(n)-6/n^2
sc1=1/log(n+9)+2.5/(n+1)


c2_1=3.23/log(n^1.25-3.1)+8.45/(n^1+1.5)

#######
cb=1/((15*0.05+0.6)^0.4/100*n+(26*0.05+0.8)^0.4)



if (n<=100) {wa=1} else {wa=1-0.7*(log10(n)-2)}  
if (n<=200) {wb=0.1} else {wb=1}
if (n<=200) {kn=0.55} else {kn=1.65}

######################3
##########################
A1=matrix(0,nrow=K,ncol=3)

for (k in 1:K){

#####################
x=rnorm(n,0,1)


result=mytsfun(x)
S3=result$S3
ZD=result$ZD
TB3=result$TB3


if (abs(ZD)>c2_1){A1[k,1]=1}
if (S3>sc3){A1[k,2]=1}
if (abs(TB3)>cb){A1[k,3]=1}
##################

################
}
B1[i,]=colMeans(A1)
}

B1


##############################3

Sys.time()-mytime

############################3
############################3
library(latex2exp)

B0=B1
f=seq(1,100,2)
m1=m
B2=B1
plot(m1,B2[,1],type="o",lty=1,pch=20,cex=1.05,col="red",ylim=c(0,0.2),ylab="Type I error rate", xlab="n")
lines(m1,B2[,2],type="o",lty=1,pch=1,cex=0.85,col="blue")
lines(m1,B2[,3],type="o",lty=1,pch=17,cex=0.85,col="black")
abline(h=0.05)


legend(350,0.15,c(TeX(r'($T_{3}$)'),TeX(r'($T_{3}^{S}$)'),TeX(r'($T_{3}^{B}$)')),lty=c(1,1,1),pch=c(20,1,17),col=c("red","blue","black"))

#save(B1,m, file="typeI300w251212_ratio323.RData")
