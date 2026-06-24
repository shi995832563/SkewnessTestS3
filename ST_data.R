set.seed(1000)
library(estmeansd)


x0=read.table("ST_data.txt",header=F)
x1=c(t(x0))
w=which(x1<0)
x=x1[-w]
n=length(x)

x2=x
c(mean(x2),sd(x2))
############################
x=log(x2)

q1=min(x)
q2=quantile(x,0.25)
q3=quantile(x,0.5)
q4=quantile(x,0.75)
q5=max(x)

S3=max(2.65*log(0.6*n)/sqrt(n)*abs(q1+q5-2*q3)/(q5-q1),abs(q4+q2-2*q3)/(q4-q2))

ZD=(q5-q4)/(q2-q1)
T3=abs(log(ZD))

kn=1.65
wa=1-0.7*(log10(n)-2)
wb=1

T3B=kn*((q1+q5-2*q3)*wa+(q4+q2-2*q3)*(1-wa))/(q5-q1)

sc3=3/sqrt(n)-40/n^3
cl=3.23/log(n^1.25-3.1)+8.45/(n+1.5)
cb=1/((15*0.05+0.6)^0.4/100*n+(26*0.05+0.8)^0.4)

c(T3,S3,T3B)
c(cl,sc3,cb)

######################################
w1=0.5*2.2/(2.2+n^0.75)
w2=0.5*(0.7-0.72/n^0.55)

r1=1/(1+0.07*n^0.6)
r2=1-r1

muhat0=w1*(q1+q5)+w2*(q2+q4)+(1-2*w1-2*w2)*q3

xi1=2*qnorm((n-0.375)/(n+0.25),0,1)
eta1=2*qnorm((0.75*n-0.125)/(n+0.25),0,1)
shat1=r1*(q5-q1)/xi1+r2*(q4-q2)/eta1
shat0=shat1/sqrt(1+0.28/(log(n))^2)


muhat1=exp(muhat0+shat0^2/2)
shat1=sqrt((exp(shat0^2)-1))*muhat1
c(muhat1,shat1)
c((muhat1-mean(x2))/mean(x2),(shat1-sd(x2))/sd(x2))
############################
x=x2
q1=min(x)
q2=quantile(x,0.25)
q3=quantile(x,0.5)
q4=quantile(x,0.75)
q5=max(x)

S3=max(2.65*log(0.6*n)/sqrt(n)*abs(q1+q5-2*q3)/(q5-q1),abs(q4+q2-2*q3)/(q4-q2))

ZD=(q5-q4)/(q2-q1)
T3=abs(log(ZD))
T3B=kn*((q1+q5-2*q3)*wa+(q4+q2-2*q3)*(1-wa))/(q5-q1)


c(T3,S3,T3B)

###########################
Blandm=((n+3)*q1+2*(n-1)*(q2+q3+q4)+(n+3)*q5)/8/n
Bs=(2*(n+3)*(q2^2+q3^2+q4^2)+2*(n-5)*(q1*q2+q2*q3+q3*q4+q4*q5)+(n+11)*(q1^2+q5^2))/16
Blands=sqrt((Bs-n*Blandm^2)/(n-1))
c(Blandm,Blands)
c((Blandm-mean(x2))/mean(x2),(Blands-sd(x2))/sd(x2))
################################
Bc=bc.mean.sd(min.val=q1,q1.val = q2, med.val = q3, q3.val =q4, max.val=q5,n = n)
muhat2=Bc$est.mean
shat2=Bc$est.sd
c(muhat2,shat2)
c((muhat2-mean(x2))/mean(x2),(shat2-sd(x2))/sd(x2))

QE=qe.mean.sd(min.val=q1,q1.val = q2, med.val = q3, q3.val =q4, max.val=q5,n = n)
muhat3=QE$est.mean
shat3=QE$est.sd
c(muhat3,shat3)
c((muhat3-mean(x2))/mean(x2),(shat3-sd(x2))/sd(x2))

MLN=mln.mean.sd(min.val=q1,q1.val = q2, med.val = q3, q3.val =q4, max.val=q5,n = n)
muhat4=MLN$est.mean
shat4=MLN$est.sd
c(muhat4,shat4)
c((muhat4-mean(x2))/mean(x2),(shat4-sd(x2))/sd(x2))
#########################
#LS
muhat0=w1*(q1+q5)+w2*(q2+q4)+(1-2*w1-2*w2)*q3

xi1=2*qnorm((n-0.375)/(n+0.25),0,1)
eta1=2*qnorm((0.75*n-0.125)/(n+0.25),0,1)
shat0=r1*(q5-q1)/xi1+r2*(q4-q2)/eta1
c(muhat0,shat0)
c((muhat0-mean(x2))/mean(x2),(shat0-sd(x2))/sd(x2))
##################################









