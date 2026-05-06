library(msm)
library(sn)
library(mixtools)
library(essHist)

a=1.35
b1=1.2
b2=1.2
b3=1.4
par(mfrow = c(2,3),mar=c(2.5,2.83,2,0),oma=c(1,1,1,1),mgp=c(1.8,0.6,0)) 

x5 <- seq(0.4, 1, length.out = 10000)
y5=dbeta(x5,shape1=20,shape2=2)
plot(x5, y5, type = "l", lwd = 2, col = "red",cex.lab=a, cex.main=b1,cex.axis =b2,cex=b3,main = "Beta(20,2)",xlab = "", ylab = "f(x)",ylim = c(0, 9))

x2 <- seq(0, 15, length.out = 10000)
y2=dweibull(x2,shape=1,scale=3)
plot(x2, y2, type = "l", lwd = 2, col = "red",cex.lab=a,cex.main=b1,cex.axis =b2,cex=b3,main = "Weibull(1,3)",xlab = "", ylab = "",ylim = c(0, 0.4))

x3 <- seq(0, 5, length.out = 10000)
y3=dtnorm(x3,0,1,lower=0, upper=Inf)
plot(x3, y3, type = "l", lwd = 2, col = "red",cex.lab=a,cex.main=b1,cex.axis =b2,cex=b3,main = "Truncated normal(0,1,0)",xlab = "", ylab = "",ylim = c(0, 0.9))


x1 <- seq(0, 5, length.out = 10000)
y1 <- dlnorm(x1, meanlog= 0, sdlog = 1)
plot(x1, y1, type = "l", lwd = 2, col = "red",cex.lab=a,cex.main=b1,cex.axis =b2,cex=b3,main = "Log-normal(0,1)",xlab = "x", ylab = "f(x)",ylim = c(0, 0.7))


x4 <- seq(-4, 2, length.out = 1000)
y4=dsn(x4,0,1,-10)
plot(x4, y4, type = "l", lwd = 2, col = "red",cex.lab=a,cex.main=b1,cex.axis =b2,cex=b3,main = "Skew-normal(0,1,-10)",xlab = "x", ylab = "",ylim = c(0, 0.8))



x6 <- seq(-4, 5, length.out = 10000)
y61=dnorm(x6,-2,1)
y62=dnorm(x6,2,1)
y6=0.1*y61+0.9*y62
plot(x6, y6, type = "l", lwd = 2, col = "red",cex.lab=a,cex.main=b1,cex.axis =b2,cex=b3,main = "0.1N(-2,1)+0.9N(2,1)",xlab = "x", ylab = "",ylim = c(0, 0.4))

dev.new()

a=1.35
b1=1.2
b2=1.2
b3=1.4
par(mfrow = c(2,3),mar=c(2.5,2.83,2,0),oma=c(1,1,1,1),mgp=c(1.8,0.6,0)) 

x1 <- seq(0, 1, length.out = 10000)
y1=dbeta(x1,shape1=3,shape2=2)
plot(x1, y1, type = "l", lwd = 2, col = "red",cex.lab=a, cex.main=b1,cex.axis =b2,cex=b3,main = "Beta(3,2)",xlab = "", ylab = "f(x)",ylim = c(0, 2))

x2 <- seq(0, 10, length.out = 10000)
y2=dweibull(x2,shape=2.5,scale=3)
plot(x2, y2, type = "l", lwd = 2, col = "red",cex.lab=a,cex.main=b1,cex.axis =b2,cex=b3,main = "Weibull(2.5,3)",xlab = "", ylab = "",ylim = c(0, 0.4))

x3 <- seq(-1.5, 5, length.out = 10000)
y3=dtnorm(x3,0,1,lower=-1.5, upper=Inf)
plot(x3, y3, type = "l", lwd = 2, col = "red",cex.lab=a,cex.main=b1,cex.axis =b2,cex=b3,main = "Truncated normal(0,1,-1.5)",xlab = "", ylab = "",ylim = c(0, 0.5))


x4 <- seq(0, 2.2, length.out = 10000)
y4 <- dlnorm(x4, meanlog= 0, sdlog = 0.1)
plot(x4, y4, type = "l", lwd = 2, col = "red",cex.lab=a,cex.main=b1,cex.axis =b2,cex=b3,main = "Log-normal(0,0.1)",xlab = "x", ylab = "f(x)",ylim = c(0, 5))


x5 <- seq(-4, 2, length.out = 1000)
y5=dsn(x5,0,1,-1.8)
plot(x5, y5, type = "l", lwd = 2, col = "red",cex.lab=a,cex.main=b1,cex.axis =b2,cex=b3,main = "Skew-normal(0,1,-1.8)",xlab = "x", ylab = "",ylim = c(0, 0.7))



x6 <- seq(-4, 5, length.out = 10000)
#y61=dnorm(x6,-2,1)
#y62=dnorm(x6,2,1)
#y6=0.5*y61+0.5*y62

mean_x=c(-2,2)
sd_x=c(1,1)
prob_x=c(0.3,0.7)
y6=dmixnorm(x6,mean=mean_x,sd=sd_x,prob=prob_x)
plot(x6, y6, type = "l", lwd = 2, col = "red",cex.lab=a,cex.main=b1,cex.axis =b2,cex=b3,main = "0.3N(-2,1)+0.7N(2,1)",xlab = "x", ylab = "",ylim = c(0, 0.35))

############
n=10000
x <- c(rnorm(ceiling(n*0.3),-2,1 ), rnorm(n-ceiling(n*0.3),2,1))
#x=x[which(x>=-4)]
#x=x[which(x<=5)]
dx=density(x,bw=0.3)
plot(dx,type = "l", lwd = 2, col = "red",cex.lab=a,cex.main=b1,cex.axis =b2,cex=b3,main = "0.3N(-2,1)+0.7N(2,1)",xlab = "x", ylab = "",ylim = c(0, 0.4),xlim=c(-4,5))
abline(h = 0, col = "white", lwd = 2,xlim=c(-4,5)) 
