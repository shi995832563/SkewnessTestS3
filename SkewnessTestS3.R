SkewnessTestS3 <- function(x,n,repe=10000)
{
  DNAME <- deparse(substitute(x))
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(n))
  stopifnot(x[5]>x[4] | x[2]>x[1])
  data <- array(rnorm(n*repe),c(n,repe))
  five_test <- apply(data, 2, fivenum)
  t3 <- log((five_test[5,]-five_test[4,])/(five_test[2,]-five_test[1,]))
  t3_test <- log((x[5]-x[4])/(x[2]-x[1]))
  pval <- sum(abs(t3_test)<abs(t3))/repe
  
  RVAL <- list(statistic = c(ST3 = t3_test), p.value = pval, 
               method = "Skewness test S3", data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
  
}

###Example###
n <- 100
data <- rnorm(n)
five_num <-fivenum(data)
SkewnessTestS3(five_num,n)