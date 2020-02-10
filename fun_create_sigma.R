
# Generate covariance matrix

create.sigma <- function(sigma2, rho, d){ 
  e <- c(rep(1,d)) 
  Sigma<-sigma2*((1-rho)*diag(d)+rho*e%*%t(e))
  return(Sigma)
} # d should be equal to J


range01 <- function(x){(x-min(x))/(max(x)-min(x))}
