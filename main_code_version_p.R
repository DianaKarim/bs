library(truncnorm)
library(quantreg)
library(statmod)
library(ggplot2)
library(reshape2)
library(smoothmest)
library(parallel)
library(mvtnorm)



### Scenario 1 ----
# Data generation 
n <- 200
p <- 10
set.seed(1515)
x <- matrix(rnorm(n*p), nrow = n,ncol = p)

beta <- c(1,1,1,0,0,0,0.5,0.5,0.5,0.2)
sigma2 <- 1

y <- vector(length = n)

for (i in 1:n){
  z <- rnorm(1, mean = x[i,]%*%beta, sd = sigma2)  
  if (z>0) 
    y[i] <- 1 else
    y[i] <-  0
}

# Estimation 
n_train <- 50
y_train <- y[c(1:n_train)]
x_train <- x[c(1:n_train), ]

imp_p <- imp.sample.p(y = y_train, x = x_train, Nsample = 10000, burn_in = 2000, n = n_train, p = 10)
ridge_p <-  ridge.sample.p(y_train, x_train, Nsample = 10000, burn_in = 2000, n = n_train, p = 10)
lasso_p <-  lasso.sample.p(y_train, x_train, Nsample = 10000, burn_in = 2000, n = n_train, p = 10)
hs_p <- hs.sample.p(y_train, x_train, Nsample = 10000, burn_in = 2000, a2 = 0.7, a4 = 0.7, n = n_train, p = 10)

# "Goodness of fit" - checking the predictions 
y_test <- y[c((n_train+1):n)]
x_test <- x[c((n_train+1):n), ]

y_imp <-vector(length = length(y_test))
for (i in 1:(n/2)){
  z <- rnorm(1, mean = x_test[i,]%*%imp_p, sd = sigma2)  
  if (z>0) 
    y_imp[i] <- 1 else
      y_imp[i] <-  0
}
sum(y_imp == y_test)


y_ridge <-vector(length = length(y_test))
for (i in 1:(n/2)){
  z <- rnorm(1, mean = x_test[i,]%*%ridge_p, sd = sigma2)  
  if (z>0) 
    y_ridge[i] <- 1 else
      y_ridge[i] <-  0
}
sum(y_ridge == y_test)


y_lasso <-vector(length = length(y_test))
for (i in 1:(n/2)){
  z <- rnorm(1, mean = x_test[i,]%*%lasso_p, sd = sigma2)  
  if (z>0) 
    y_lasso[i] <- 1 else
      y_lasso[i] <-  0
}
sum(y_lasso == y_test)


y_hs <-vector(length = length(y_test))
for (i in 1:(n/2)){
  z <- rnorm(1, mean = x_test[i,]%*%hs_p, sd = sigma2)  
  if (z>0) 
    y_hs[i] <- 1 else
      y_hs[i] <-  0
}
sum(y_hs == y_test)


### Scenario 2 ----
# Data generation 
n <- 200
p <- 10
set.seed(1515)
x <- matrix(rnorm(n*p), nrow = n,ncol = p)

beta <- c(1,1,1,0,0,0,0.5,0.5,0.5,0.2)
sigma2 <- 1

y <- vector(length = n)

for (i in 1:n){
  z <- rnorm(1, mean = x[i,]%*%beta, sd = sigma2)  
  if (z>0) 
    y[i] <- 1 else
      y[i] <-  0
}



