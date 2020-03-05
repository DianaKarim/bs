
library(mvtnorm)
library(truncnorm)
library(quantreg)
library(statmod)
library(extraDistr)
#library(progress)
#library(coda)


lasso.new <- function (S, R, X, Sigma, Nsample = 500, burn_in = 200, r = 5, a = 0.5, b = 20){
  
  # Define dimensions
  M <- length(R) # number of observations
  N <- dim(as.matrix(X[[1]]))[1] # number of actors
  P <- dim(as.matrix(X[[1]]))[2]  # number of predictors
  
  # Define the vectors to store the results
  beta_lasso <- matrix(nrow = Nsample+1, ncol = P)
  beta_lasso[1, ] <- rep(1, P)
  
  delta <- rep(0, Nsample+1)
  delta[1] <- rinvgamma(1, alpha = a, beta = b)
  
  
  lambda2 <- rep(0,Nsample+1)
  lambda2[1] <- rgamma(1, shape = r, scale = delta[1])
  
  Z_saved <- list()
  
  tau <- rep(100, P) # prior variance of beta
  D_inv <- diag(tau^(-1), nrow = P) ##Diana: matrix D_tau^{-1}

  # Define the pre-iteration variables 
  
  inv_sigma <- solve(Sigma)
  Sigma11 <- Sigma[-1,-1] #note that this is the matrix we have to work with because the sender does not have a Z
  Sigma121 <- Sigma11[1,-1]%*%solve(Sigma11[-1,-1]) #needed in every step
  sigma2_j <- c(Sigma11[1,1] - t(Sigma11[1,-1])%*%solve(Sigma11[-1,-1])%*%Sigma11[-1,1]) #needed in every step
  
  XSX <- lapply(1:M,function(i){ #S added in name to better capture the meaning of the object
    X[[i]] <- as.matrix(X[[i]])
    t(X[[i]][-S[i],])%*%inv_sigma[-S[i],-S[i]]%*%X[[i]][-S[i],]	
  })
  XSX_sum <- Reduce("+", XSX) ## sum_i^M(X_i^T Sigma^(-1) X_i)
  
  
  # Initial values of Z's
  
  Z_old <- matrix(0,M,N)
  for (i in 1:M){ #JM: would be better to use lapply in case Z is large
    Z_old[i,S[i]] <- NA
    Z_old[i,-c(R[i],S[i])] <- rtruncnorm(N-2,a=rep(-Inf,N-2),b=rep(0,N-2),mean=0,sd=1)
  }
  
  
  # Sampler iterations
  
  for (t in 1:Nsample){
    # 1. Generate Z
    Z <- generate.latent.lapply(S, R, X, beta=beta_lasso[t,], Z=Z_old, Sigma121, sigma2_j)
    
    # 2. generate beta
    XSZ <- lapply(1:M,function(i){
      X[[i]] <- as.matrix(X[[i]])
      return(t(X[[i]][-S[i],])%*%inv_sigma[-S[i],-S[i]]%*%Z[i,-S[i]]) 
    })
    XSZ_sum <- Reduce("+", XSZ)
    
    var_beta <- solve(XSX_sum + D_inv)
    mu_beta <- var_beta %*% XSZ_sum
    
    beta_lasso[t+1, ] <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))
    
    # 3. Generate  1/tau^2 
    mu_tau <- sqrt((lambda2[t])/beta_lasso[t+1, ]^2)
    tau_sample <- unlist(lapply(1:P,function(p){ 
      rinvgauss(1, mean = mu_tau[p], shape = lambda2[t])
    }))
    
    D_inv <- diag(tau_sample, nrow = P)  
    
    # 4.  Generate lambda^2
    lambda2[t+1] <- rgamma(1, shape = r+P, scale = 1/delta[t] + 1/2*sum(1/tau_sample))
    
    # 5.  Generate delta
    
    delta[t+1] <- rinvgamma(1, alpha  = r + a, beta = lambda2[t+1]+b)
    
    # 4. Save the Z's for the next iteration:
    Z_saved[[t]] <- Z[1, ]
    Z_old <- Z  
  }
  
  list(Beta = beta_lasso, Lambda2 = lambda2, Z = Z_saved, Delta = delta) # return the list of the objects 
  #[c((burn_in+1):Nsample+1)]
}

#################################################################################################################
N <- 10 # number of actors
P <- 1 # enumber of covariates 
D <- 2 # number of messages for every potential receiver in most extreme case (maybe use more (less) depending on result)
M <- (N-1) * D #number of events 

S <- rep(1,M) #actor 1 is always sender (say)
X <- lapply(1:M,function(i){ #can be constructed in different ways
  X1 <- rep(0,N)
  X1[ceiling(i / D) + 1] <- 1
  Xout <- X1
  Xout[1] <- NA
  Xout
})
Sigma <- create.sigma(2.5, 0.5, d = N)

R1 <- rep(2,M) #one extreme where actor 2 is always receiver
R2 <- rep(2:N,each=D) #other extreme where receiver is exactly as expected by predictor variable X
#actual receivers go from extreme R1 to other extreme R2

step <- 5 #step can be D (for R1) ... N (for R2)
R <- c(R2[1:step],rep(2,M-step))
#####################################################################################################################

Nsample_test <- 1000
burn_in_test <- 200
r_test <- 5
a_test <- 0.5
b_test <- 20
res <- lasso.new(S,R,X,Nsample = Nsample_test, burn_in = burn_in_test,
                 Sigma, r = r_test, a = a_test, b = b_test)
# Lamdba with F 
xx <- seq(0, max(res$Lambda2), length.out = length(res$Lambda2)) 
yy <- df(xx, df1 = 2*r_test, df2 = 2*a_test, ncp = b_test) #prior 
dns_lambda <- density(res$Lambda2) 
hist(res$Lambda2, freq = F, ylim = c(0, max(dns_lambda$y, yy)),
     main = "", xlab = "", ylab = "")
title(main = "Lambda2", outer = F, 
      sub = paste0("r = ",  r_test))
lines(xx,yy, type = "l", col = "blue", ylim = c(min(yy), max(yy)))
lines(density(res$Lambda2), col = "red", lwd = 1.2)
legend("topright", 
       legend = c("Prior", "Posterior"), 
       col = c("blue", "red"), 
       text.col = "black",
       lwd = c(1.1,1.1),
       horiz = F)

# Delta

xx <- seq(0, max(res$Delta), length.out = Nsample_test)
yy <- dinvgamma(xx, alpha = a_test, b = b_test)
dns_delta <- density(res$Delta) 
hist(res$Delta, freq = F, ylim = c(0, max(dns_delta$y, yy)),
     main = "", xlab = "", ylab = "")
title(main = "Delta", outer = F, 
      sub = paste0("a = ",  a_test, ", b = ", b_test))
lines(xx,yy, type = "l", col = "blue" )
lines(density(res$Delta), col = "red", lwd = 1.2)
legend("topright", 
       legend = c("Prior", "Posterior"), 
       col = c("blue","red"), 
       #text.col = "black",
       lwd = c(1,1),
       inset = c(0.001, 0.001))


####################################################################

library(bayesplot)
mcmc_trace(data.frame(beta_lasso))
mcmc_trace(data.frame(res$Lambda2))

plot(res$Beta, type = "l")
hist(res$Delta)
hist(res$Lambda2)


