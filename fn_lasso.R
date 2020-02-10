
library(mvtnorm)
library(truncnorm)
library(quantreg)
library(statmod)
#library(progress)
#library(coda)


lasso <- function (S, R, X, Sigma, Nsample = 500, burn_in = 200, r = 2, delta = 1/2){

  # Define dimensions
  M <- length(R) # number of observations
  N <- dim(as.matrix(X[[1]]))[1] # number of actors
  P <- dim(as.matrix(X[[1]]))[2]  # number of predictors
  
  # Define the vectors to store the results
  beta_lasso <- matrix(nrow = Nsample+1, ncol = P)
  beta_lasso[1, ] <- rep(1, P)
  
  lambda2 <- rep(0,Nsample+1)
  lambda2[1] <- rgamma(1, r, delta)
  
  Z_saved <- list()
  
  tau <- rep(100, P) # prior variance of beta
  D_inv <- diag(tau^(-1), nrow = P) ##Diana: matrix D_tau^{-1}
  #Z - ???
  
  
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
    lambda2[t+1] <- rgamma(1, shape = r+P, rate = delta + 1/2*sum(1/tau_sample))
    
    # 4. Save the Z's for the next iteration:
    Z_saved[[t]] <- Z[1, ]
    Z_old <- Z  
    }
  
  list(Beta = beta_lasso, Lambda2 = lambda2, Z = Z_saved) # return the list of the objects 
  
}

library(bayesplot)
mcmc_trace(data.frame(beta_lasso))
mcmc_trace(data.frame(lambda2))

plot(beta_lasso[1,])

hist(lambda2)


res <- lasso(S,R,X, Sigma, Nsample = 100, burn_in = 20)
res$Beta
res$Lambda2



