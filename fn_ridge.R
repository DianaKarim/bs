
library(truncnorm)
library(mvtnorm)


ridge <-  function(S, R, X, Sigma, Nsample = 1000, burn_in = 200, r = 2, delta = 1/2){
  
  # Define dimensions
  M <- length(R) # number of observations
  N <- dim(as.matrix(X[[1]]))[1] # number of actors
  P <- dim(as.matrix(X[[1]]))[2]  # number of predictors
  
  # Define the vectors to store the results
  beta_ridge <- matrix(nrow = Nsample+1, ncol = P)
  beta_ridge[1, ] <- rep(1, P)
  
  lambda <- rep(0,Nsample)
  lambda[1] <- rgamma(1, r, delta)
  Z_saved <- list()
  
  #variance -???
  
  # Pre-calculated values 
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
  
  
  # Sampler ----
  for (t in 1:Nsample){
    # 1. Generate Z
    Z <- generate.latent.lapply(S, R, X, beta=beta_ridge[t,], Z=Z_old, Sigma121, sigma2_j)
    
    # 2. generate beta
    XSZ <- lapply(1:M,function(i){
      X[[i]] <- as.matrix(X[[i]])
      return(t(X[[i]][-S[i],])%*%inv_sigma[-S[i],-S[i]]%*%Z[i,-S[i]]) 
    })
    XSZ_sum <- Reduce("+", XSZ)
    
    var_beta <- solve(XSX_sum + diag(lambda[t],nrow = P))
    mu_beta <- var_beta %*% XSZ_sum
    
    beta_ridge[t+1, ] <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))
    
    
    # 4.  Generate lambda
  
    lambda[t+1] <- rgamma(1, shape = r+P/2, rate = delta + sum(beta_ridge[t+1,]^2)/2)
    
    # 4. Save the Z's for the next iteration:
    Z_saved[[t]] <- Z[1, ]
    Z_old <- Z  
  }
  
  list(Beta = beta_ridge, Lambda = lambda, Z = Z_saved)
  
}


