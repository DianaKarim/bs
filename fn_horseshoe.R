
library(truncnorm)
library(extraDistr)
library(mvtnorm)

horseshoe <- function(S, R, X, Sigma, Nsample = 500, burn_in = 200, a2 = 0.7, a4 = 0.7){
  
  # Define dimensions
  M <- length(R) # number of observations
  N <- dim(as.matrix(X[[1]]))[1] # number of actors
  P <- dim(as.matrix(X[[1]]))[2]  # number of predictors
  
  # Define the vectors to store the results
  beta_hs <- matrix(nrow = Nsample+1, ncol = P)
  beta_hs[1, ] <- rep(1, P)
  lambda2 <- rep(0, Nsample +1)
  
  Z_saved <- list()
  
  
  tau2 <- rep(100, P) # prior variance of beta
  D_inv <- diag(tau2^(-1), nrow = P) ##Diana: matrix D_tau^{-1}
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
  
  # Hyperparameters 
  ## stable
  a1 <- 1/2
  a3 <- 1/2
  b <- 1
  ## changing over iterations
  gma <- 1
  psi2 <- 1
  
  # Srinkage parameter

  lambda2[1] <- rgamma(1, a3, psi2)

  
  # Sampler iterations
  
  for (t in 1:Nsample){
    # 1. Generate Z
    Z <- generate.latent.lapply(S, R, X, beta=beta_hs[t,], Z=Z_old, Sigma121, sigma2_j)
    
    # 2. generate beta
    XSZ <- lapply(1:M,function(i){
      X[[i]] <- as.matrix(X[[i]])
      return(t(X[[i]][-S[i],])%*%inv_sigma[-S[i],-S[i]]%*%Z[i,-S[i]]) 
    })
    XSZ_sum <- Reduce("+", XSZ)
    
    var_beta <- solve(XSX_sum + D_inv)
    mu_beta <- var_beta %*% XSZ_sum
    
    beta_hs[t+1, ] <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))
    
    # 3. Generate  tau^2 
    tau2 <- unlist(lapply(1:P,function(p){ 
      rinvgamma(1, 1/2 + a1, (beta_hs[t+1,p]^2+gma)/2 )
    }))
    
    D_inv <- diag(1/tau2, nrow = P)  

    
    # 4.  Generate gma
    gma <- rgamma(1, shape = a1 + a2, rate =1/2*(1/lambda2[t] +sum(1/tau2)))
    
    # 5. Generate lambda2
    
    lambda2[t+1] <- rinvgamma(1, a2+a3, psi2 + gma/2)
    
    # 6. Generate psi
    psi2 <- rgamma(1, shape = a3 + a4, rate = 1/b +1/lambda2[t+1])
    
    # 4. Save the Z's for the next iteration:
    Z_saved[[t]] <- Z[1, ]
    Z_old <- Z  
  }
  
  list(Beta = beta_hs, Lambda2 = lambda2, Z = Z_saved) # return the list of the objects 
  
}
