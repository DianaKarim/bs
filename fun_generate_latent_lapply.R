

# New improved version from JM
generate.latent.lapply <- function(S, R, X, beta, Z, Sigma121, sigma2_j){
  
  # initialization of dimensions
  M <- length(R) # number of events
  N <- dim(as.matrix(X[[1]]))[1] # number of actors
  P <- dim(as.matrix(X[[1]]))[2] # number of covariates
  
  Z_out <- t(matrix(unlist(lapply(1:M,function(i){
    X[[i]] <- as.matrix(X[[i]]) # to avoid errors on P=1 case
    mu <- X[[i]]%*%beta
    Z_i <- Z[i,]
    
    for(j in (1:N)[-c(R[i],S[i],1)]){ #note that the first Z is fixed to 0
      
      mu_j <- c(mu[j] + Sigma121%*%(Z_i[-c(S[i],j)] - mu[-c(S[i],j)]))
      Z_i[j] <- rtruncnorm(1,mean=mu_j,sd=sqrt(sigma2_j),a=-Inf,b=Z_i[R[i]])
    }
    
    mu_Ri <- c(mu[R[i]] + Sigma121%*%(Z_i[-c(S[i],R[i])] - mu[-c(S[i],R[i])]))
    
    Z_i[R[i]] <- rtruncnorm(1,mean=mu_j,sd=sqrt(sigma2_j),a=max(Z_i[-c(S[i],R[i])]),b=Inf)*as.integer(R[i]!=1)
    return(Z_i)
  }
  
  )),nrow=N))
  
  return(Z_out)
  
}




