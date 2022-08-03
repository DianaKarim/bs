

#functions for model fitting
generate.indlatent.lapply.dyad.stack <- function(Events, Xstack, beta, Z_max){
  
  # initialization of dimensions
  M <- nrow(Events) # number of events
  N <- nrow(Xstack)/M # number of all possible dyads
  P <- length(beta) # number of covariates
  
  muStack <- matrix(c(Xstack%*%beta),byrow=TRUE,ncol=N)
  
  #draw z's for nonactive events
  Z_out <- t(matrix(unlist(lapply(1:M,function(i){
    
    nonactive <- (1:N)[-Events[i,3]]
    active <- Events[i,3]
    
    #mu <- na.omit(X[i,,])%*%beta
    Z_i <- rep(0,N)
    Z_i[nonactive] <- rtruncnorm(N-1,mean=muStack[i,nonactive],sd=1,a=-Inf,b=Z_max[i])
    Z_i[1] <- 0
    Z_i[active] <- rtruncnorm(1,mean=muStack[i,active],sd=1,a=max(Z_i[nonactive]),b=Inf)
    Z_i[1] <- 0
    return(Z_i)
  })),nrow=N))
  return(Z_out)
}

# Gibbs sampler under horseshoe prior with global / local part
horseshoe.dyad2 <- function(Events, X, Nsample = 500, a1 = .5, a2 = 0.5, a3 = .5, a4 = 0.5,
                            b1 = 1, b2 = 1, store = 10, burnin = 10){
  # D
  print("Define pre-iteration variables ...")
  # Define dimensions
  M <- nrow(Events) # number of events
  N <- dim(X)[2] # number of all possible dyads
  P <- dim(X)[3] # number of covariates without intercept
  effect_names <- dimnames(X)[[3]]
  
  # Define the vectors to store the results
  beta_STORE <- matrix(0,nrow = Nsample/store, ncol = P)
  lambda2_STORE <- psi2_STORE <- matrix(0,nrow = Nsample/store, ncol = P)
  tau2_STORE <- gamma2_STORE <- matrix(0,nrow = Nsample/store, ncol = 1)
  predcheck1_STORE <- predcheck2_STORE <- matrix(0,nrow = Nsample/store, ncol = M)
  
  # initial values
  tau2 <- gamma2 <- 1
  lambda2 <- psi2 <- rep(1, P)
  D <- diag(tau2*lambda2)
  D_inv <- diag((tau2*lambda2)^(-1))
  beta <- rep(0,P)
  
  #initial computation for posterior cov matrix beta
  # time0 <- Sys.time()
  # XSX <- lapply(1:M,function(i){ #S added in name to better capture the meaning of the object
  #   t(X[i,,])%*%X[i,,]	
  # })
  # XSX_sum <- Reduce("+", XSX)
  # Sys.time() - time0
  
  #a stacked covariates matrix seems faster to compute with than the array
  Xstack <- lapply(1:M,function(i){X[i,,]})
  Xstack <- do.call(rbind,Xstack)
  #time0 <- Sys.time()
  XSX_sum <- t(Xstack)%*%Xstack
  #Sys.time() - time0
  
  # Initial values of Z's which are NOT USED
  Z <- matrix(1,M,N)
  for (i in 1:M){ #JM: would be better to use lapply in case Z is large
    Z[i,-Events[i,3]]<- rtruncnorm(N-1,a=rep(-Inf,N-1),b=rep(0,N-1),mean=0,sd=1)
    Z[,1] <- 0  
  }
  
  print("Start burnin ... ")
  # Sampler iterations ----
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = burnin, clear = F, width = 80)
  for (t in 1:burnin){
    # 1. Sample Z
    #time0 <- Sys.time()
    Z <- generate.indlatent.lapply.dyad.stack(Events, Xstack, beta, Z_max=apply(Z,1,max))
    #Sys.time() - time0
    
    # 2. sample beta
    #time0 <- Sys.time()
    XSZ_sum <- c(t(c(t(Z))) %*% Xstack)
    #Sys.time() - time0
    
    var_beta <- solve(XSX_sum + D_inv)
    mu_beta <- var_beta %*% XSZ_sum
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))
    
    # 3. Sample tau2 
    tau2 <- rinvgamma(1, a1 + P/2, gamma2 + sum(beta^2/(2*lambda2)) )
    
    # 4. Sample lambda2
    lambda2 <- rinvgamma(P, a3 + .5, psi2 + beta**2/(2*tau2))
    D_inv <- diag(1/(tau2*lambda2))
    
    # 5.  Sample gamma2
    gamma2 <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/tau2)
    
    # 6. Sample psi2
    psi2 <- stats::rgamma(P, shape = a3 + a4, rate = 1/b2 +1/lambda2)
    
    pb$tick()
    Sys.sleep(1/burnin)
  }
  
  print("Start iterations ... ")
  # Sampler iterations ----
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = Nsample, clear = F, width = 80)
  storecount <- 0
  for (t in 1:Nsample){
    # 1. Sample Z
    Z <- generate.indlatent.lapply.dyad.stack(Events, Xstack, beta, Z_max=apply(Z,1,max))
    
    # 2. sample beta
    XSZ_sum <- c(t(c(t(Z))) %*% Xstack)
    var_beta <- solve(XSX_sum + D_inv)
    mu_beta <- var_beta %*% XSZ_sum
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))
    
    # 3. Sample tau2 
    tau2 <- rinvgamma(1, a1 + P/2, gamma2 + sum(beta^2/(2*lambda2)) )
    
    # 4. Sample lambda2
    lambda2 <- rinvgamma(P, a3 + .5, psi2 + beta**2/(2*tau2))
    D_inv <- diag(1/(tau2*lambda2))
    
    # 5.  Sample gamma2
    gamma2 <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/tau2)
    
    # 6. Sample psi2
    psi2 <- stats::rgamma(P, shape = a3 + a4, rate = 1/b2 +1/lambda2)
    
    # save each 'store' iterations 
    if (t%%store == 0){
      #update store counter
      storecount <- storecount + 1
      #store draws
      beta_STORE[storecount,] <- beta
      lambda2_STORE[storecount,] <- lambda2
      tau2_STORE[storecount,] <- tau2
      gamma2_STORE[storecount,] <- gamma2
      psi2_STORE[storecount,] <- psi2
      
      # prediction performance checks
      Z_predcheck1 <- unlist(lapply(1:M,function(i){
        Z_i <- c(X[i,,]%*%beta + rnorm(N))
        Z_i[1] <- 0
        #rank of actually observed event in this predicted data
        (N + 1 - rank(Z_i))[Events[i,3]]
      }))
      # below the random part of Z is omitted.
      Z_predcheck2 <- unlist(lapply(1:M,function(i){
        Z_i <- c(X[i,,]%*%beta)
        Z_i[1] <- 0
        #rank of actually observed event in this predicted data
        (N + 1 - rank(Z_i))[Events[i,3]]
      }))
      
      predcheck1_STORE[storecount,] <- Z_predcheck1
      predcheck2_STORE[storecount,] <- Z_predcheck2
      
      #a13_horseshoe_model18 <- beta_hs[c(1:t),]
      #colnames(a13_horseshoe_model18) <- effect_names
      #save(a13_horseshoe_model18, file = paste0("a13_hs_model18_chunk_", t, ".RData"))
      #print("Saved chunk .. ")
    }
    ###################
    
    pb$tick()
    Sys.sleep(1/Nsample)
  }
  
  colnames(beta_STORE) <- dimnames(X)[[3]]
  
  return(list(beta = beta_STORE, predcheck = list(predcheck1_STORE,predcheck2_STORE), 
              lambda2 = lambda2_STORE, gamma2 = gamma2_STORE, psi2 = psi2_STORE, 
              tau2 = tau2_STORE,Z = Z))
}

# Gibbs sampler under flat prior
flat.dyad2 <- function(Events, X, Nsample = 500, store = 10, burnin = 10){
  # D
  print("Define pre-iteration variables ...")
  # Define dimensions
  M <- nrow(Events) # number of events
  N <- dim(X)[2] # number of all possible dyads
  P <- dim(X)[3] # number of covariates without intercept
  effect_names <- dimnames(X)[[3]]
  
  # Define the vectors to store the results
  beta_STORE <- matrix(0,nrow = Nsample/store, ncol = P)
  predcheck1_STORE <- predcheck2_STORE <- matrix(0,nrow = Nsample/store, ncol = M)
  
  # initial values
  beta <- rep(0,P)
  
  #initial computation for posterior cov matrix beta
  # time0 <- Sys.time()
  # XSX <- lapply(1:M,function(i){ #S added in name to better capture the meaning of the object
  #   t(X[i,,])%*%X[i,,]	
  # })
  # XSX_sum <- Reduce("+", XSX)
  # Sys.time() - time0
  
  #a stacked covariates matrix seems faster to compute with than the array
  Xstack <- lapply(1:M,function(i){X[i,,]})
  Xstack <- do.call(rbind,Xstack)
  #time0 <- Sys.time()
  XSX_sum <- t(Xstack)%*%Xstack
  #Sys.time() - time0
  
  # Initial values of Z's which are NOT USED
  Z <- matrix(1,M,N)
  for (i in 1:M){ #JM: would be better to use lapply in case Z is large
    Z[i,-Events[i,3]]<- rtruncnorm(N-1,a=rep(-Inf,N-1),b=rep(0,N-1),mean=0,sd=1)
    Z[,1] <- 0  
  }
  
  print("Start burnin ... ")
  # Sampler iterations ----
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = burnin, clear = F, width = 80)
  for (t in 1:burnin){
    # 1. Sample Z
    #time0 <- Sys.time()
    Z <- generate.indlatent.lapply.dyad.stack(Events, Xstack, beta, Z_max=apply(Z,1,max))
    #Sys.time() - time0
    
    # 2. sample beta
    #time0 <- Sys.time()
    XSZ_sum <- c(t(c(t(Z))) %*% Xstack)
    #Sys.time() - time0
    
    var_beta <- solve(XSX_sum)
    mu_beta <- var_beta %*% XSZ_sum
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))
    
    pb$tick()
    Sys.sleep(1/burnin)
  }
  
  print("Start iterations ... ")
  # Sampler iterations ----
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = Nsample, clear = F, width = 80)
  storecount <- 0
  for (t in 1:Nsample){
    # 1. Sample Z
    Z <- generate.indlatent.lapply.dyad.stack(Events, Xstack, beta, Z_max=apply(Z,1,max))
    
    # 2. sample beta
    XSZ_sum <- c(t(c(t(Z))) %*% Xstack)
    var_beta <- solve(XSX_sum)
    mu_beta <- var_beta %*% XSZ_sum
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))
    
    # save each 'store' iterations 
    if (t%%store == 0){
      #update store counter
      storecount <- storecount + 1
      #store draws
      beta_STORE[storecount,] <- beta
      
      # prediction performance checks
      Z_predcheck1 <- unlist(lapply(1:M,function(i){
        Z_i <- c(X[i,,]%*%beta + rnorm(N))
        Z_i[1] <- 0
        #rank of actually observed event in this predicted data
        (N + 1 - rank(Z_i))[Events[i,3]]
      }))
      # below the random part of Z is omitted.
      Z_predcheck2 <- unlist(lapply(1:M,function(i){
        Z_i <- c(X[i,,]%*%beta)
        Z_i[1] <- 0
        #rank of actually observed event in this predicted data
        (N + 1 - rank(Z_i))[Events[i,3]]
      }))
      
      predcheck1_STORE[storecount,] <- Z_predcheck1
      predcheck2_STORE[storecount,] <- Z_predcheck2
      
      #a13_horseshoe_model18 <- beta_hs[c(1:t),]
      #colnames(a13_horseshoe_model18) <- effect_names
      #save(a13_horseshoe_model18, file = paste0("a13_hs_model18_chunk_", t, ".RData"))
      #print("Saved chunk .. ")
    }
    ###################
    
    pb$tick()
    Sys.sleep(1/Nsample)
  }
  
  colnames(beta_STORE) <- dimnames(X)[[3]]
  
  return(list(beta = beta_STORE, predcheck = list(predcheck1_STORE,predcheck2_STORE), 
              Z = Z))
}

# Gibbs sampler under ridge prior
ridge.dyad2 <- function(Events, X, Nsample = 500, a1 = .5, a2 = 0.5,
                        b1 = 1, store = 10, burnin = 10){
  # D
  print("Define pre-iteration variables ...")
  # Define dimensions
  M <- nrow(Events) # number of events
  N <- dim(X)[2] # number of all possible dyads
  P <- dim(X)[3] # number of covariates without intercept
  effect_names <- dimnames(X)[[3]]
  
  # Define the vectors to store the results
  beta_STORE <- matrix(0,nrow = Nsample/store, ncol = P)
  #lambda2_STORE <- psi2_STORE <- matrix(0,nrow = Nsample/store, ncol = P)
  tau2_STORE <- gamma2_STORE <- matrix(0,nrow = Nsample/store, ncol = 1)
  predcheck1_STORE <- predcheck2_STORE <- matrix(0,nrow = Nsample/store, ncol = M)
  
  # initial values
  tau2 <- gamma2 <- 1
  lambda2 <- psi2 <- rep(1, P)
  D <- diag(tau2*lambda2)
  D_inv <- diag((tau2*lambda2)^(-1))
  beta <- rep(0,P)
  
  #initial computation for posterior cov matrix beta
  # time0 <- Sys.time()
  # XSX <- lapply(1:M,function(i){ #S added in name to better capture the meaning of the object
  #   t(X[i,,])%*%X[i,,]	
  # })
  # XSX_sum <- Reduce("+", XSX)
  # Sys.time() - time0
  
  #a stacked covariates matrix seems faster to compute with than the array
  Xstack <- lapply(1:M,function(i){X[i,,]})
  Xstack <- do.call(rbind,Xstack)
  #time0 <- Sys.time()
  XSX_sum <- t(Xstack)%*%Xstack
  #Sys.time() - time0
  
  # Initial values of Z's which are NOT USED
  Z <- matrix(1,M,N)
  for (i in 1:M){ #JM: would be better to use lapply in case Z is large
    Z[i,-Events[i,3]]<- rtruncnorm(N-1,a=rep(-Inf,N-1),b=rep(0,N-1),mean=0,sd=1)
    Z[,1] <- 0  
  }
  
  print("Start burnin ... ")
  # Sampler iterations ----
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = burnin, clear = F, width = 80)
  for (t in 1:burnin){
    # 1. Sample Z
    #time0 <- Sys.time()
    Z <- generate.indlatent.lapply.dyad.stack(Events, Xstack, beta, Z_max=apply(Z,1,max))
    #Sys.time() - time0
    
    # 2. sample beta
    #time0 <- Sys.time()
    XSZ_sum <- c(t(c(t(Z))) %*% Xstack)
    #Sys.time() - time0
    
    var_beta <- solve(XSX_sum + D_inv)
    mu_beta <- var_beta %*% XSZ_sum
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))
    
    # 3. Sample tau2 
    tau2 <- rinvgamma(1, a1 + P/2, gamma2 + sum(beta^2/(2*lambda2)) )
    
    # 4. Sample lambda2
    #lambda2 <- rinvgamma(P, a3 + .5, psi2 + beta**2/(2*tau2))
    D_inv <- diag(1/(tau2*lambda2))
    
    # 5.  Sample gamma2
    gamma2 <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/tau2)
    
    # 6. Sample psi2
    #psi2 <- stats::rgamma(P, shape = a3 + a4, rate = 1/b2 +1/lambda2)
    
    pb$tick()
    Sys.sleep(1/burnin)
  }
  
  print("Start iterations ... ")
  # Sampler iterations ----
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = Nsample, clear = F, width = 80)
  storecount <- 0
  for (t in 1:Nsample){
    # 1. Sample Z
    Z <- generate.indlatent.lapply.dyad.stack(Events, Xstack, beta, Z_max=apply(Z,1,max))
    
    # 2. sample beta
    XSZ_sum <- c(t(c(t(Z))) %*% Xstack)
    var_beta <- solve(XSX_sum + D_inv)
    mu_beta <- var_beta %*% XSZ_sum
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))
    
    # 3. Sample tau2 
    tau2 <- rinvgamma(1, a1 + P/2, gamma2 + sum(beta^2/(2*lambda2)) )
    
    # 4. Sample lambda2
    #lambda2 <- rinvgamma(P, a3 + .5, psi2 + beta**2/(2*tau2))
    D_inv <- diag(1/(tau2*lambda2))
    
    # 5.  Sample gamma2
    gamma2 <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/tau2)
    
    # 6. Sample psi2
    #psi2 <- stats::rgamma(P, shape = a3 + a4, rate = 1/b2 +1/lambda2)
    
    # save each 'store' iterations 
    if (t%%store == 0){
      #update store counter
      storecount <- storecount + 1
      #store draws
      beta_STORE[storecount,] <- beta
      #lambda2_STORE[storecount,] <- lambda2
      tau2_STORE[storecount,] <- tau2
      gamma2_STORE[storecount,] <- gamma2
      #psi2_STORE[storecount,] <- psi2
      
      # prediction performance checks
      Z_predcheck1 <- unlist(lapply(1:M,function(i){
        Z_i <- c(X[i,,]%*%beta + rnorm(N))
        Z_i[1] <- 0
        #rank of actually observed event in this predicted data
        (N + 1 - rank(Z_i))[Events[i,3]]
      }))
      # below the random part of Z is omitted.
      Z_predcheck2 <- unlist(lapply(1:M,function(i){
        Z_i <- c(X[i,,]%*%beta)
        Z_i[1] <- 0
        #rank of actually observed event in this predicted data
        (N + 1 - rank(Z_i))[Events[i,3]]
      }))
      
      predcheck1_STORE[storecount,] <- Z_predcheck1
      predcheck2_STORE[storecount,] <- Z_predcheck2
      
      #a13_horseshoe_model18 <- beta_hs[c(1:t),]
      #colnames(a13_horseshoe_model18) <- effect_names
      #save(a13_horseshoe_model18, file = paste0("a13_hs_model18_chunk_", t, ".RData"))
      #print("Saved chunk .. ")
    }
    ###################
    
    pb$tick()
    Sys.sleep(1/Nsample)
  }
  
  colnames(beta_STORE) <- dimnames(X)[[3]]
  
  return(list(beta = beta_STORE, predcheck = list(predcheck1_STORE,predcheck2_STORE), 
              #lambda2 = lambda2_STORE, 
              gamma2 = gamma2_STORE, 
              #psi2 = psi2_STORE, 
              tau2 = tau2_STORE,Z = Z))
}

# Gibbs sampler under ridge prior
lasso.dyad2 <- function(Events, X, Nsample = 500, a1 = .5, a2 = 0.5,
                        b1 = 1, store = 10, burnin = 10){
  # D
  print("Define pre-iteration variables ...")
  # Define dimensions
  M <- nrow(Events) # number of events
  N <- dim(X)[2] # number of all possible dyads
  P <- dim(X)[3] # number of covariates without intercept
  effect_names <- dimnames(X)[[3]]
  
  # Define the vectors to store the results
  beta_STORE <- matrix(0,nrow = Nsample/store, ncol = P)
  lambda2_STORE <- gamma2_STORE <- matrix(0,nrow = Nsample/store, ncol = 1)
  tau2_STORE <- matrix(0,nrow = Nsample/store, ncol = P)
  predcheck1_STORE <- predcheck2_STORE <- matrix(0,nrow = Nsample/store, ncol = M)
  
  # initial values
  tau2 <- rep(1, P)
  lambda2 <- gamma2 <- 1
  D <- diag(tau2*lambda2)
  D_inv <- diag(1/(tau2*lambda2))
  beta <- rep(0,P)
  
  #initial computation for posterior cov matrix beta
  # time0 <- Sys.time()
  # XSX <- lapply(1:M,function(i){ #S added in name to better capture the meaning of the object
  #   t(X[i,,])%*%X[i,,]	
  # })
  # XSX_sum <- Reduce("+", XSX)
  # Sys.time() - time0
  
  #a stacked covariates matrix seems faster to compute with than the array
  Xstack <- lapply(1:M,function(i){X[i,,]})
  Xstack <- do.call(rbind,Xstack)
  #time0 <- Sys.time()
  XSX_sum <- t(Xstack)%*%Xstack
  #Sys.time() - time0
  
  # Initial values of Z's which are NOT USED
  Z <- matrix(1,M,N)
  for (i in 1:M){ #JM: would be better to use lapply in case Z is large
    Z[i,-Events[i,3]]<- rtruncnorm(N-1,a=rep(-Inf,N-1),b=rep(0,N-1),mean=0,sd=1)
    Z[,1] <- 0  
  }
  
  print("Start burnin ... ")
  # Sampler iterations ----
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = burnin, clear = F, width = 80)
  for (t in 1:burnin){
    # 1. Sample Z
    #time0 <- Sys.time()
    Z <- generate.indlatent.lapply.dyad.stack(Events, Xstack, beta, Z_max=apply(Z,1,max))
    #Sys.time() - time0
    
    # 2. sample beta
    #time0 <- Sys.time()
    XSZ_sum <- c(t(c(t(Z))) %*% Xstack)
    #Sys.time() - time0
    
    var_beta <- solve(XSX_sum + D_inv)
    mu_beta <- var_beta %*% XSZ_sum
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))
    
    # 3. Sample 1/tau^2 
    mu_tau <- sqrt(lambda2/beta^2)
    tau2 <- 1/rinvgauss(P, mean = mu_tau, shape = 1) 
    
    # 4. Sample lambda2 & gamma2
    lambda2 <- rinvgamma(1, a1 + P/2, gamma2 + sum(beta^2/(2*tau2)) )
    gamma2 <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/lambda2)
    
    # update conditional prior covariance matrix
    D_inv <- diag(1/(tau2 * lambda2))
    
    pb$tick()
    Sys.sleep(1/burnin)
  }
  
  print("Start iterations ... ")
  # Sampler iterations ----
  pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                         total = Nsample, clear = F, width = 80)
  storecount <- 0
  for (t in 1:Nsample){
    # 1. Sample Z
    Z <- generate.indlatent.lapply.dyad.stack(Events, Xstack, beta, Z_max=apply(Z,1,max))
    
    # 2. sample beta
    XSZ_sum <- c(t(c(t(Z))) %*% Xstack)
    var_beta <- solve(XSX_sum + D_inv)
    mu_beta <- var_beta %*% XSZ_sum
    beta <- c(rmvnorm(1, mean = mu_beta, sigma = var_beta))
    
    # 3. Sample 1/tau^2 
    mu_tau <- sqrt(lambda2/beta^2)
    tau2 <- 1/rinvgauss(P, mean = mu_tau, shape = 1)
    
    # 4. Sample lambda2 & gamma2
    lambda2 <- rinvgamma(1, a1 + P/2, gamma2 + sum(beta^2/(2*tau2)) )
    gamma2 <- stats::rgamma(1, shape = a1 + a2, rate = 1/b1 + 1/lambda2)
    
    # update conditional prior covariance matrix
    D_inv <- diag(1/(tau2 * lambda2))
    
    # save each 'store' iterations 
    if (t%%store == 0){
      #update store counter
      storecount <- storecount + 1
      #store draws
      beta_STORE[storecount,] <- beta
      lambda2_STORE[storecount,] <- lambda2
      tau2_STORE[storecount,] <- tau2
      gamma2_STORE[storecount,] <- gamma2
      #psi2_STORE[storecount,] <- psi2
      
      # prediction performance checks
      Z_predcheck1 <- unlist(lapply(1:M,function(i){
        Z_i <- c(X[i,,]%*%beta + rnorm(N))
        Z_i[1] <- 0
        #rank of actually observed event in this predicted data
        (N + 1 - rank(Z_i))[Events[i,3]]
      }))
      # below the random part of Z is omitted.
      Z_predcheck2 <- unlist(lapply(1:M,function(i){
        Z_i <- c(X[i,,]%*%beta)
        Z_i[1] <- 0
        #rank of actually observed event in this predicted data
        (N + 1 - rank(Z_i))[Events[i,3]]
      }))
      
      predcheck1_STORE[storecount,] <- Z_predcheck1
      predcheck2_STORE[storecount,] <- Z_predcheck2
      
      #a13_horseshoe_model18 <- beta_hs[c(1:t),]
      #colnames(a13_horseshoe_model18) <- effect_names
      #save(a13_horseshoe_model18, file = paste0("a13_hs_model18_chunk_", t, ".RData"))
      #print("Saved chunk .. ")
    }
    ###################
    
    pb$tick()
    Sys.sleep(1/Nsample)
  }
  
  colnames(beta_STORE) <- dimnames(X)[[3]]
  
  return(list(beta = beta_STORE, predcheck = list(predcheck1_STORE,predcheck2_STORE), 
              lambda2 = lambda2_STORE, 
              gamma2 = gamma2_STORE, 
              #psi2 = psi2_STORE, 
              tau2 = tau2_STORE,Z = Z))
}



