# Sampler for Bayesian Lasso model with fixed lasso parameter 

library("mvtnorm")
library("truncnorm")
library("quantreg")
library("statmod")
library("progress")
library("coda")

## Set wd to source file location
rm(list=ls()) #remove all
#rm(list = setdiff(ls(), lsf.str())) #remove all except function 


#############################################################################################################
# Dimensions of the network ----
#############################################################################################################

N <- 1000   # number of observations [1000]
J <- 100  # number of actors (network size) [100]
K <- 20  # number of covariates (excuding intercept) [20]

#############################################################################################################
####### Data generation ----
#############################################################################################################

Sigma <- create.sigma(sigma2 = 2.5, rho = 0.5, d = J)
inv_sigma <- solve(Sigma)

S <- vector()
R <- vector()
set.seed(9999)
for (i in (1:N)) {
  S[i] <- sample(1:J, 1, replace = T)} 

set.seed(9999)
R[1] <- sample(c(1:J)[-S[1]], 1)

#beta <- c(1,rep(1, 3), 1,1,1,1,0, 1,1,1,0,0, 0,0,0,0,0)
#beta <- c(1,rep(1, 3), 1,1,1,1,1, 0,0,0,0,0, 1,1,1,0,0)
#beta <- c(1, 1,1,0, 1,1,1,1,1, rep(0, 5), 1,1,1,0,0)

beta <- c(1, 1,1,0, 1,1,1, rep(0, 7), 1,1,0,0,0,0,0)



set.seed(9999)
my_data <- gen.data.with.stats.v2(S,R,beta,Sigma, J, N)

S <- my_data$senders
R <- my_data$receivers
X <- my_data$cov

#############################################################################################################
####### Initials ----
#############################################################################################################
# Prior variance of beta
tau0 <- 10^6 # variance of intercept
tau <- c(tau0, rep(100, K))
D_inv <- diag(tau^(-1))

# Initial values of beta
set.seed(2323)
#beta0 <-rep(0, times = K+1)
beta0 <- runif(K+1, 0, 1)

# Initial value for Z 
Z_old <- matrix(0,N,J)
for (i in 1:N){
  Z_old[i,S[i]] <- NA
  Z_old[i,R[i]] <- ifelse(R[i]==1,0,1)
  Z_old[i,-c(R[i],S[i])] <- rtruncnorm(J-2,a=rep(-Inf,J-2),b=rep(Z_old[R[i]],J-2),mean=0,sd=1)
  Z_old[,1] <- 0
}



#############################################################################################################
####### Preliminary calculations ----
#############################################################################################################

for (i in 1:N){
  X[[i]] <- as.matrix(X[[i]])
  X[[i]][S[i],] <- NA 
}

XX <- list()
for (i in 1:N) {
  XX[[i]] <- t(X[[i]][-S[i],])%*%inv_sigma[-S[i],-S[i]]%*%X[[i]][-S[i],]}

XX_sum <- Reduce("+", XX) ## sum(tilde(Sigma))_i
xsigmaz <- list() 

# Helper function's index
ind <- seq(1,N,1)


###############
#JM: The function below will be faster and correctly incorporates that S[i] cannot be the receiver of message i.
#initial stuff for generate.latent.lapply
Sigma11 <- Sigma[-1,-1] #note that this is the matrix we have to work with because the sender does not have a Z
Sigma121 <- Sigma11[1,-1]%*%solve(Sigma11[-1,-1]) #needed in every step
sigma2_j <- c(Sigma11[1,1] - t(Sigma11[1,-1])%*%solve(Sigma11[-1,-1])%*%Sigma11[-1,1]) #needed in every step


#############################################################################################################
# Sampler   ---- 
#############################################################################################################

Nsample <- 5000
burn_in <- 2000

# Hyperparameters
delta <- vector(length = Nsample+1)
delta[1] <- 0.5 #(rate)
r <- 3 #shape



# Beta initials 
beta_sample3 <- matrix(nrow = Nsample+1, ncol = K+1)
beta_sample3[1,] <- beta0 # fixing initial value of beta

# Lasso parameter
lambda2 <- rep(0, Nsample+1)
lambda2[1] <- 1
# inverse_D_sample <- list() # list for storing the presicion matrix of beta
# inverse_D_sample[[1]] <- D_inv

tau_sample <- rep(0, K+1) 


pb <- progress_bar$new(format = "(:spin) [:bar] [:percent]",
                       total = Nsample, clear = F, width = 80)
system.time(
  for (t in 1:Nsample){
    # 1. Generate (sample) Z
    beta_current <- beta_sample3[t,]
    Z <- generate.latent.lapply(S, R, X, beta=beta_current, Z=Z_old, Sigma121, sigma2_j)
    
    # 2. generate (sample) beta
    #var_beta <- solve(XX_sum + inverse_D_sample[[t]])
    var_beta <- solve(XX_sum + D_inv)
    var_beta[lower.tri(var_beta)] = t(var_beta)[lower.tri(var_beta)]
    xsigmaz <- lapply(ind, xsigmaz.lasso)
    sum_xsigmaz <- Reduce("+",xsigmaz)
    mu_beta <- var_beta%*%sum_xsigmaz
    
    beta1 <- rmvnorm(1, mean = mu_beta, sigma = var_beta, method = "eigen")
    
    # 3. Generate  1/tau^2 
    
    mu_tau <- sqrt((lambda2[t])/beta1[-1]^2)
    
    tau_sample <- unlist(lapply(1:K,function(k){ 
      rinvgauss(1, mean = mu_tau[k], shape = lambda2[t])
    }))
    
    #inverse_D_sample[[t+1]] <- diag(c(tau0^(-1), tau_sample))  
    D_inv <-  diag(c(tau0^(-1), tau_sample))  
    
    # 4.  Generate lambda^2
    
    delta[t+1] <- delta[1] + 1/2*sum(1/tau_sample)
    lambda2[t+1] <- rgamma(1, shape = r+K, rate = delta[t+1])
    
    
    # 4. Save the intermediary values
    beta_sample3[t+1,] <- beta1
    #Z_sample[[t]] <- Z
    Z_old <- Z  
    pb$tick()
    Sys.sleep(1/Nsample)})


c_names <- colnames(X[[1]])
c_names[2] <- "popularity"
c_names[3] <- "reciprocity"
c_names[4] <- "inertia"
colnames(beta_sample3) <- c_names

save(beta_sample3, file = "beta_sample3_6.09_newb.RData")


#############################################################################################################
# Results and Graphs   ---- 
#############################################################################################################

# Dynamics of the sample

jpeg("Fullbayes_dynamics.jpg", width = 1200, height = 1000)
par(mfrow=c(4,5), oma = c(3,1,3,1))
for (i in 2:(K+1)){
  plot(beta_sample3[,i], type = "l", xlab = paste0(c_names[i], "=", beta[i]), cex.lab = 2.5, cex.axis = 1.5, ylab= "", main = "", 
       ylim =  c(min(beta_sample3[,i], beta[i]),max(beta_sample3[,i], beta[i])))
  
  #(h=beta[i], col = "red1", lwd = 0.5) #true value of coefficients
  abline(h = 0)
  abline(v=2000, type="l", lty=2)}
par(mfrow=c(1,1))
#title(main = "Lasso Prior Model", outer = T)
dev.off()



jpeg("Fullbayes_dynamics_intercept.jpg", width = 900, height = 700)
par(oma = c(5,1,3,1))
  plot(beta_sample3[,1], type = "l", main = "Intercept", xlab = "", ylab="", 
       ylim =  c(min(beta_sample3[,1], beta[1]),max(beta_sample3[,1], beta[1])) )
  abline(h=beta[1], col = "red1", lwd = 0.5) #true value of coefficients
  #  legend("center", legend = "True value", lty = 1, col = "red", box.lty = 0)
title(main = "Full Bayes Model, Sample dynamics",
      sub = paste0(" Prior parameters: r = ", r, ", delta = ", delta[1]), outer = T)
dev.off()


# Reduced sample (after burn-in)

jpeg("Fullbayes_dynamics_reduced.jpg", width = 1200, height = 1000)
par(mfrow=c(4,5), oma = c(5,1,5,1))
for (i in 2:(K+1)){
  plot(beta_reduced3[,i], type = "l", main = c_names[i], 
       ylim = c(min(beta0, beta_reduced3[,i]),max(beta0, beta_reduced3[,i])), 
       xlab = "", ylab = "")
  abline(h=beta[i], col = "red1", lwd = 0.5) #true value of coefficients
  abline(h=beta0[i], col = "green", lwd = 0.5) #initial values of the sample
}
par(mfrow=c(1,1))
title(main = "Full Bayes Model, Sample dynamics after burn-in", 
      sub = paste0(" Prior parameters: r = ", r, ", delta = ", delta[1]), outer = T )
dev.off()

jpeg("Fullbayes_dynamics_reduced_intercept.jpg", width = 900, height = 700)
par(oma = c(5,1,5,1))
plot(beta_reduced3[,1], type = "l", main = "Intercept", 
       ylim = c(min(beta0[1], beta_reduced3[,1]),max(beta0[1], beta_reduced3[,1])), 
       xlab = "", ylab = "")
  abline(h=beta[1], col = "red1", lwd = 0.5) #true value of coefficients
  abline(h=beta0[1], col = "green", lwd = 0.5) #initial values of the sample
title(main = "Full Bayes Model, Sample dynamics after burn-in", 
      sub = paste0(" Prior parameters: r = ", r, ", delta = ", delta[1]), outer = T )
dev.off()



# Sampled posterior densities 

jpeg("Fullbayes_densities.jpg", width = 1200, height = 1000)
par(mfrow=c(4,5), oma = c(8, 1, 5, 1))
for (i in 2:(K+1)){
  plot(density(beta_reduced3[,i]), type = "l", xlab = c_names[i],
       xlim = c(min(beta[i], beta0[i], beta_reduced3[,i]), max(beta[i], beta0[i], beta_reduced3[,i])), 
       ylab = '', main = "")
  q <- seq(min(beta[i], beta0[i], beta_reduced3[,i]), max(beta[i], beta0[i], beta_reduced3[,i]), length.out = 500)
  points(q, dnorm(q, mean = 0, sd = tau[i]), type = "l", col = "blue")
  abline(v = beta[i], col = "red1") #true value
  abline(v = beta0[i], col = "green") #initial value of the sampler
}
par(mfrow=c(1,1))
title(main = paste0("Full Bayes Model with lasso prior parameters: r = ", r, ", delta = ", delta[1]), 
      sub = "Posterior density(black), prior density (blue), true value(red) and starting value(green)", outer = T)
dev.off()

# Intercept
jpeg("Fullbayes_densities_intercept.jpg", width = 1000, height = 800)
par(oma = c(8, 1, 5, 1))
plot(density(beta_reduced3[,1]), type = "l", xlab = "Intercept",
       xlim = c(min(beta[1], beta0[1], beta_reduced3[,1]), max(beta[1], beta0[1], beta_reduced3[,1])), 
       ylab = '', main = "")
  q <- seq(min(beta[1], beta0[1], beta_reduced3[,1]), max(beta[i], beta0[1], beta_reduced3[,1]), length.out = 500)
  points(q, dnorm(q, mean = 0, sd = tau[i]), type = "l", col = "blue")
  abline(v = beta[1], col = "red1") #true value
  abline(v = beta0[1], col = "green") #initial value of the sampler
title(main = paste0("Full Bayes Model with lasso prior parameters: r = ", r, ", delta = ", delta[1]), 
      sub = "Posterior density(black), prior density (blue), true value(red) and starting value(green)", outer = T)
dev.off()



## Lambda 

jpeg("Fullbayes_lamdba_densities.jpg", width = 600, height = 400)
plot(density(lambda2[c((burn_in+1):Nsample)]), type = "l", lwd = 2, ylim = c(0,1),
     xlab="", ylab = "", main = "" )
q <- seq(0, 1.1*max(lambda2), length.out = 500)
points(q, dgamma(q, shape = r, rate = delta[1]), type = "l", lty = 'dashed',lwd = 2, col = "darkblue")
legend("topright", legend = c("Prior","Posterior"), lty = c(2,1), col = c("darkblue","black"))
#title(sub =  paste0(" Prior parameters: r = ", r, ", delta = ", delta[1]))
dev.off()

# Sample dynamics
jpeg("Fullbayes_lamdba_sample.jpg", width = 600, height = 800)
par(mfrow = c(2,1))
plot(lambda2, type  = "l", main = "Sample values of Lasso parameter", ylab = "", xlab = "")
abline(h = lambda2[1], lty = "dashed", col = "green", lwd = 2)
plot(lambda2[c((burn_in+1):Nsample)], type  = "l", lwd = 1, 
     main = "Sample values of Lasso parameters, after burn-in", ylab = "", xlab = "" )
abline(h = lambda2[1], lty = "dashed", col = "green", lwd = 2)
#legend("topright", legend = "Starting value of the sample", col = "green", lty = 2, box.lty = 0, outer = T)
par(mfrow = c(1,1))
dev.off()

mean(lambda2[c((burn_in+1):Nsample)]) # posterior mean 
mean(delta[c((burn_in+1):Nsample)]) # posterior rate

#------------------------------------------------------------------------------#

effectiveSize(beta_sample3)

unlist(lapply(1:(K+1),function (i) {
  return(effectiveSize(beta_sample3[,i]))
}))
      

effectiveSize(beta_reduced3)
