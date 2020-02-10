
# Code for analysis of latent variables Z
library(ggplot2)
library(tidyr)

# Data generation: -----

N <- 6 # number of actors
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
step <- 6 #step can be D (for R1) ... N (for R2)
R <- c(R2[1:step],rep(2,M-step))


# Analyse the behaviour of latent variables Z 
imp_res1 <- imp(S, R, X, Sigma, Nsample = 100, burn_in = 20)
Z_temp <- matrix(unlist(imp_res1$Z), ncol = N, byrow = T)
head(Z_temp, 10)


Z_temp <- as.data.frame(Z_temp)
Z_temp <- cbind(iteration = c(1:100), Z_temp)

Z_temp

Z_temp %>%
  gather(key, value, V1, V2, V3, V4, V5, V6) %>%
  ggplot(aes(x=iteration, y=value, colour=key)) +
  geom_line()

