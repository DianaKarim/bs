library(truncnorm)
library(extraDistr)
library(statmod)
library(mvtnorm)

library(parallel)
library(reshape2)
library(ggplot2)

# Data generation:

N <- 15 # number of actors
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

step <- 7 #step can be D (for R1) ... N (for R2)
R <- c(R2[1:step],rep(2,M-step))
#R

#Estimation

Nsample_test <- 10000
burn_in_test <- 2000

#########################################################################################################
#### All 4 models in parallel ----
#########################################################################################################


cl<-makeCluster(3)
clusterExport(cl, c("R2", "S", "R", "X", "D", "M", "imp", "lasso.new", "ridge","ridge.new" , "horseshoe", "Nsample_test", 
                    "burn_in_test", "Sigma", "generate.latent.lapply"))
clusterEvalQ(cl, {library(truncnorm)
  library(extraDistr)
  library(statmod)
  library(mvtnorm)})

imp_res <-unlist(parLapply(cl, D:M, function(step){
  R <- c(R2[1:step],rep(2,M-step))
  res <- imp(S, R, X, Sigma, Nsample = Nsample_test, burn_in = burn_in_test)
  mean(res$Beta)
} ))

ridge_res <-unlist(parLapply(cl, D:M, function(step){
  R <- c(R2[1:step],rep(2,M-step))
  res <- ridge.new(S, R, X, Sigma, Nsample = Nsample_test, burn_in = burn_in_test)
  mean(res$Beta)
} ))


lasso_res <-unlist(parLapply(cl, D:M, function(step){
  R <- c(R2[1:step],rep(2,M-step))
  res <- lasso.new(S, R, X, Sigma, Nsample = Nsample_test, burn_in = burn_in_test)
  mean(res$Beta)
} ))


hs_res <-unlist(parLapply(cl, D:M, function(step){
  R <- c(R2[1:step],rep(2,M-step))
  res <- horseshoe(S, R, X, Sigma, Nsample = Nsample_test, burn_in = burn_in_test)
  mean(res$Beta)
} ))

stopCluster(cl)

# Save the results ----

save(list = c("imp_res", "ridge_res", "lasso_res", "hs_res"), file = "results_simple_exmpl_correct.Rdata")


# Vizz ----

plot(imp_res, type = "l")
lines(lasso_res, type = "l", col = "green")
lines(hs_res, col = "red")
lines(ridge_res, type = "l", col = "blue")

plot(lasso_res, type = "l", col = "green")


plot(imp_res, lasso_res, type = "l", col = "green")
lines(imp_res, ridge_res, type = "l", col = "blue")
lines(imp_res, hs_res, type = "l", col = "red")


### ggplot ---- 
estim_all <- data.frame(index = c(1:length(c(M:D))),improper = imp_res, ridge = ridge_res, 
                        lasso = lasso_res, horseshoe=hs_res)

estim_all_dt  <-  melt(estim_all, id=c("index"))

#final plot with all estimates:
ggplot(estim_all_dt) + geom_line(aes(x=index, y=value, colour=variable)) +
  scale_colour_manual(values=c("black","blue","green", "red")) +
  labs( x = "Dataset Index", y = "Estimates")

# All estimates vs improper 
df <- melt(estim_all[,-1], id = c("improper"))
ggplot(df) + geom_line(aes(x=improper, y=value, colour=variable)) +
  scale_colour_manual(values=c("blue","green", "red")) +
  labs(x = "Improper prior estimates", y = "Estimates")








