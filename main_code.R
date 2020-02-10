# Comparison of the priors
library(truncnorm)
library(quantreg)
library(statmod)
library(ggplot2)
library(reshape2)
library(smoothmest)
library(parallel)
library(tidyverse)
library(extraDistr)


########################################################################################
# 1. Generate output Y ------ 
########################################################################################

# Output is a vector with 75 observations with different beta, that is considered as prob-ty 
#from 0.5 to approx-tly 0.85

n <- 70
sigma2 <- 1 # known variance of data 
y <- matrix(nrow = 25, ncol= n)
for(i in 1:25){
  y[i,] <- c(rep(1, i+34), rep(0,n-i-34))}
apply(y, 1 , sum)
dim(y)

# make a random sample 
for (i in 1:dim(y)[1]) {
  y[i,] <- sample(y[i,], size  = dim(y)[2], replace = F)}

# vector of observed probabilities
apply(y, 1, sum)/n


########################################################################################
# 2. Estimation -----
########################################################################################

## use parallel computation 

y_list <- list()
for (i in 1:dim(y)[1]) {
  y_list[[i]] <- y[i,]
}

# All 4 estimates ()

cl<-makeCluster(3)
clusterExport(cl, c("n", "sigma2"))
clusterEvalQ(cl, {library(truncnorm)
  library(extraDistr)
  library(statmod)})
estim_imp <-  unlist(parLapply(cl, y_list, imp.sample.1, n = 70,  Nsample = 50000, burn_in = 10000))
estim_ridge <-  unlist(parLapply(cl, y_list, ridge.sample.1, n = 70, Nsample = 50000, burn_in = 10000))
estim_lasso <-  unlist(parLapply(cl, y_list, lasso.sample.1, n = 70, Nsample = 50000, burn_in = 10000))
estim_hs <-  unlist(parLapply(cl, y_list, hs.sample.1, n = 70, Nsample = 50000, burn_in = 10000))
stopCluster(cl)



########################################################################################
#3. Plots ----
########################################################################################

# results as data frame 
estim_all <- data.frame(beta = c(1:dim(y)[1]),improper = estim_imp, 
                        ridge = estim_ridge, lasso = estim_lasso, horseshoe=estim_hs)
estim_all_dt  <-  melt(estim_all, id=c("beta"))

#final plot with all estimates:
ggplot(estim_all_dt) + geom_line(aes(x=beta, y=value, colour=variable)) +
  scale_colour_manual(values=c("black","red","blue", "green")) +
  labs( x = "Dataset Index", y = "Estimates")

#ggsave("4estim.pdf", width = 8, height = 6)


## Plot estimates vs improper estimates 

df <- melt(estim_all[,-1], id = c("improper"))
ggplot(df) + geom_line(aes(x=improper, y=value, colour=variable)) +
  scale_colour_manual(values=c("red","blue", "green")) +
  labs(x = "Improper prior estimates", y = "Estimates")
#ggsave("estim_vs_imp.pdf", width = 8, height = 6)



# 4. Compute estimates separately ----- 

# improper
cl<-makeCluster(3)
clusterExport(cl, c("n", "sigma2"))
clusterEvalQ(cl, library(truncnorm))
estim_imp <-  parLapply(cl, y_list, imp.sample.1, n = 70)
stopCluster(cl)
estim_imp <- unlist(estim_imp)
########################################################

# ridge
cl<-makeCluster(3)
clusterExport(cl, c("n", "sigma2"))
clusterEvalQ(cl, library(truncnorm) )
estim_ridge <-  parLapply(cl, y_list, ridge.sample.1, n=70)
stopCluster(cl)
estim_ridge <- unlist(estim_ridge)
########################################################

#lasso
cl<-makeCluster(3)
clusterExport(cl, c("n", "sigma2"))
clusterEvalQ(cl, {library(truncnorm)
  library(statmod) })
estim_lasso <-  parLapply(cl, y_list, lasso.sample.1, n=70)
stopCluster(cl)
estim_lasso <- unlist(estim_lasso)
########################################################

#horseshoe 
cl<-makeCluster(3)
clusterExport(cl, c("n", "sigma2"))
clusterEvalQ(cl, {library(truncnorm)
  library(extraDistr) })
estim_hs <-  parLapply(cl, y_list, hs.sample.1, n = 70)
stopCluster(cl)
estim_hs <- unlist(estim_hs)
########################################################






