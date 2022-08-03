library(MCMCglmm)

# In sample prediction based on posterior mode: 
get.insample.pred <- function(beta, Events, stat, M){
  # Get posterior mode
  hat_beta <- apply(beta, 2, posterior.mode)
  
  # Get temporary data.frame with event ranks and indicators for ranks in the 
  # set
  temp <- matrix(unlist(lapply(1:M, function(i){
    Z <- stat[i,,]%*%hat_beta 
    Z_ordered <- order(Z, decreasing = T)
    return(cbind(which(Z_ordered == Events[i, 3]), which.max(Z) == Events[i, 3],
                 sum(Z_ordered[1:floor(0.05*length(Z))] == Events[i, 3]), 
                 sum(Z_ordered[1:floor(0.1*length(Z))] == Events[i, 3]), 
                 sum(Z_ordered[1:floor(0.2*length(Z))] == Events[i, 3])))
  })), nrow = M, byrow = T)
  
  # Summarize the temporary data.frame
  res = data.frame("strict" = sum(temp[,2])/M,
                   "five" = sum(temp[,3])/M, 
                   "ten" = sum(temp[,4])/M,
                   "twenty" = sum(temp[,5])/M)
  
  # Output
  list(rank_observed_mean = mean(temp[,1]),
       rank_observed_vector = temp[,1],
       scores = res) 
}

# In-sample prediction based on posterior draws: 
get.insample.pred.draws <- function(beta, Events, stat, M){
  
  Nsample <- nrow(beta)
  res <- matrix(0, nrow = Nsample, ncol = 5)
  
  for (j in 1:Nsample){
    temp <- matrix(unlist(lapply(1:M, function(i){
      Z <- stat[i,,]%*%beta[j,] 
      Z_ordered <- order(Z, decreasing = T)
      return(cbind(which(Z_ordered == Events[i, 3]), 
                   which.max(Z) == Events[i, 3],
                   sum(Z_ordered[1:floor(0.05*length(Z))] == Events[i, 3]), 
                   sum(Z_ordered[1:floor(0.1*length(Z))] == Events[i, 3]), 
                   sum(Z_ordered[1:floor(0.2*length(Z))] == Events[i, 3])))
    })), nrow = M, byrow = T)
    res[j,] <- c(mean(temp[,1]), sum(temp[,2]), sum(temp[,3]), sum(temp[,4]), sum(temp[,5]))
  }
  
  rs <- data.frame("strict" = mean(res[,2])/M, 
                   "five" = mean(res[,3])/M,
                   "ten" = mean(res[,4])/M,
                   "twenty" = mean(res[,5]/M))
  
  return(list(rank_observed_mean = mean(res[,1]),
              scores = rs))  
  
}

# Out of sample prediction: 
# Note: M is the number of events the posterior draws are based on, n
# is the next number of events the out of sample prediction is computed for
get.outofsample.pred <- function(beta, Events, stat, M, n){
  # Get posterior mode
  hat_beta <- apply(beta, 2, posterior.mode)
  
  temp <- matrix(unlist(lapply(1:n, function(i){
    X_i <- stat[(M+i),,]
    Z <- X_i%*%hat_beta 
    Z_ordered <- order(Z, decreasing = T)
    return(cbind(which(Z_ordered == Events[M+i, 3]), 
                 which.max(Z) == Events[M+i, 3],
                 sum(Z_ordered[1:floor(0.05*length(Z))] == Events[M+i, 3]), 
                 sum(Z_ordered[1:floor(0.1*length(Z))] == Events[M+i, 3]), 
                 sum(Z_ordered[1:floor(0.2*length(Z))] == Events[M+i, 3])))
  })), nrow = n, byrow = T)
  
  
  res = data.frame("strict" = round(sum(temp[,2]), 4),
                   "five" = round(sum(temp[,3]), 4), 
                   "ten" = round(sum(temp[,4]), 4),
                   "twenty" = round(sum(temp[,5]), 4))
  
  return(list(rank_observed_mean = mean(temp[,1]),
              rank_observed_vector = temp[,1],
              scores = res/n))  
}

# Out of sample prediction: 
# Note: M is the number of events the posterior draws are based on, n
# is the next number of events the out of sample prediction is computed for
get.outofsample.pred.trend <- function(beta, Events, stat, M, n){
  # Get posterior mode
  hat_beta <- apply(beta, 2, posterior.mode)
  
  # Running score
  scores <- sapply(10:n, function(j) {
    nT <- j
    
    temp <- matrix(unlist(lapply(1:nT, function(i){
      X_i <- stat[(M+i),,]
      Z <- X_i%*%hat_beta 
      Z_ordered <- order(Z, decreasing = T)
      return(cbind(which(Z_ordered == Events[M+i, 3]), 
                   which.max(Z) == Events[M+i, 3],
                   sum(Z_ordered[1:floor(0.05*length(Z))] == Events[M+i, 3]), 
                   sum(Z_ordered[1:floor(0.1*length(Z))] == Events[M+i, 3]), 
                   sum(Z_ordered[1:floor(0.2*length(Z))] == Events[M+i, 3])))
    })), nrow = nT, byrow = T)
    
    
    res <- data.frame("strict" = round(sum(temp[,2]), 4),
                      "five" = round(sum(temp[,3]), 4), 
                      "ten" = round(sum(temp[,4]), 4),
                      "twenty" = round(sum(temp[,5]), 4))
    
    data.frame(res/nT)
  })
  
  scores <- matrix(unlist(scores), byrow = T, ncol = 4)
  names(scores) <- c("strict", "five", "ten", "twenty")
  scores
}

# Out of sample prediction by posterior draws
# n next events to predict 
get.outofsample.pred.draws <- function(beta, Events, stat, M, n){
  
  Nsample <- nrow(beta)
  res <- matrix(0, nrow = Nsample, ncol = 5)
  
  for (j in 1:Nsample){
    
    temp <- matrix(unlist(lapply(1:n, function(i){
      Z <- stat[(M+i),,]%*%beta[j,] 
      Z_ordered <- order(Z, decreasing = T)
      return(cbind(which(Z_ordered == Events[M+i, 3]), 
                   which.max(Z) == Events[M+i, 3],
                   sum(Z_ordered[1:floor(0.05*length(Z))] == Events[M+i, 3]), 
                   sum(Z_ordered[1:floor(0.1*length(Z))] == Events[M+i, 3]), 
                   sum(Z_ordered[1:floor(0.2*length(Z))] == Events[M+i, 3])))
    })), nrow = n, byrow = T)
    
    res[j,] <- cbind(rank_observed_mean = mean(temp[,1])/n,
                     prediction_strict = sum(temp[,2])/n,
                     prediction_5 = sum(temp[,3])/n, 
                     prediction_10 = sum(temp[,4])/n,
                     prediction_20 = sum(temp[,5])/n)
  }
  
  rs = data.frame("strict" = mean(res[,2]),
                  "five" = mean(res[,3]), 
                  "ten" = mean(res[,4]),
                  "twenty" = mean(res[,5]))
  
  return(list(rank_observed_mean = mean(res[,1]),
              scores = rs))  
  
}

barplots <- function(flat, ridge, lasso, horseshoe, plot_title = element_blank()) {
  # Add prior
  flat$prior <- "flat"
  ridge$prior <- "ridge"
  lasso$prior <- "lasso"
  horseshoe$prior <- "horseshoe"
  
  # Prepare data.frame for plotting
  bars_data <- rbind(flat, ridge, lasso, horseshoe)
  bars_data$prior <- factor(bars_data$prior, levels = c("flat", "ridge", "lasso", "horseshoe"))
  bars_data <- reshape2::melt(bars_data, id = "prior")
  bars_data$percent <- bars_data$value * 100
  bars_data <- bars_data[bars_data$variable %in% c("five", "ten", "twenty"),]
 
  
  # Check whether to plot 5%
  if(all(bars_data$value[bars_data$variable == "strict"] == 
         bars_data$value[bars_data$variable == "five"]))
  {
    bars_data <- subset(bars_data, variable != "strict")
  }
  
  
  # Plot
  ggplot(data = bars_data, aes(x = variable, y = percent, fill = prior, group = prior)) +
    geom_col(position = "dodge") + 
    geom_text(aes(x = variable, y = percent, label = sprintf("%.1f", round(percent, digits = 1)), group = prior), 
              position = position_dodge(width = 1), hjust = -0.1, size = 4, angle = 90) +
    scale_x_discrete(breaks = c("five", "ten", "twenty"),
                     labels = c("5%", "10%", "20%")) +
    scale_fill_grey(name = "Prior", labels = c("Flat", "Ridge", "Lasso", "Horseshoe")) +
    labs(x = "", y = "Percent", title = plot_title) + 
    coord_cartesian(ylim = c(0,100)) +
    theme_classic() +
    theme(legend.position = c(0,1), legend.justification = c(0,1),
          legend.background = element_blank(),
          legend.direction = 'horizontal',
          text = element_text(size = 16), legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
}






# Legend at the bottom
barplots <- function(flat, ridge, lasso, horseshoe, plot_title = element_blank()) {
  # Add prior
  flat$prior <- "flat"
  ridge$prior <- "ridge"
  lasso$prior <- "lasso"
  horseshoe$prior <- "horseshoe"
  
  # Prepare data.frame for plotting
  bars_data <- rbind(flat, ridge, lasso, horseshoe)
  bars_data$prior <- factor(bars_data$prior, levels = c("flat", "ridge", "lasso", "horseshoe"))
  bars_data <- reshape2::melt(bars_data, id = "prior")
  bars_data$percent <- bars_data$value * 100
  bars_data <- bars_data[bars_data$variable %in% c("five", "ten", "twenty"),]
  
  
  # Check whether to plot 5%
  if(all(bars_data$value[bars_data$variable == "strict"] == 
         bars_data$value[bars_data$variable == "five"]))
  {
    bars_data <- subset(bars_data, variable != "strict")
  }
  
  
  # Plot
  ggplot(data = bars_data, aes(x = variable, y = percent, fill = prior, group = prior)) +
    geom_col(position = "dodge") + 
    geom_text(aes(x = variable, y = percent, label = sprintf("%.1f", round(percent, digits = 1)), group = prior), 
              position = position_dodge(width = 1), hjust = -0.1, size = 4, angle = 90) +
    scale_x_discrete(breaks = c("five", "ten", "twenty"),
                     labels = c("5%", "10%", "20%")) +
    scale_fill_grey(name = "Prior", labels = c("Flat", "Ridge", "Lasso", "Horseshoe")) +
    labs(x = "", y = "Percent", title = plot_title) + 
    coord_cartesian(ylim = c(0,102)) +
    theme_classic() +
    # theme(legend.position = c(0,1), legend.justification = c(0,1),
    #     legend.background = element_blank(),
    #     legend.direction = 'horizontal',
    #     text = element_text(size = 16), legend.text = element_text(size = 12),
    #     legend.title = element_text(size = 12))
    theme(legend.position = "bottom",
          legend.background = element_blank(),
          legend.direction = 'horizontal',
          text = element_text(size = 16), legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))
  
}

