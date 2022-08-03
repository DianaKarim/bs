library(truncnorm)
library(extraDistr)
library(mvtnorm)
library(progress)
library(matlib)
library(statmod)

source("all_functions_dyad.R")

# Load data
load("apollo_estimation_data2.RData")

# Filename
filename <- "model2_subset"

# Stats
stats <- stats2

# Subset
events <- 1:(nrow(Events_apollo)-500)

# Nsample
Nsample <- 10000

# burnin
burnin <- 1000

# Flat prior 
system.time(
  apollo_result_flat <- tryCatch(
    flat.dyad2(Events_apollo[events,], stats[events,,], 
               Nsample, store = 10, burnin = burnin),
    error = function(cond) {
      message(cond)
      return(NA)
    },
    warning = function(cond) {
      message(cond)
      return(NA)
    }
  )
)

# Save results
save(apollo_result_flat,   
     file = paste0("apollo_results_", filename, ".Rdata"))


# Ridge 
system.time(
  apollo_result_ridge <- tryCatch(
    ridge.dyad2(Events_apollo[events,], stats[events,,],
                Nsample, a1 = .5, a2 = .5, b1 = 1, store = 10, burnin = burnin),
    error = function(cond) {
      message(cond)
      return(NA)
    },
    warning = function(cond) {
      message(cond)
      return(NA)
    }
  )
)

# Save results
save(apollo_result_flat, apollo_result_ridge, 
     file = paste0("apollo_results_", filename, ".Rdata"))


# Lasso
system.time(
  apollo_result_lasso <- tryCatch(
    lasso.dyad2(Events_apollo[events,], stats[events,,],
                Nsample = Nsample, a1 = .5, a2 = .5, b1 = 1, store = 10, burnin = burnin),
    error = function(cond) {
      message(cond)
      return(NA)
    },
    warning = function(cond) {
      message(cond)
      return(NA)
    }
  )
)

# Save results
save(apollo_result_flat, apollo_result_ridge, apollo_result_lasso,  
     file = paste0("apollo_results_", filename, ".Rdata"))


# horseshoe
system.time(
  apollo_result_horseshoe <- tryCatch(
    horseshoe.dyad2(Events_apollo[events,], stats[events,,],
                    Nsample, a1 = .5, a2 = .5, a3 = .5, a4 = .5, b1 = 1, b2 = 1, store = 10, burnin = burnin),
    error = function(cond) {
      message(cond)
      return(NA)
    },
    warning = function(cond) {
      message(cond)
      return(NA)
    }
  )
)

# Save results
save(apollo_result_flat, apollo_result_ridge, apollo_result_lasso, apollo_result_horseshoe,  
     file = paste0("apollo_results_", filename, ".Rdata"))
