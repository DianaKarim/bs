setwd("E:/Users/mlbosman/Desktop/Revision Diana/run_models_team4")

source("all_functions_dyad.R")

library(truncnorm)
library(extraDistr)
library(mvtnorm)
library(progress)
library(matlib)
library(statmod)

# Load the team4 data
load("team4_estimation_data.RData")

# Statistics in the model
statsi <- 1:dim(stats)[[3]]

# Filename
filename <- "team4_config1"

# Nsample
Nsample <- 100000

# burnin
burnin <- 10000

# Flat prior 
system.time(
	team4_result_flat <- tryCatch(
		flat.dyad2(Events_team4, std_stats[,,statsi], 
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

# Ridge 
system.time(
	team4_result_ridge <- tryCatch(
		ridge.dyad2(Events_team4, std_stats[,,statsi], 
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

# Lasso
system.time(
	team4_result_lasso <- tryCatch(
		lasso.dyad2(Events_team4, std_stats[,,statsi], 
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

# horseshoe
system.time(
	team4_result_horseshoe <- tryCatch(
		horseshoe.dyad2(Events_team4, std_stats[,,statsi], 
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

# Relevent
MLEfit <- relevent::rem(Events_team4$index, std_stats, estimator = "MLE", timing = "ordinal")

# Save results
save(team4_result_flat, team4_result_ridge, team4_result_lasso, team4_result_horseshoe, MLEfit, 
		 file = paste0("team4_results_", filename, ".Rdata"))