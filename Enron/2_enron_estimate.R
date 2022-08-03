setwd("E:/Users/mlbosman/Desktop/Revision Diana")
source("all_functions_actor.R")
#setwd("E:/Users/mlbosman/Desktop/Revision Diana/run_models")

# Load the enron data
load("enron_estimation_data.RData")

# Events
maxEvents <- 2000  # number of event in the subset of the data, used for the estimation

# Statistics
# all main endogenous statistics plus the "tie" exogenous stats
statsi <- c(1:48, 60:75)
filename <- "main_config1_2000events" 

# Nsample
Nsample <- 10000

# burnin
burnin <- 1000

# Flat prior 
system.time(
	enron_result_flat <- tryCatch(
		flat.actor(Events_enron[1:maxEvents,], std_stats[1:maxEvents,,statsi], 
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
	enron_result_ridge <- tryCatch(
		ridge.actor(Events_enron[1:maxEvents,], std_stats[1:maxEvents,,statsi], 
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
	enron_result_lasso <- tryCatch(
		lasso.actor(Events_enron[1:maxEvents,], std_stats[1:maxEvents,,statsi], 
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
	enron_result_horseshoe <- tryCatch(
		horseshoe.actor(Events_enron[1:maxEvents,], std_stats[1:maxEvents,,statsi], 
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
save(enron_result_flat, enron_result_ridge, enron_result_lasso, enron_result_horseshoe,
		 file = paste0("enron_results_", filename, ".Rdata"))