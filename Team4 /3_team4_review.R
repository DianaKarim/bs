setwd("E:/Users/mlbosman/Desktop/Revision Diana/run_models_team4")

library(MCMCglmm)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(plotly)
library(stringr)

# Load results
filename <- "team4_config1"
load(paste0("team4_results_", filename, ".Rdata"))

# Check samples
check.samples <- function(x) {
	beta <- x$beta
	par(mfrow=c(4,4))
	sapply(1:ncol(beta), function(i) {
		plot(beta[,i], type = "l", ylab = paste("stat", i))
		traceMean <- sapply(1:NROW(beta[,i]), function(j) mean(beta[1:j,i]))
		lines(traceMean, col = "red")
	})
	par(mfrow=c(1,1))
}

pdf(file = "checksamples.pdf")
check.samples(team4_result_flat)
check.samples(team4_result_ridge)
check.samples(team4_result_lasso)
check.samples(team4_result_horseshoe) 
dev.off()

# Posterior densities
post.dens <- function(x1, x2, x3, x4) {
	beta1 <- x1$beta
	beta2 <- x2$beta
	beta3 <- x3$beta
	beta4 <- x4$beta
	par(mfrow=c(4,4))
	sapply(1:ncol(beta1), function(i) {
		d1 <- density(beta1[,i])
		d2 <- density(beta2[,i])
		d3 <- density(beta3[,i])
		d4 <- density(beta4[,i])
		xlim <- c(min(d1$x, d2$x, d3$x, d4$x), max(d1$x, d2$x, d3$x, d4$x))
		ylim <- c(min(d1$y, d2$y, d3$y, d4$y), max(d1$y, d2$y, d3$y, d4$y))
		plot(density(beta1[,i]), ylab = paste("stat", i), xlim = xlim, ylim = ylim)
		abline(v = mean(beta1[,i]))
		lines(density(beta2[,i]), col = "red")
		abline(v = mean(beta2[,i]), col = "red")
		lines(density(beta3[,i]), col = "green")
		abline(v = mean(beta3[,i]), col = "green")
		lines(density(beta4[,i]), col = "blue")
		abline(v = mean(beta4[,i]), col = "blue")
	})
	par(mfrow=c(1,1))
}

pdf(file = "postdensities.pdf")
post.dens(team4_result_flat, team4_result_ridge, 
					team4_result_lasso, team4_result_horseshoe)
dev.off()

# Obtain posterior mode, upper and lower 95% CI and significance
get.results <- function(x) {
	results <- apply(x, 2, function(y) {
		data.frame(postMode = posterior.mode(y), 
							 lowerCI = quantile(y, 0.025),
							 upperCI = quantile(y, 0.975),
							 lowerCI99 = quantile(y, 0.005),
							 upperCI99 = quantile(y, 0.995))
	})
	results <- do.call(rbind, results)
	results$sig <- apply(results, 1, function(y) {
		!dplyr::between(0, y[2], y[3])
	})
	results$sig99 <- apply(results, 1, function(y) {
		!dplyr::between(0, y[4], y[5])
	})
	results$variable <- 1:nrow(results)
	results$variableName <- rownames(results)
	results <- results[,c(8,9,1:7)]
	rownames(results) <- NULL
	results
}

resflat <- get.results(team4_result_flat$beta)
resridge <- get.results(team4_result_ridge$beta)
reslasso <- get.results(team4_result_lasso$beta)
reshorseshoe <- get.results(team4_result_horseshoe$beta)

# Combine
resflat$prior <- "flat"
resridge$prior <- "ridge"
reslasso$prior <- "lasso"
reshorseshoe$prior <- "horseshoe"
res <- rbind(resflat, resridge, reslasso, reshorseshoe)
res$width <- res$upperCI-res$lowerCI

# Descriptives
# Significance
table(res$sig, res$prior)
table(res$sig99, res$prior)

# Plotting preparations 
# Determine and label the facets 
res$facet <- c(rep(c(rep(1, 4), rep(2, 8), rep(4, 8), rep(6, 8), rep(3, 8), rep(5, 8)),4))
res$facet <- factor(res$facet, labels = c("Endogenous effects", "Exogenous effects",
																					"Interactions with 'degree'", "Interactions with 'inertia'", 
																					"Interactions with 'recency'", "Interactions with 'shared partners'"))

# Make sure variable is a vactor variable
res$variable <- factor(res$variable)

# Replace variable names
res$variableName <- gsub("tie.", "", res$variableName)
res$variableName <- gsub("totaldegreeDyad.x.", "", res$variableName)
res$variableName <- gsub("inertia.x.", "", res$variableName)
res$variableName <- gsub("recencyContinue.x.", "", res$variableName)
res$variableName <- gsub("sp.x.", "", res$variableName)
res$variableName <- gsub("same.gender", "SS", res$variableName)
res$variableName <- gsub("difference.age", "AD", res$variableName)
res$variableName <- gsub("same.hqual", "SE", res$variableName)
res$variableName <- gsub("difference.tenure", "TD", res$variableName)
res$variableName <- gsub("diffSocial", "SD", res$variableName)
res$variableName <- gsub("diffAdvice", "ASD", res$variableName)
res$variableName <- gsub("meanSocial", "SA", res$variableName)
res$variableName <- gsub("meanAdvice", "ASA", res$variableName)
res$variableName <- gsub("recencyContinue", "R", res$variableName)
res$variableName <- gsub("totaldegreeDyad", "D", res$variableName)
res$variableName <- gsub("inertia", "I", res$variableName)
res$variableName <- gsub("sp", "SP", res$variableName)
res$variable[res$variableName == "D"] <- 1
res$variable[res$variableName == "I"] <- 2
res$variable[res$variableName == "R"] <- 3
res$variable[res$variableName == "SP"] <- 4


# Add the 99% flat as a different model
# resflat99 <- subset(res, prior == "flat")
# resflat99$sig <- resflat99$sig99
# resflat99$prior <- "flat 99%"
# resflat99$lowerCI <- resflat99$lowerCI99
# resflat99$upperCI <- resflat99$upperCI99
# res$prior <- gsub("flat", "flat 95%", res$prior)
# res <- rbind(res, resflat99)
# res$prior <- factor(res$prior)
# res$prior <- factor(res$prior, levels = c("flat 95%", "flat 99%", "ridge", "lasso", "horseshoe"))

# Plot results
pdf(file = paste0("team4estimates", ".pdf"), width = 7.5, height = 0.75*11)

a <- ggplot(res, aes(x = postMode, y = variable, group = prior)) +
	geom_errorbar(aes(xmin = lowerCI, xmax = upperCI, col = sig), 
								width = 0.75, position = "dodge") +
	geom_point(aes(shape = prior, color = sig), position = ggstance::position_dodgev(height = 0.75)) +
	geom_vline(xintercept = 0, linetype = "dashed") +
	scale_x_continuous(name = "Posterior mode") +
	scale_y_discrete(name = element_blank(), breaks = (res$variable), 
									 labels = (res$variableName), limits = rev) +
	scale_color_manual(name = "Significance", values = c("darkgrey", "black"), guide = "none") +
	scale_shape_manual(name = "Prior", labels = c("Flat", "Ridge", "Lasso", "Horsehoe"),
							values = c(16, 17, 15, 4)) +
	facet_wrap(~ facet, scales = "free", ncol = 2) +
	theme_bw() +
	theme(legend.position = "bottom", axis.text.y = element_text(angle = 0),
				text = element_text(size = 12),
				rect = element_rect(fill = "white")); a

dev.off()

# Load predictive performance functions
source("predictive_performance_functions.R")

# Load data
load("team4_estimation_data.RData")

# In-sample prediction 
inpred_flat <- get.insample.pred(team4_result_flat$beta, Events_team4, 
																 std_stats, M = nrow(Events_team4))$scores
inpred_ridge <- get.insample.pred(team4_result_ridge$beta, Events_team4, 
																	std_stats, M = nrow(Events_team4))$scores
inpred_lasso <- get.insample.pred(team4_result_lasso$beta, Events_team4, 
																	std_stats, M = nrow(Events_team4))$scores
inpred_horseshoe <- get.insample.pred(team4_result_horseshoe$beta, Events_team4, 
																			std_stats, M = nrow(Events_team4))$scores

# Plot
a <- barplots(inpred_flat, inpred_ridge, inpred_lasso, inpred_horseshoe)

# In-sample prediction based on posterior draws
inpredd_flat <- get.insample.pred.draws(team4_result_flat$beta, Events_team4, 
																 std_stats, M = nrow(Events_team4))$scores
inpredd_ridge <- get.insample.pred.draws(team4_result_ridge$beta, Events_team4, 
																	std_stats, M = nrow(Events_team4))$scores
inpredd_lasso <- get.insample.pred.draws(team4_result_lasso$beta, Events_team4, 
																	std_stats, M = nrow(Events_team4))$scores
inpredd_horseshoe <- get.insample.pred.draws(team4_result_horseshoe$beta, Events_team4, 
																			std_stats, M = nrow(Events_team4))$scores

# Save results
save(inpredd_flat, inpredd_ridge, inpredd_lasso, inpredd_horseshoe, file = "in_sample_prediction_draws_team4.RData")

# Plot
b <- barplots(inpredd_flat, inpredd_ridge, inpredd_lasso, inpredd_horseshoe);b

# Out-sample prediction
# Load data
load("team4_out_of_sample_prediction_data.RData")

# Number of events to predict
n <- 500

# Out-of-sample-prediction
outpred_flat <- get.outofsample.pred(team4_result_flat$beta, Events_team4, 
																 std_stats, M = length(secondDay), n = n)$scores
outpred_ridge <- get.outofsample.pred(team4_result_ridge$beta, Events_team4, 
																	std_stats, M = length(secondDay), n = n)$scores
outpred_lasso <- get.outofsample.pred(team4_result_lasso$beta, Events_team4, 
																	std_stats, M = length(secondDay), n = n)$scores
outpred_horseshoe <- get.outofsample.pred(team4_result_horseshoe$beta, Events_team4, 
																			std_stats, M = length(secondDay), n = n)$scores

# Plot
c <- barplots(outpred_flat, outpred_ridge, outpred_lasso, outpred_horseshoe)

# out-sample prediction based on posterior draws
outpredd_flat <- get.outofsample.pred.draws(team4_result_flat$beta, Events_team4, 
																				std_stats, M = length(secondDay), n = n)$scores
outpredd_ridge <- get.outofsample.pred.draws(team4_result_ridge$beta, Events_team4, 
																				 std_stats, M = length(secondDay), n = n)$scores
outpredd_lasso <- get.outofsample.pred.draws(team4_result_lasso$beta, Events_team4, 
																				 std_stats, M = length(secondDay), n = n)$scores
outpredd_horseshoe <- get.outofsample.pred.draws(team4_result_horseshoe$beta, Events_team4, 
																						 std_stats, M = length(secondDay), n = n)$scores

# Save results
save(outpredd_flat, outpredd_ridge, outpredd_lasso, outpredd_horseshoe, file = "out_sample_prediction_draws_team4.RData")

# Plot
d <- barplots(outpredd_flat, outpredd_ridge, outpredd_lasso, outpredd_horseshoe);d

grid.arrange(a,b,c,d,ncol = 2)

# Plots for in paper
pdf(file = "team4_insample.pdf", width = 11/2, height = 8.5/2)
b
dev.off()
pdf(file = "team4_outsample.pdf", width = 11/2, height = 8.5/2)
d
dev.off()

# Plots for in supmat
pdf(file = "team4_insample_postmode.pdf", width = 11/2, height = 8.5/2)
a
dev.off()
pdf(file = "team4_outsample_postmode.pdf", width = 11/2, height = 8.5/2)
c
dev.off()
