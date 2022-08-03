library(MCMCglmm)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(plotly)
library(stringr)

# Load results
filename <- "model2_subset"
load(paste0("apollo_results_", filename, ".Rdata"))



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
check.samples(apollo_result_flat)
check.samples(apollo_result_ridge)
check.samples(apollo_result_lasso)
check.samples(apollo_result_horseshoe) 
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
post.dens(apollo_result_flat, apollo_result_ridge, 
          apollo_result_lasso, apollo_result_horseshoe)
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

resflat <- get.results(apollo_result_flat$beta)
resridge <- get.results(apollo_result_ridge$beta)
reslasso <- get.results(apollo_result_lasso$beta)
reshorseshoe <- get.results(apollo_result_horseshoe$beta)

# Combine
resflat$prior <- "flat"
resridge$prior <- "ridge"
reslasso$prior <- "lasso"
reshorseshoe$prior <- "horseshoe"
res <- rbind(resflat, resridge, reslasso, reshorseshoe)
res$width <- res$upperCI-res$lowerCI

# Make sure prior is a factor
res$prior <- factor(res$prior)
res$prior <- factor(res$prior, levels = c("flat", "ridge", "lasso", "horseshoe"))

# Descriptives
# Significance
table(res$sig, res$prior)
table(res$sig99, res$prior)

# Plotting preparations 
# Determine and label the facets (model 1)
res$facet <- c(rep(1, 12), rep(2:8, times = 13))
res$facet <- factor(res$facet)
res$facet <- factor(res$facet, labels = c("Ground to ground", "Air to air", 
                                          "Ground to CAPCOM", "Air to CAPCOM",
                                          "CAPCOM to ground", "CAPCOM to air", 
                                          "ground to FLIGHT", "FLIGHT to ground"))

# Make sure variable is a factor variable
res$variable <- factor(res$variable)

# Plotting variableNames
res$variableName <- gsub(".x.tie.air_to_air", "", res$variableName)
res$variableName <- gsub(".x.tie.ground_to_CAPCOM", "", res$variableName)
res$variableName <- gsub(".x.tie.air_to_CAPCOM", "", res$variableName)
res$variableName <- gsub(".x.tie.CAPCOM_to_air", "", res$variableName)
res$variableName <- gsub(".x.tie.CAPCOM_to_ground", "", res$variableName)
res$variableName <- gsub(".x.tie.FLIGHT_to_ground", "", res$variableName)
res$variableName <- gsub(".x.tie.ground_to_FLIGHT", "", res$variableName)
res$variableName <- gsub("tie.", "", res$variableName)
res$variableName <- gsub("inertia", "I", res$variableName)
res$variableName <- gsub("reciprocity", "R", res$variableName)
res$variableName <- gsub("indegree", "ID", res$variableName)
res$variableName <- gsub("outdegree", "OD", res$variableName)
res$variableName <- gsub("CAPCOM", "CC", res$variableName)
res$variableName <- gsub("air", "A", res$variableName)
res$variableName <- gsub("ground", "G", res$variableName)
res$variableName <- gsub("Receiver", "Rec", res$variableName)
res$variableName <- gsub("Receive", "Rec", res$variableName)
res$variableName <- gsub("Sender", "Snd", res$variableName)
res$variableName <- gsub("Send", "Snd", res$variableName)
res$variableName <- gsub("FLIGHT", "F", res$variableName)

# Plot results
pdf(file = paste0("apollo_estimates", ".pdf"), width = 7.5, height = 0.8*11)

ggplot(res, aes(x = postMode, y = variable, group = prior)) +
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
        rect = element_rect(fill = "white"))

dev.off()



res_part1 <- res[res$facet %in% c("Ground to ground", "Air to air", "Ground to CAPCOM", "Air to CAPCOM"),]

pdf(file = paste0("apollo_estimates_part1", ".pdf"), width = 7.5, height = 0.8*11)

ggplot(res_part1, aes(x = postMode, y = variable, group = prior)) +
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
        rect = element_rect(fill = "white"))

dev.off()


res_part2 <- res[res$facet %in% c("CAPCOM to ground", "CAPCOM to air", "ground to FLIGHT", "FLIGHT to ground"),]
pdf(file = paste0("apollo_estimates_part2", ".pdf"), width = 7.5, height = 0.8*11)

ggplot(res_part2, aes(x = postMode, y = variable, group = prior)) +
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
        rect = element_rect(fill = "white"))

dev.off()


# Load predictive performance functions
source("predictive_performance_functions.R")

# Load data
load("apollo_estimation_data2.RData")

# In-sample prediction 
inpred_flat <- get.insample.pred(apollo_result_flat$beta, Events_apollo, 
                                 stats2, M = nrow(Events_apollo)-500)$scores
inpred_ridge <- get.insample.pred(apollo_result_ridge$beta, Events_apollo, 
                                  stats2, M = nrow(Events_apollo)-500)$scores
inpred_lasso <- get.insample.pred(apollo_result_lasso$beta, Events_apollo, 
                                  stats2, M = nrow(Events_apollo)-500)$scores
inpred_horseshoe <- get.insample.pred(apollo_result_horseshoe$beta, Events_apollo, 
                                      stats2, M = nrow(Events_apollo)-500)$scores

# Plot
pdf(file = 'apollo_insample_mode.pdf', height = 6, width = 8)
a <- barplots(inpred_flat, inpred_ridge, inpred_lasso, inpred_horseshoe)
a 
dev.off()


# Out-of-sample-prediction
outpred_flat <- get.outofsample.pred(apollo_result_flat$beta, Events_apollo, 
                                     stats2, M = nrow(Events_apollo)-500, n = 500)$scores
outpred_ridge <- get.outofsample.pred(apollo_result_ridge$beta, Events_apollo, 
                                      stats2, M = nrow(Events_apollo)-500, n = 500)$scores
outpred_lasso <- get.outofsample.pred(apollo_result_lasso$beta, Events_apollo, 
                                      stats2, M = nrow(Events_apollo)-500, n = 500)$scores
outpred_horseshoe <- get.outofsample.pred(apollo_result_horseshoe$beta, Events_apollo, 
                                          stats2, M = nrow(Events_apollo)-500, n = 500)$scores


#plot 

pdf(file = 'apollo_outofsample_mode.pdf', height = 6, width = 8)
b <- barplots(outpred_flat, outpred_ridge, outpred_lasso, outpred_horseshoe)
b 
dev.off()






#Draws

inpred_flat_draws <- get.insample.pred.draws(apollo_result_flat$beta, Events_apollo, 
                                 stats2, M = nrow(Events_apollo)-500)$scores
inpred_ridge_draws <- get.insample.pred.draws(apollo_result_ridge$beta, Events_apollo, 
                                  stats2, M = nrow(Events_apollo)-500)$scores
inpred_lasso_draws <- get.insample.pred.draws(apollo_result_lasso$beta, Events_apollo, 
                                  stats2, M = nrow(Events_apollo)-500)$scores
inpred_horseshoe_draws <- get.insample.pred.draws(apollo_result_horseshoe$beta, Events_apollo, 
                                      stats2, M = nrow(Events_apollo)-500)$scores

pdf(file = 'apollo_insample_draws.pdf', height = 6, width = 8)
c <- barplots(inpred_flat_draws, inpred_ridge_draws, inpred_lasso_draws, inpred_horseshoe_draws)
c 
dev.off()


outpred_flat_draws <- get.outofsample.pred.draws(apollo_result_flat$beta, Events_apollo, 
                                     stats2, M = nrow(Events_apollo)-500, n = 500)$scores
outpred_ridge_draws <- get.outofsample.pred.draws(apollo_result_ridge$beta, Events_apollo, 
                                      stats2, M = nrow(Events_apollo)-500, n = 500)$scores
outpred_lasso_draws <- get.outofsample.pred.draws(apollo_result_lasso$beta, Events_apollo, 
                                      stats2, M = nrow(Events_apollo)-500, n = 500)$scores
outpred_horseshoe_draws <- get.outofsample.pred.draws(apollo_result_horseshoe$beta, Events_apollo, 
                                          stats2, M = nrow(Events_apollo)-500, n = 500)$scores

pdf(file = 'apollo_outofsample_draws.pdf', height = 6, width = 8)
d <- barplots(outpred_flat_draws, outpred_ridge_draws, outpred_lasso_draws, outpred_horseshoe_draws)
d 
dev.off()


save(outpred_flat, outpred_ridge, outpred_lasso, outpred_horseshoe, 
     outpred_flat_draws, outpred_ridge_draws, outpred_lasso_draws, outpred_horseshoe_draws, 
     file = "apollo_outpred_scores_50K.RData")

save(inpred_flat, inpred_ridge, inpred_lasso, inpred_horseshoe,
     inpred_flat_draws, inpred_ridge_draws, inpred_lasso_draws, inpred_horseshoe_draws,
     file = "apollo_inpred_scores_50K.RData")

