#setwd("E:/Users/mlbosman/Desktop/Revision Diana/run_models")

library(ggplot2)
library(ggthemes)
library(gridExtra)
library(plotly)
library(stringr)
library(MCMCglmm)

# Load results
filename <- "main_config1_2000events_100K"
load(paste0("enron_results_", filename, ".RData"))

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

check.samples(enron_result_flat)
check.samples(enron_result_ridge)
check.samples(enron_result_lasso)
check.samples(enron_result_horseshoe)

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

post.dens(enron_result_flat, enron_result_ridge, 
          enron_result_lasso, enron_result_horseshoe)

# Obtain posterior mean, upper and lower 95% CI and significance
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
  results$variable <- 1:nrow(results)
  results$variableName <- rownames(results)
  results <- results[,c(5,6,1:4)]
  rownames(results) <- NULL
  results
}

resflat <- get.results(enron_result_flat$beta)
resridge <- get.results(enron_result_ridge$beta)
reslasso <- get.results(enron_result_lasso$beta)
reshorseshoe <- get.results(enron_result_horseshoe$beta)

# Combine
resflat$prior <- "flat"
resridge$prior <- "ridge"
reslasso$prior <- "lasso"
reshorseshoe$prior <- "horseshoe"
res <- rbind(resflat, resridge, reslasso, reshorseshoe)
res$width <- res$upperCI-res$lowerCI
res$prior <- factor(res$prior)
res$prior <- factor(res$prior, levels = c("flat", "ridge", "lasso", "horseshoe"))

# Descriptives
# Significance
table(res$sig, res$prior)

# Save the result file ?

#res_original <- res
#save(res_original, file = "enron_res_original.RData")

# Confidence intervals
summary(res$width[res$prior == "flat"])
summary(res$width[res$prior == "ridge"])
summary(res$width[res$prior == "lasso"])
summary(res$width[res$prior == "horseshoe"]) 

# Plotting preparations 
# Determine and label the facets 
res$facet <- rep(1:8, each = 8)
res$facet <- factor(res$facet, labels = c("1 day", "2 days",
                                          "1 week", "2 weeks", 
                                          "1 month", "3 months",
                                          "Exogenous effects (1)", "Exogenous effects (2)"))

# Make sure variable is a factor variable
res$variable <- factor(res$variable)

# Replace variable names
res$variableName <- gsub("tie.", "", res$variableName)
res$variableName <- gsub("24", "", res$variableName)
res$variableName <- gsub("48", "", res$variableName)
res$variableName <- gsub("168", "", res$variableName)
res$variableName <- gsub("336", "", res$variableName)
res$variableName <- gsub("744", "", res$variableName)
res$variableName <- gsub("2184", "", res$variableName)
res$variableName <- gsub("outdegreeReceiver", "popularity", res$variableName)
res$variableName <- gsub("indegreeReceiver", "activity", res$variableName)
res$variableName <- gsub("inertia", "I", res$variableName)
res$variableName <- gsub("reciprocity", "R", res$variableName)
res$variableName <- gsub("activity", "A", res$variableName)
res$variableName <- gsub("popularity", "P", res$variableName)
res$variableName <- gsub("otp", "OTP", res$variableName)
res$variableName <- gsub("osp", "OSP", res$variableName)
res$variableName <- gsub("itp", "ITP", res$variableName)
res$variableName <- gsub("isp", "ISP", res$variableName)
res$variableName <- gsub("LegalLegal", "LL", res$variableName)
res$variableName <- gsub("TradingLegal", "TL", res$variableName)
res$variableName <- gsub("JuniorLegal", "JL", res$variableName)
res$variableName <- gsub("FemaleLegal", "FL", res$variableName)
res$variableName <- gsub("LegalTrading", "LT", res$variableName)
res$variableName <- gsub("TradingTrading", "TT", res$variableName)
res$variableName <- gsub("FemaleTrading", "FT", res$variableName)
res$variableName <- gsub("JuniorTrading", "JT", res$variableName)
res$variableName <- gsub("LegalJunior", "LJ", res$variableName)
res$variableName <- gsub("TradingJunior", "TJ", res$variableName)
res$variableName <- gsub("JuniorJunior", "JJ", res$variableName)
res$variableName <- gsub("FemaleJunior", "FJ", res$variableName)
res$variableName <- gsub("LegalFemale", "LF", res$variableName)
res$variableName <- gsub("TradingFemale", "TF", res$variableName)
res$variableName <- gsub("FemaleFemale", "FF", res$variableName)
res$variableName <- gsub("JuniorFemale", "JF", res$variableName)

# Plot results
pdf(file = "enron_estimates.pdf", height = 10.5, width = 7.5)

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

# Prediction performance:

# Get predictive performance (within sample)
source("predictive_performance_functions.R")

# In-sample 

#load the stats for config 1 
load("enron_std_stats_config1.RData") 
stats_std <- enron_std_stats_config1
load("Events_enron.RData")
# load here the extended stats


# In sapmle 
in_sample_flat <- get.insample.pred(enron_result_flat$beta, Events_enron, stats_std, M = 2000)$scores
in_sample_ridge <- get.insample.pred(enron_result_ridge$beta, Events_enron, stats_std, M = 2000)$scores
in_sample_lasso <- get.insample.pred(enron_result_lasso$beta, Events_enron, stats_std, M = 2000)$scores
in_sample_horseshoe <- get.insample.pred(enron_result_horseshoe$beta, Events_enron, stats_std, M = 2000)$scores


pdf(file = 'enron_insample_mode.pdf', height = 6, width = 8)
a <- barplots(in_sample_flat, in_sample_ridge, in_sample_lasso, 
              in_sample_horseshoe, 
              plot_title = "")
a
dev.off()


# In-sample prediction based on posterior draws: 

system.time(
  in_sample_draws_flat <- get.insample.pred.draws(enron_result_flat$beta, Events_enron, stats_std, M = 2000)$scores
)

in_sample_draws_ridge <- get.insample.pred.draws(enron_result_ridge$beta, Events_enron, stats_std, M = 2000)$scores
in_sample_draws_lasso <- get.insample.pred.draws(enron_result_lasso$beta, Events_enron, stats_std, M = 2000)$scores
in_sample_draws_horseshoe <- get.insample.pred.draws(enron_result_horseshoe$beta, Events_enron, stats_std, M = 2000)$scores

pdf(file = "enron_insample_draws.pdf", height = 6, width = 8)
b <- barplots(in_sample_draws_flat, in_sample_draws_ridge, in_sample_draws_lasso, 
              in_sample_draws_horseshoe, 
              plot_title = "") 
b
dev.off()


# Out of sample prediction: 


out_of_sample_flat <- get.outofsample.pred(enron_result_flat$beta, Events_enron, stats_std, M = 2000, n = 2000)$scores
out_of_sample_ridge <- get.outofsample.pred(enron_result_ridge$beta, Events_enron, stats_std, M = 2000, n = 2000)$scores
out_of_sample_lasso <- get.outofsample.pred(enron_result_lasso$beta, Events_enron, stats_std, M = 2000, n = 2000)$scores
out_of_sample_horseshoe <- get.outofsample.pred(enron_result_horseshoe$beta, Events_enron, stats_std, M = 2000, n = 2000)$scores

#Plot
pdf(file = 'enron_outofsample_sample_mode.pdf', height = 6, width = 8)
c <- barplots(out_of_sample_flat, out_of_sample_ridge, out_of_sample_lasso, out_of_sample_horseshoe,
              plot_title = "")
c
dev.off()
# Out of sample prediction by posterior draws
# n next events to predict 

out_of_sample_draws_flat <- get.outofsample.pred.draws(enron_result_flat$beta, Events_enron, stats_std, M=2000, n=2000)$scores
out_of_sample_draws_ridge <- get.outofsample.pred.draws(enron_result_ridge$beta, Events_enron, stats_std, M=2000, n=2000)$scores
out_of_sample_draws_lasso <- get.outofsample.pred.draws(enron_result_lasso$beta, Events_enron, stats_std, M=2000, n=2000)$scores
out_of_sample_draws_horseshoe <- get.outofsample.pred.draws(enron_result_horseshoe$beta, Events_enron, stats_std, M=2000, n=2000)$scores

# plot: 

pdf(file = "enron_outofsample_draws.pdf", height = 6, width = 8)
d <- barplots(out_of_sample_draws_flat, out_of_sample_draws_ridge, 
              out_of_sample_draws_lasso, out_of_sample_draws_horseshoe,
              plot_title = "")
d
dev.off()

# Only 5% prediction score 
modepred <- rbind(c(in_sample_flat$five, in_sample_ridge$five, 
                    in_sample_lasso$five, in_sample_horseshoe$five),
                  c(out_of_sample_flat$five, out_of_sample_ridge$five, 
                    out_of_sample_lasso$five, out_of_sample_horseshoe$five))
rownames(modepred) <- c("insample", "outofsample")
colnames(modepred) <- c("flat", "ridge", "lasso", "horseshoe")

drawspred <- rbind(c(in_sample_draws_flat$five, in_sample_draws_ridge$five, 
                     in_sample_draws_lasso$five, in_sample_draws_horseshoe$five),
                   c(out_of_sample_draws_flat$five, out_of_sample_draws_ridge$five, 
                     out_of_sample_draws_$five, out_of_sample_draws_horseshoe$five))
rownames(trendpred) <- c("insample", "outofsample")
colnames(trendpred) <- c("flat", "ridge", "lasso", "horseshoe")

round(modepred*100, 2)
round(drawspred*100, 2)


# 99% credible interval 

get.results99 <- function(x) {
  results <- apply(x, 2, function(y) {
    data.frame(postMode = posterior.mode(y), 
               lowerCI = quantile(y, 0.005),
               upperCI = quantile(y, 0.995))
  })
  results <- do.call(rbind, results)
  results$sig <- apply(results, 1, function(y) {
    !dplyr::between(0, y[2], y[3])
  })
  results$variable <- 1:nrow(results)
  results$variableName <- rownames(results)
  results <- results[,c(5,6,1:4)]
  rownames(results) <- NULL
  results
}


resflat99 <- get.results99(enron_result_flat$beta)
resridge99 <- get.results99(enron_result_ridge$beta)
reslasso99 <- get.results99(enron_result_lasso$beta)
reshorseshoe99 <- get.results99(enron_result_horseshoe$beta)

# Combine
resflat99$prior <- "flat99"
resridge99$prior <- "ridge99"
reslasso99$prior <- "lasso99"
reshorseshoe99$prior <- "horseshoe99"
res99 <- rbind(resflat99, resridge99, reslasso99, reshorseshoe99)
res99$width <- res99$upperCI-res$lowerCI
res99$prior <- factor(res99$prior)
res99$prior <- factor(res99$prior, levels = c("flat99", "ridge99", "lasso99", "horseshoe99"))

# Significance
table(res99$sig, res99$prior)


res <- rbind(res, res99)

tab <- table(res$variable, res$sig)
dt  <- data.frame(variable = c(1:64), nsig = tab[,2])
var_order <- dt[order(dt$nsig, decreasing = T), 1]

ord <- var_order
for (i in c(1:7)) {
  ord <- c(ord, c(var_order+64*i))
}

res$prior <- factor(res$prior, levels = c("flat", "flat99", 
                                            "ridge", "ridge99", 
                                             "lasso", "lasso99", 
                                            "horseshoe", "horseshoe99"))

temp <- res[ord, ]
temp$variable <- as.character(temp$variable)



pdf( file = 'dotplot95_99_reduced.pdf', height = 7, width = 5.5)

ggplot(data = temp[temp$sig == T,], aes(x = prior, y = variable)) + 
  geom_point(shape = 19) +
  scale_y_discrete(breaks = as.character(var_order), labels = as.character(var_order)) + 
  scale_x_discrete(breaks=c("flat", "flat99", 
                            "ridge", "ridge99", 
                            "lasso", "lasso99", 
                            "horseshoe", "horseshoe99"),
                   labels=c("flat, 95% (27)", "flat, 99% (16)", 
                            "ridge, 95% (24)", "ridge, 99% (16)", 
                            "lasso, 95% (22)", "lasso, 99% (15)", 
                            "horseshoe, 95% (18)", "horseshoe, 99% (14)")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5), 
        axis.text.y = element_text(size = 7)) + ylab("Effect") + xlab("")

dev.off()










