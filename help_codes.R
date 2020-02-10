# 5. Illustration of priors ----

xx <- seq(-4,4,length.out = 1000)
y_rigde <- dnorm(xx)
y_lasso <- ddoublex(xx)
y_cauchy <- dcauchy(xx)

df <- data.frame(xx, y_rigde, y_lasso, y_cauchy)
df_melt  <-  melt(df, id=c("xx"))

ggplot(df_melt) + geom_line(aes(x=xx, y=value, colour=variable)) +
  scale_color_manual(values = c("red", "blue", "purple"))
#ggsave("2priors.jpeg", width = 6, height = 4)





#########################################################################################
#scaled plots
plot(pnorm(estim_imp), type = "l" )
lines(pnorm(estim_ridge), col = "red")
lines(pnorm(estim_lasso), col = "blue")
lines(pnorm(estim_hs), col = "green")
#lines(pnorm(estim_cauchy), col = "purple")

estim_all_scaled <- data.frame(beta = c(1:dim(output)[1]),improper = pnorm(estim_imp), 
                               ridge = pnorm(estim_ridge), lasso = pnorm(estim_lasso),
                               horseshoe = pnorm(estim_hs))
estim_all_scaled  <-  melt(estim_all_scaled, id=c("beta"))

ggplot(estim_all_scaled) + geom_line(aes(x=beta, y=value, colour=variable)) +
  scale_colour_manual(values=c("black","red","blue", "green"))
ggsave("scaled1.pdf")

ggplot(estim_all_scaled) + geom_smooth(aes(x=beta, y=value, colour=variable)) +
  scale_colour_manual(values=c("black","red","blue", "green"))
ggsave("scaled2.pdf")



plot(estim_imp, estim_ridge, type = "l")
lines(estim_imp, estim_lasso, col = "red")
lines(estim_imp, estim_hs, col = "green")


plot(estim_ridge, type = "l")
plot(estim_lasso, type = "l")
plot(estim_hs, type = "l")
#plot(estim_cauchy, type = "l")



# Sampling the parameters within 1 sample 
plot(beta_sample_hs, type = "l")
plot(tau2, type = "l")
plot(gma, type = "l")
plot(lambda2, type = "l")
plot(psi2, type = "l")



