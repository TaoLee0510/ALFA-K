# Load required package
library(fields)

# Set seed for reproducibility
set.seed(123)

# Generate 1D data
n <- 30
x_obs <- sort(runif(n, 0, 10))  # Observed locations (sorted for clarity)
y_obs <- sin(x_obs) + rnorm(n, sd = 0.2)  # Observed values with noise

# Convert x_obs to a matrix (Krig expects a matrix)
x_obs <- matrix(x_obs, ncol = 1)

# Define uncertainty in x (normally distributed)
x_sd <- rep(0.2, n)  # Standard deviation of each x location

# Fit the initial Kriging model
krig_fit <- Krig(x_obs, y_obs)

# Define new locations for prediction
x_new <- matrix(seq(0, 10, length.out = 100), ncol = 1)

# Get standard predictions and standard errors
mean_pred <- predict(krig_fit, x_new)
se_pred <- predictSE(krig_fit, x_new)

# Monte Carlo sampling to propagate input uncertainty
n_sim <- 1000  # Number of samples
posterior_samples <- matrix(NA, nrow = n_sim, ncol = length(x_new))

for (i in 1:n_sim) {
  # Sample perturbed input locations
  x_perturbed <- matrix(x_obs + rnorm(n, mean = 0, sd = x_sd), ncol = 1)
  
  # Fit Kriging model for this sample
  krig_fit_mc <- Krig(x_perturbed, y_obs)
  
  # Predict at new locations
  posterior_samples[i, ] <- predict(krig_fit_mc, x_new)
}

# Compute posterior mean and 95% credible interval
posterior_summary <- apply(posterior_samples, 2, function(x) c(
  mean = mean(x), 
  lower = quantile(x, 0.025), 
  upper = quantile(x, 0.975)
))

# Convert to data frame for visualization
posterior_df <- data.frame(x = x_new, t(posterior_summary))

# Plot results
plot(x_obs, y_obs, pch = 16, col = "blue", xlab = "x", ylab = "Predicted value", 
     main = "Kriging with Input Uncertainty")
lines(x_new, posterior_df$mean, col = "red", lwd = 2)
lines(x_new, posterior_df$lower, col = "black", lty = 2)
lines(x_new, posterior_df$upper, col = "black", lty = 2)
legend("topright", legend = c("Predicted Mean", "95% CI"), 
       col = c("red", "black"), lty = c(1, 2), lwd = c(2, 1))



