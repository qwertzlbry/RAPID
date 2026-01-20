
################################################################################
# Real Data Application: UCI Adult Dataset
# Using the new RAPID implementation with normalized gain (RCS_marginal)
################################################################################

# Load the RAPID package
devtools::load_all()

# Load required packages
library(synthpop)      # For synthetic data generation (CART)
library(ranger)        # For Random Forest
library(boot)          # For bootstrapped confidence intervals
library(dplyr)         # For data wrangling
library(ggplot2)       # For plotting
library(tidyr)         # For data reshaping

################################################################################
# 1. Load and prepare UCI Adult dataset
################################################################################

url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data"
df <- read.csv(url, header = FALSE, strip.white = TRUE)

# Set column names as per UCI documentation
colnames(df) <- c(
  "age", "workclass", "fnlwgt", "education", "education.num", "marital.status",
  "occupation", "relationship", "race", "sex", "capital.gain", "capital.loss",
  "hours.per.week", "native.country", "income"
)

# Convert character columns to factors
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)], factor)

# Drop unneeded or redundant variables
df <- df %>% select(-fnlwgt, -education)

# Remove rows with missing values in the sensitive attribute
df <- df[complete.cases(df$income), ]

cat("Dataset loaded:", nrow(df), "rows,", ncol(df), "columns\n")

################################################################################
# 2. Generate synthetic datasets using CART
################################################################################

cat("Generating 5 synthetic datasets using CART...\n")
synths <- synthpop::syn(df, m = 5, method = "cart")
synthetic_list <- synths$syn

################################################################################
# 3. Define function to compute RAPID using new implementation
################################################################################

compute_rapid_new <- function(original_data,
                              synthetic_data,
                              y = "income",
                              tau = 0.3,        # normalized gain threshold
                              model = "rf",
                              n_boot = 500) {

  quasi_ids <- setdiff(colnames(original_data), y)

  # Compute RAPID using RCS_marginal method (normalized gain)
  rapid_result <- rapid(
    original_data = original_data,
    synthetic_data = synthetic_data,
    quasi_identifiers = quasi_ids,
    sensitive_attribute = y,
    model_type = model,
    cat_eval_method = "RCS_marginal",  # Use normalized gain method
    cat_tau = tau,
    trace = FALSE
  )

  # Bootstrap for confidence intervals
  boot_fn <- function(data, indices) {
    tmp_orig <- data[indices, ]
    tmp_result <- rapid(
      original_data = tmp_orig,
      synthetic_data = synthetic_data,
      quasi_identifiers = quasi_ids,
      sensitive_attribute = y,
      model_type = model,
      cat_eval_method = "RCS_marginal",
      cat_tau = tau,
      trace = FALSE
    )
    return(tmp_result$risk$confidence_rate)
  }

  boot_vals <- boot::boot(data = original_data, statistic = boot_fn, R = n_boot)
  ci <- boot::boot.ci(boot_vals, type = "perc")$percent[4:5]

  return(list(
    rapid = rapid_result$risk$confidence_rate,
    ci = ci,
    n_at_risk = rapid_result$risk$n_at_risk
  ))
}

################################################################################
# 4. Compute RAPID across range of tau values for threshold curve
################################################################################

cat("Computing RAPID across tau values for threshold curve...\n")

taus <- seq(0.0, 1.0, by = 0.15) # Reduced for speed; decrease to 0.05 for final
R <- n_boot_curve <- 10  # Reduced for speed; increase to 500 for final
mymodel <- "rf" # Used for speed; use "rf" for final.

# Compute RAPID for each tau and each synthetic replicate
rapid_matrix <- matrix(NA, nrow = length(taus), ncol = length(synthetic_list))

for (i in seq_along(taus)) {
  tau_val <- taus[i]
  cat("  tau =", tau_val, "\n")

  for (m in seq_along(synthetic_list)) {
    result <- compute_rapid_new(
      original_data = df,
      synthetic_data = synthetic_list[[m]],
      y = "income",
      tau = tau_val,
      model = mymodel,
      n_boot = R  # minimal bootstrap for speed
    )
    rapid_matrix[i, m] <- result$rapid
  }
}

# Create summary data frame
rapid_df <- data.frame(
  tau = taus,
  rapid_mean = rowMeans(rapid_matrix),
  rapid_q05 = apply(rapid_matrix, 1, quantile, probs = 0.05),
  rapid_q95 = apply(rapid_matrix, 1, quantile, probs = 0.95),
  rapid_min = apply(rapid_matrix, 1, min),
  rapid_max = apply(rapid_matrix, 1, max)
)

################################################################################
# 5. Create threshold curve plot
################################################################################

cat("Creating threshold curve plot...\n")

p <- ggplot(rapid_df, aes(x = tau, y = rapid_mean)) +
  geom_ribbon(aes(ymin = rapid_min, ymax = rapid_max),
              fill = "darkorange", alpha = 0.3) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "red", linewidth = 0.5) +
  annotate("text", x = 0.32, y = max(rapid_df$rapid_mean) * 0.9,
           label = expression(tau == 0.3), hjust = 0, color = "red") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(
    title = expression("RAPID vs Normalized Gain Threshold (" * tau * ")"),
    subtitle = "Mean and range across 5 synthetic replicates (UCI Adult dataset)",
    x = expression(paste("Normalized Gain Threshold ", tau)),
    y = "RAPID (proportion at risk)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# Save plot
output_path <- "/Users/matthias/workspace26/RAPID/inst/doc/paper/rapid_curve_quantile.pdf"
pdf(output_path, width = 8, height = 6)
print(p)
dev.off()

cat("Plot saved to:", output_path, "\n")

# Also display
print(p)

################################################################################
# 6. Summary statistics at default threshold
################################################################################

R <- 5 # increase to at least 500
cat("\n=== Summary at default threshold tau = 0.3 ===\n")

default_tau <- 0.3
default_row <- rapid_df[rapid_df$tau == default_tau, ]

cat("Mean RAPID:", round(default_row$rapid_mean, 4), "\n")
cat("Range:", round(default_row$rapid_min, 4), "-", round(default_row$rapid_max, 4), "\n")

# Compute with proper bootstrap at default threshold
cat("\nComputing RAPID with full bootstrap at tau = 0.3...\n")

rapid_at_default <- lapply(synthetic_list, function(syn_df) {
  compute_rapid_new(
    original_data = df,
    synthetic_data = syn_df,
    y = "income",
    tau = 0.3,
    model = "rf",
    n_boot = R
  )
})

# Summary table
summary_table <- data.frame(
  Replicate = 1:5,
  RAPID = sapply(rapid_at_default, function(x) round(x$rapid, 4)),
  CI_Lower = sapply(rapid_at_default, function(x) round(x$ci[1], 4)),
  CI_Upper = sapply(rapid_at_default, function(x) round(x$ci[2], 4)),
  N_at_risk = sapply(rapid_at_default, function(x) x$n_at_risk)
)

cat("\nRAPID estimates at tau = 0.3:\n")
print(summary_table)

cat("\nMean RAPID across replicates:",
    round(mean(summary_table$RAPID), 4), "\n")

################################################################################
# 7. Save results
################################################################################

save(rapid_df, summary_table, rapid_matrix,
     file = "/Users/matthias/workspace26/RAPID/dev/real_data_results_new_rapid.RData")

cat("\nResults saved to: dev/real_data_results_new_rapid.RData\n")
cat("Done!\n")
