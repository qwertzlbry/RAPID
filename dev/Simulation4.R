################################################################################
## SIMULATION 4: ATTACKER COMPARISON
## Research Question: Which attacker model is most effective?
## Fixed: κ = 10, τ = 0.3, CART synthesizer
## Varied: Attacker models (RF, CART, GBM, LOGIT)
################################################################################

library(devtools)
load_all()
library(synthpop)
library(dplyr)
library(ggplot2)
library(tidyr)

set.seed(2025)

################################################################################
## DATA GENERATION (Same as before)
################################################################################

simulate_microdata <- function(n = 1000, seed = NULL, dep = 1,
                               age_mean = 45, age_sd = 15,
                               age_min = 18, age_max = 90) {
  if (!is.null(seed)) set.seed(seed)

  signal_weight <- sqrt(dep / (1 + dep))
  noise_weight  <- sqrt(1 / (1 + dep))

  rtruncnorm <- function(n, mean, sd, lo, hi) {
    x <- rnorm(n, mean, sd)
    while (any(out <- x < lo | x > hi)) {
      x[out] <- rnorm(sum(out), mean, sd)
    }
    x
  }

  age <- rtruncnorm(n, mean = age_mean, sd = age_sd, lo = age_min, hi = age_max)
  age_z <- scale(age)[, 1]
  latent_SES <- rnorm(n)

  edu_latent <- signal_weight * (0.8 * latent_SES - 0.4 * age_z) + noise_weight * rnorm(n)
  education_num <- cut(edu_latent, breaks = c(-Inf, -0.3, 0.7, Inf), labels = c(0,1,2)) |> as.integer() - 1
  education <- factor(c("low", "medium", "high")[education_num + 1], levels = c("low", "medium", "high"))

  lin_income <- signal_weight * (0.5 * latent_SES + 0.3 * age_z + 0.25 * education_num) + noise_weight * rnorm(n)
  income <- exp(10 + lin_income)

  lin_health <- signal_weight * (0.6 * latent_SES - 0.5 * age_z + 0.2 * education_num + 0.2 * scale(log(income))[,1]) + noise_weight * rnorm(n)
  health_score <- 100 / (1 + exp(-lin_health))

  health_z <- scale(health_score)[, 1]
  log_income_z <- scale(log(income))[, 1]

  eta_diab <- -1.5 + dep * (0.8 * age_z - 0.3 * log_income_z - 0.2 * education_num)
  eta_hyp <- -1.3 + dep * (1.0 * age_z - 0.2 * log_income_z - 0.1 * education_num)
  eta_healthy <- 0

  exp_eta <- cbind(healthy = exp(eta_healthy), diabetic = exp(eta_diab), hypertensive = exp(eta_hyp))
  prob_mat <- exp_eta / rowSums(exp_eta)

  draw_multinom <- function(p, levels) {
    idx <- apply(p, 1, function(prob) sample(seq_along(prob), size = 1, prob = prob))
    factor(levels[idx], levels = levels)
  }
  disease_status <- draw_multinom(prob_mat, c("healthy", "diabetic", "hypertensive"))

  eta_gender <- signal_weight * (0.3 * latent_SES - 0.2 * age_z + 0.2 * education_num)
  prob_male <- 1 / (1 + exp(-eta_gender))
  gender <- factor(rbinom(n, 1, prob_male), labels = c("female", "male"))

  data.frame(gender = gender, age = age, education = education,
             income = income, health_score = health_score, disease_status = disease_status)
}

################################################################################
## SIMULATION FUNCTION
################################################################################

run_simulation_attacker <- function(attacker, kappa = 10, n = 1000, tau = 0.3,
                                    synthesizer = "cart", sim = 1, seed_base = 6025) {

  cat(sprintf("Attacker: %s | Sim %d\n", attacker, sim))

  # Generate data
  orig_data <- simulate_microdata(n = n, dep = kappa, seed = seed_base + sim)

  tryCatch({
    syn_obj <- syn(orig_data, method = synthesizer, print.flag = FALSE,
                   seed = seed_base + sim + 1000)
    syn_data <- syn_obj$syn

    # Compute RAPID
    rapid_result <- rapid(
      original_data = orig_data,
      synthetic_data = syn_data,
      quasi_identifiers = c("gender", "age", "education", "income", "health_score"),
      sensitive_attribute = "disease_status",
      model_type = attacker,
      cat_eval_method = "RCS_marginal",
      cat_tau = tau,
      seed = seed_base + sim + 2000,
      trace = FALSE
    )

    data.frame(
      attacker = attacker,
      sim = sim,
      RAPID = rapid_result$risk$confidence_rate,
      n_at_risk = rapid_result$risk$n_at_risk,
      accuracy = rapid_result$metrics$accuracy
    )

  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    data.frame(
      attacker = attacker, sim = sim,
      RAPID = NA, n_at_risk = NA, accuracy = NA
    )
  })
}

################################################################################
## RUN SIMULATION 4
################################################################################

cat("\n=== RUNNING SIMULATION 4: ATTACKER COMPARISON ===\n")

attackers <- c("rf", "cart", "gbm")
kappa_fixed <- 10
tau_fixed <- 0.3
n_sim <- 50

sim4_results <- do.call(rbind, lapply(attackers, function(att) {
  do.call(rbind, lapply(1:n_sim, function(s) {
    run_simulation_attacker(
      attacker = att,
      kappa = kappa_fixed,
      tau = tau_fixed,
      sim = s
    )
  }))
}))

save(sim4_results, file = "simulation4_attacker_comparison.RData")
cat("\n=== Results Saved ===\n")

################################################################################
## ANALYZE
################################################################################

# Summary statistics
sim4_summary <- sim4_results %>%
  group_by(attacker) %>%
  summarise(
    RAPID_mean = mean(RAPID, na.rm = TRUE),
    RAPID_sd = sd(RAPID, na.rm = TRUE),
    RAPID_se = sd(RAPID, na.rm = TRUE) / sqrt(n()),
    accuracy_mean = mean(accuracy, na.rm = TRUE),
    accuracy_sd = sd(accuracy, na.rm = TRUE),
    n_obs = sum(!is.na(RAPID)),
    .groups = "drop"
  ) %>%
  arrange(desc(RAPID_mean))

cat("\n=== ATTACKER COMPARISON SUMMARY ===\n")
print(sim4_summary)

# Format for nice labels
sim4_results <- sim4_results %>%
  mutate(attacker_label = factor(
    attacker,
    levels = c("gbm", "rf", "cart", "logit"),
    labels = c("GBM", "Random Forest", "CART", "Logistic Regression")
  ))

sim4_summary <- sim4_summary %>%
  mutate(attacker_label = factor(
    attacker,
    levels = c("gbm", "rf", "cart", "logit"),
    labels = c("GBM", "Random Forest", "CART", "Logistic Regression")
  ))

################################################################################
## VISUALIZE
################################################################################

# Plot 1: RAPID comparison (violin + boxplot)
p1 <- ggplot(sim4_results, aes(x = reorder(attacker_label, -RAPID), y = RAPID, fill = attacker_label)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.alpha = 0.3, notch = TRUE) +
  scale_fill_manual(
    values = c("GBM" = "darkred",
               "Random Forest" = "steelblue",
               "CART" = "darkorange",
               "Logistic Regression" = "gray50")
  ) +
  labs(
    title = "RAPID by Attacker Model",
    subtitle = sprintf("κ = %d, τ = %.1f | %d simulations per attacker",
                       kappa_fixed, tau_fixed, n_sim),
    x = "Attacker Model",
    y = "RAPID (Proportion at Risk)",
    caption = "Higher RAPID = More dangerous attacker"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

ggsave("sim4_rapid_by_attacker.pdf", p1, width = 8, height = 6)
print(p1)

# Plot 2: Accuracy comparison
p2 <- ggplot(sim4_results, aes(x = reorder(attacker_label, -accuracy), y = accuracy, fill = attacker_label)) +
  geom_violin(alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.alpha = 0.3, notch = TRUE) +
  scale_fill_manual(
    values = c("GBM" = "darkred",
               "Random Forest" = "steelblue",
               "CART" = "darkorange",
               "Logistic Regression" = "gray50")
  ) +
  labs(
    title = "Attacker Model Accuracy",
    subtitle = sprintf("κ = %d | %d simulations per attacker", kappa_fixed, n_sim),
    x = "Attacker Model",
    y = "Classification Accuracy",
    caption = "Higher accuracy enables higher confidence predictions"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

ggsave("sim4_accuracy_by_attacker.pdf", p2, width = 8, height = 6)
print(p2)

# Plot 3: RAPID vs Accuracy scatter
p3 <- ggplot(sim4_results, aes(x = accuracy, y = RAPID, color = attacker_label)) +
  geom_point(alpha = 0.4, size = 2) +
  stat_ellipse(level = 0.95, size = 1) +
  scale_color_manual(
    values = c("GBM" = "darkred",
               "Random Forest" = "steelblue",
               "CART" = "darkorange",
               "Logistic Regression" = "gray50")
  ) +
  labs(
    title = "RAPID vs Accuracy by Attacker Model",
    subtitle = "95% confidence ellipses | Higher accuracy generally enables higher RAPID",
    x = "Classification Accuracy",
    y = "RAPID (Proportion at Risk)",
    color = "Attacker Model"
  ) +
  scale_x_continuous(limits = c(0.7, 1), breaks = seq(0.7, 1, 0.05)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("sim4_rapid_vs_accuracy.pdf", p3, width = 8, height = 6)
print(p3)

# Plot 4: Bar chart with error bars (clean summary)
p4 <- ggplot(sim4_summary, aes(x = reorder(attacker_label, -RAPID_mean), y = RAPID_mean, fill = attacker_label)) +
  geom_col(alpha = 0.7, width = 0.7) +
  geom_errorbar(aes(ymin = RAPID_mean - RAPID_sd, ymax = RAPID_mean + RAPID_sd),
                width = 0.2, size = 1) +
  geom_text(aes(label = sprintf("%.3f", RAPID_mean)),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(
    values = c("GBM" = "darkred",
               "Random Forest" = "steelblue",
               "CART" = "darkorange",
               "Logistic Regression" = "gray50")
  ) +
  labs(
    title = "Mean RAPID by Attacker Model",
    subtitle = sprintf("κ = %d, τ = %.1f | Error bars: ±1 SD", kappa_fixed, tau_fixed),
    x = "Attacker Model",
    y = "RAPID (Mean over 50 simulations)"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

ggsave("sim4_mean_rapid_bar.pdf", p4, width = 8, height = 6)
print(p4)

################################################################################
## TABLE FOR PAPER
################################################################################

table_sim4 <- sim4_summary %>%
  mutate(
    `RAPID (SD)` = sprintf("%.3f (%.3f)", RAPID_mean, RAPID_sd),
    `Accuracy (SD)` = sprintf("%.3f (%.3f)", accuracy_mean, accuracy_sd),
    Rank = row_number()
  ) %>%
  select(Rank, attacker_label, `RAPID (SD)`, `Accuracy (SD)`) %>%
  rename(Attacker = attacker_label)

cat("\n=== TABLE FOR PAPER ===\n")
print(table_sim4)
write.csv(table_sim4, "sim4_table.csv", row.names = FALSE)

################################################################################
## STATISTICAL TESTS
################################################################################

# Pairwise comparisons (Wilcoxon rank-sum tests)
cat("\n=== PAIRWISE COMPARISONS (Wilcoxon tests) ===\n")

attacker_pairs <- combn(unique(sim4_results$attacker), 2, simplify = FALSE)

pairwise_tests <- do.call(rbind, lapply(attacker_pairs, function(pair) {
  data1 <- sim4_results %>% filter(attacker == pair[1]) %>% pull(RAPID)
  data2 <- sim4_results %>% filter(attacker == pair[2]) %>% pull(RAPID)

  test <- wilcox.test(data1, data2)

  data.frame(
    comparison = paste(pair[1], "vs", pair[2]),
    p_value = test$p.value,
    significant = ifelse(test$p.value < 0.05, "Yes", "No")
  )
}))

print(pairwise_tests)
write.csv(pairwise_tests, "sim4_pairwise_tests.csv", row.names = FALSE)

################################################################################
## SUMMARY
################################################################################

cat("\n=== SIMULATION 4 COMPLETE ===\n")
cat("Files saved:\n")
cat("  - simulation4_attacker_comparison.RData\n")
cat("  - sim4_rapid_by_attacker.pdf (MAIN RESULT)\n")
cat("  - sim4_accuracy_by_attacker.pdf\n")
cat("  - sim4_rapid_vs_accuracy.pdf\n")
cat("  - sim4_mean_rapid_bar.pdf\n")
cat("  - sim4_table.csv\n")
cat("  - sim4_pairwise_tests.csv\n")

cat("\n=== ATTACKER RANKING (by mean RAPID) ===\n")
print(sim4_summary %>% select(attacker_label, RAPID_mean, accuracy_mean))
