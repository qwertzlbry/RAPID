################################################################################
## SIMULATION 2: THRESHOLD SENSITIVITY (MULTIPLE KAPPA)
## Research Question: How does RAPID vary with τ across different κ levels?
## Varied: τ = 0.05, 0.10, ..., 0.95 AND κ = 0, 5, 10, 20, 50
################################################################################

library(devtools)
load_all()
library(synthpop)
library(dplyr)
library(ggplot2)
library(tidyr)
source("dev/simulate_microdata.R")
set.seed(2025)


################################################################################
## SIMULATION FUNCTION
################################################################################

run_simulation_tau_kappa <- function(tau, kappa, n = 1000,
                                     synthesizer = "cart", attacker = "rf",
                                     n_sim = 10, seed_base = 3025) {
  results <- list()

  for (sim in 1:n_sim) {
    cat(sprintf("κ = %d, τ = %.2f | Sim %d/%d\n", kappa, tau, sim, n_sim))

    orig_data <- simulate_microdata(n = n, dep = kappa, seed = seed_base + kappa*1000 + sim)

    tryCatch({
      syn_obj <- syn(orig_data, method = synthesizer, print.flag = FALSE,
                     seed = seed_base + kappa*1000 + sim + 10000)
      syn_data <- syn_obj$syn

      rapid_result <- rapid(
        original_data = orig_data,
        synthetic_data = syn_data,
        quasi_identifiers = c("gender", "age", "education", "income", "health_score"),
        sensitive_attribute = "disease_status",
        model_type = attacker,
        cat_eval_method = "RCS_marginal",
        cat_tau = tau,
        seed = seed_base + kappa*1000 + sim + 20000,
        trace = FALSE
      )

      results[[length(results) + 1]] <- data.frame(
        kappa = kappa, tau = tau, sim = sim,
        RAPID = rapid_result$risk$confidence_rate,
        n_at_risk = rapid_result$risk$n_at_risk,
        accuracy = rapid_result$metrics$accuracy
      )
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      results[[length(results) + 1]] <- data.frame(
        kappa = kappa, tau = tau, sim = sim,
        RAPID = NA, n_at_risk = NA, accuracy = NA
      )
    })
  }

  do.call(rbind, results)
}

################################################################################
## RUN SIMULATION 2
################################################################################

tau_values <- seq(0.05, 0.95, by = 0.05)  # 19 points
kappa_values <- c(0, 5, 10, 20, 50)       # 5 dependency levels
n_sim <- 10

# Run for all combinations of tau and kappa
sim2_results <- do.call(rbind, lapply(kappa_values, function(k) {
  do.call(rbind, lapply(tau_values, function(t) {
    run_simulation_tau_kappa(
      tau = t,
      kappa = k,
      n = 1000,
      synthesizer = "cart",
      attacker = "rf",
      n_sim = n_sim
    )
  }))
}))

save(sim2_results, file = "simulation2_threshold_sensitivity.RData")

################################################################################
## ANALYZE & VISUALIZE
################################################################################

sim2_summary <- sim2_results %>%
  group_by(kappa, tau) %>%
  summarise(
    RAPID_mean = mean(RAPID, na.rm = TRUE),
    RAPID_sd = sd(RAPID, na.rm = TRUE),
    RAPID_se = sd(RAPID, na.rm = TRUE) / sqrt(n()),
    n_at_risk_mean = mean(n_at_risk, na.rm = TRUE),
    accuracy_mean = mean(accuracy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(kappa_label = factor(paste0("κ = ", kappa),
                              levels = paste0("κ = ", c(0, 5, 10, 20, 50))))

print(sim2_summary)

# Plot 1: Multiple lines (main result)
p1 <- ggplot(sim2_summary, aes(x = tau, y = RAPID_mean, color = kappa_label, group = kappa_label)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c("κ = 0" = "gray50",
               "κ = 5" = "steelblue",
               "κ = 10" = "dodgerblue",
               "κ = 20" = "firebrick",
               "κ = 50" = "darkred")
  ) +
  labs(
    #title = "RAPID vs Confidence Threshold for Different Dependency Strengths",
    #subtitle = sprintf("Mean over %d simulations per (κ, τ) combination", n_sim),
    x = "Confidence Threshold (τ)",
    y = "RAPID (Proportion at Risk)",
    color = "Dependency\nStrength",
    #caption = "Synthesizer: CART | Attacker: Random Forest"
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("sim2_rapid_vs_tau_multiple_kappa.pdf", p1, width = 10, height = 6)
print(p1)

# Plot 2: With ribbons for one kappa (e.g., κ=10)
sim2_kappa10 <- sim2_summary %>% filter(kappa == 10)

p2 <- ggplot(sim2_kappa10, aes(x = tau, y = RAPID_mean)) +
  geom_line(size = 1.2, color = "dodgerblue") +
  geom_point(size = 2, color = "dodgerblue") +
  geom_ribbon(aes(ymin = RAPID_mean - RAPID_sd, ymax = RAPID_mean + RAPID_sd),
              alpha = 0.2, fill = "dodgerblue") +
  labs(
    title = "RAPID vs Confidence Threshold",
    subtitle = "Fixed dependency strength κ = 10 | Mean ± 1 SD",
    x = "Confidence Threshold (τ)",
    y = "RAPID (Proportion at Risk)"
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("sim2_rapid_vs_tau_kappa10.pdf", p2, width = 8, height = 6)
print(p2)

# Plot 3: Faceted by kappa
p3 <- ggplot(sim2_summary, aes(x = tau, y = RAPID_mean)) +
  geom_line(size = 1, color = "steelblue") +
  geom_point(size = 1.5, color = "steelblue") +
  facet_wrap(~kappa_label, ncol = 5) +
  labs(
    title = "RAPID vs Confidence Threshold by Dependency Strength",
    x = "Confidence Threshold (τ)",
    y = "RAPID"
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave("sim2_faceted_by_kappa.pdf", p3, width = 12, height = 4)
print(p3)

# Table for paper
table_sim2 <- sim2_summary %>%
  filter(tau %in% c(0.2, 0.3, 0.4, 0.5)) %>%
  select(kappa, tau, RAPID_mean, RAPID_sd) %>%
  mutate(
    RAPID = sprintf("%.3f (%.3f)", RAPID_mean, RAPID_sd)
  ) %>%
  select(kappa, tau, RAPID) %>%
  pivot_wider(names_from = tau, values_from = RAPID, names_prefix = "τ = ")

print(table_sim2)
write.csv(table_sim2, "sim2_table.csv", row.names = FALSE)

