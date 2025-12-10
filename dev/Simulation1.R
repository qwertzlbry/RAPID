################################################################################
## SIMULATION 1: DEPENDENCY STRENGTH (ORIGINAL VERSION)
## κ = 0,::100
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

run_simulation <- function(kappa, n = 1000, synthesizer = "cart", attacker = "rf",
                           tau = 0.3, n_sim = 10, seed_base = 2025) {
  results <- list()

  for (sim in 1:n_sim) {
    cat(sprintf("κ = %d | Sim %d/%d\n", kappa, sim, n_sim))

    orig_data <- simulate_microdata(n = n, dep = kappa, seed = seed_base + sim)

    tryCatch({
      syn_obj <- syn(orig_data, method = synthesizer, print.flag = FALSE, seed = seed_base + sim + 1000)
      syn_data <- syn_obj$syn

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

      results[[length(results) + 1]] <- data.frame(
        kappa = kappa, sim = sim,
        RAPID = rapid_result$risk$confidence_rate,
        n_at_risk = rapid_result$risk$n_at_risk,
        accuracy = rapid_result$metrics$accuracy
      )
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
      results[[length(results) + 1]] <- data.frame(
        kappa = kappa, sim = sim, RAPID = NA, n_at_risk = NA, accuracy = NA
      )
    })
  }

  do.call(rbind, results)
}

################################################################################
## RUN SIMULATION
################################################################################

kappa_values <- seq(0, 100, by = 1)
n_sim <- 10

sim1_results <- do.call(rbind, lapply(kappa_values, function(k) {
  run_simulation(kappa = k, n = 1000, synthesizer = "cart", attacker = "rf",
                 tau = 0.3, n_sim = n_sim)
}))

save(sim1_results, file = "simulation1_dependency_strength.RData")

################################################################################
## ANALYZE & VISUALIZE
################################################################################

sim1_summary <- sim1_results %>%
  group_by(kappa) %>%
  summarise(
    RAPID_mean = mean(RAPID, na.rm = TRUE),
    RAPID_sd = sd(RAPID, na.rm = TRUE),
    accuracy_mean = mean(accuracy, na.rm = TRUE),
    accuracy_sd = sd(accuracy, na.rm = TRUE),
    .groups = "drop"
  )

print(head(sim1_summary, 20))
print(tail(sim1_summary, 20))

# Plot: Combined (like your old one)
sim1_long <- sim1_summary %>%
  select(kappa, RAPID_mean, RAPID_sd, accuracy_mean, accuracy_sd) %>%
  pivot_longer(
    cols = c(RAPID_mean, accuracy_mean),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    sd_value = ifelse(metric == "RAPID_mean", RAPID_sd, accuracy_sd),
    metric = factor(metric,
                    levels = c("RAPID_mean", "accuracy_mean"),
                    labels = c("Confidence Rate", "Accuracy"))
  )

p <- ggplot(sim1_long, aes(x = kappa, y = value, color = metric, group = metric)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = value - sd_value, ymax = value + sd_value, fill = metric),
              alpha = 0.2, color = NA) +
  scale_color_manual(values = c("Confidence Rate" = "steelblue",
                                "Accuracy" = "firebrick")) +
  scale_fill_manual(values = c("Confidence Rate" = "steelblue",
                               "Accuracy" = "firebrick")) +
  labs(
   # title = "Attribute Disclosure Risk and Accuracy vs Dependency Strength",
    # subtitle = sprintf("Mean ± 1 SD over repeated simulations (n = %d per dep)", n_sim),
    x = "κ (dependency strength)",
    y = "RAPID (Proportion at Risk)"
  ) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave("sim1_combined.pdf", p, width = 10, height = 6)
print(p)

# Table
table_sim1 <- sim1_summary %>%
  filter(kappa %in% c(0, 25, 50, 75, 100)) %>%
  mutate(
    Kappa = kappa,
    RAPID = sprintf("%.3f (%.3f)", RAPID_mean, RAPID_sd),
    Accuracy = sprintf("%.3f (%.3f)", accuracy_mean, accuracy_sd)
  ) %>%
  select(Kappa, RAPID, Accuracy)

print(table_sim1)
write.csv(table_sim1, "sim1_table.csv", row.names = FALSE)


