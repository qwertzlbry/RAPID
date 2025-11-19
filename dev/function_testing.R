library(tidyverse)
library(progressr)
library(synthpop)

# edge cases -------------------------------------------------------------------
ng <- function(gi, bi){
  out <- (gi-bi)/(1-bi)
  round(out,3)
}

# Create sequences
gi_vals <- seq(0, 1, by = 0.1)
bi_vals <- seq(0.1, 1, by = 0.1)

# Expand grid of all combinations
grid <- expand.grid(gi = gi_vals, bi = bi_vals)

grid$ng <- mapply(ng, grid$gi, grid$bi)

# Make matrix (rows = gi, columns = bi)
ng_matrix <- with(grid, tapply(ng, list(gi, bi), mean))

grid$at_risk <- grid$ng >= 0.3


ggplot(grid, aes(x = bi, y = gi, fill = ng)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(ng, 2)), size = 2.5) +
  scale_fill_gradientn(
    colors = c("darkgreen", "green", "yellow", "orange", "red", "darkred"),
    values = scales::rescale(c(-9, -2, 0, 0.3, 0.6, 1)),
    limits = c(-9, 1)
  ) +
  geom_contour(aes(z = ng), breaks = 0.3, color = "black", linewidth = 0.5) +
  labs(
    title = "Normalized Gain: Edge Case Analysis",
    subtitle = expression(
      "normalized_gain = " * frac(g[i] - b[i], 1 - b[i]) *
        " | Black line = " * tau * " = 0.3 | Green = safe, Red = at risk"
    ),
    x = expression("Baseline (" * b[i] * ")"),
    y = expression("Predicted Probability (" * g[i] * ")"),
    fill = "Normalized\nGain"
  ) +
  theme_minimal(base_size = 14) +
  coord_fixed()

# tau sensitiviy --------------------------------------------------------------

test_tau_sensitivity <- function() {

  tau_values <- seq(0.1, 0.7, by = 0.1)
  dep_values <- seq(0, 100, by = 10)
  n_rep <- 3  # Wiederholungen pro Kombination

  results <- with_progress({
    p <- progressor(steps = length(tau_values) * length(dep_values) * n_rep)

    expand.grid(
      tau = tau_values,
      dep = dep_values,
      rep = 1:n_rep
    ) %>%
      rowwise() %>%
      mutate(
        data = {
          p(sprintf("tau=%.2f, dep=%d, rep=%d", tau, dep, rep))

          # Deine echte Simulation!
          truth <- simulate_microdata(dep = dep, n = 1000, seed = rep * 1000 + dep)
          synth <- synthpop::syn(truth, m = 1, print.flag = FALSE, seed = rep * 1000 + dep + 1)$syn

          rapid_result <- rapid(
            original_data = truth,
            synthetic_data = synth,
            quasi_identifiers = c("age", "income", "education"),
            sensitive_attribute = "disease_status",
            model_type = "rf",
            tau = tau,
            cat_eval_method = "RCS_marginal",
            trace = FALSE
          )

          tibble(
            confidence_rate = rapid_result$risk$confidence_rate,
            accuracy = rapid_result$metrics$accuracy,
            mean_gain = mean(rapid_result$normalized_gain, na.rm = TRUE),
            sd_gain = sd(rapid_result$normalized_gain, na.rm = TRUE)
          )
        }
      ) %>%
      unnest(data)
  })

  return(results)
}

# Führe aus
tau_dep_results <- test_tau_sensitivity()

# ===== Visualisierung =====
library(ggplot2)

# Plot 1: Confidence Rate über dependency, farbig nach tau
ggplot(tau_dep_results, aes(x = dep, y = confidence_rate, color = factor(tau))) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 2, alpha = 0.5) +
  scale_color_viridis_d(option = "plasma") +
  labs(
    title = "Tau Sensitivity: Confidence Rate vs Dependency Strength",
    subtitle = "How does threshold τ affect risk measurement across data qualities?",
    x = "Dependency Strength (κ)",
    y = "Confidence Rate (proportion at risk)",
    color = "Threshold τ"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

# Plot 2: Heatmap
summary_results <- tau_dep_results %>%
  group_by(tau, dep) %>%
  summarise(mean_conf = mean(confidence_rate), .groups = "drop")

ggplot(summary_results, aes(x = dep, y = tau, fill = mean_conf)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", mean_conf)), size = 3) +
  scale_fill_gradient2(
    low = "green", mid = "yellow", high = "red",
    midpoint = 0.5,
    limits = c(0, 1)
  ) +
  labs(
    title = "Sensitivity Heatmap: τ × Dependency",
    x = "Dependency Strength (κ)",
    y = "Threshold τ",
    fill = "Confidence\nRate"
  ) +
  theme_minimal(base_size = 14)

# Plot 3: Verschiedene Metriken ------------------------------------------------
tau_dep_results %>%
  pivot_longer(cols = c(confidence_rate, accuracy, mean_gain),
               names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = dep, y = value, color = factor(tau))) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_color_viridis_d(option = "plasma") +
  labs(
    title = "Multiple Metrics: Tau Sensitivity Analysis",
    x = "Dependency Strength (κ)",
    y = "Metric Value",
    color = "Threshold τ"
  ) +
  theme_minimal(base_size = 14)

# test number of classes

# ===== Test: Verschiedene Anzahl Klassen =====
test_n_classes <- function() {

  n_classes_values <- c(2, 3, 5, 10, 20, 50)
  tau <- 0.3

  results <- map_dfr(n_classes_values, function(K) {

    map_dfr(1:10, function(rep) {  # 10 Wiederholungen

      n <- 1000

      # Generiere K Klassen mit unterschiedlichen Frequenzen
      # Zipf-like distribution: häufige und seltene Klassen
      freqs <- (1:K)^(-1.5)
      freqs <- freqs / sum(freqs)

      # Original data
      true_class <- sample(1:K, n, replace = TRUE, prob = freqs)

      # Predicted probabilities (simuliert)
      predicted_probs <- matrix(runif(n * K), ncol = K)
      predicted_probs <- predicted_probs / rowSums(predicted_probs)

      # Normalized gain berechnen
      g_i <- predicted_probs[cbind(1:n, true_class)]

      marginal_freq <- prop.table(table(factor(true_class, levels = 1:K)))
      b_i <- marginal_freq[true_class]

      normalized_gain <- (g_i - b_i) / (1 - b_i)

      # Metriken
      tibble(
        n_classes = K,
        rep = rep,
        confidence_rate = mean(normalized_gain > tau, na.rm = TRUE),
        mean_gain = mean(normalized_gain, na.rm = TRUE),
        sd_gain = sd(normalized_gain, na.rm = TRUE),
        min_baseline = min(marginal_freq),
        max_baseline = max(marginal_freq),
        entropy = -sum(marginal_freq * log(marginal_freq + 1e-10))
      )
    })
  })

  return(results)
}

# Führe Test aus
class_results <- test_n_classes()

# ===== Visualisierung: Multiple Plots =====
library(patchwork)

# Plot 1: Confidence Rate vs. Anzahl Klassen
p1 <- ggplot(class_results, aes(x = n_classes, y = confidence_rate)) +
  geom_jitter(alpha = 0.3, width = 0.1) +
  stat_summary(fun = mean, geom = "line", color = "red", size = 1.2) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 3) +
  labs(
    title = "Confidence Rate vs. Number of Classes",
    x = "Number of Classes (K)",
    y = "Confidence Rate"
  ) +
  theme_minimal()

# Plot 2: Mean Gain vs. Anzahl Klassen
p2 <- ggplot(class_results, aes(x = n_classes, y = mean_gain)) +
  geom_jitter(alpha = 0.3, width = 0.1) +
  stat_summary(fun = mean, geom = "line", color = "blue", size = 1.2) +
  stat_summary(fun = mean, geom = "point", color = "blue", size = 3) +
  labs(
    title = "Mean Normalized Gain vs. Number of Classes",
    x = "Number of Classes (K)",
    y = "Mean Normalized Gain"
  ) +
  theme_minimal()

# Plot 3: Baseline Range vs. Anzahl Klassen
p3 <- class_results %>%
  select(n_classes, rep, min_baseline, max_baseline) %>%
  pivot_longer(cols = c(min_baseline, max_baseline), names_to = "type", values_to = "baseline") %>%
  ggplot(aes(x = n_classes, y = baseline, color = type)) +
  geom_jitter(alpha = 0.3, width = 0.1) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  scale_color_manual(
    values = c("min_baseline" = "purple", "max_baseline" = "orange"),
    labels = c("Minimum", "Maximum")
  ) +
  labs(
    title = "Baseline Range vs. Number of Classes",
    x = "Number of Classes (K)",
    y = "Baseline Frequency",
    color = "Baseline"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Kombiniere Plots
(p1 | p2) / p3 +
  plot_annotation(
    title = "Scalability Analysis: How does the metric behave with many classes?",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )
