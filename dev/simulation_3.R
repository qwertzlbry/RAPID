################################################################################
## SIMULATION 3: QI ATTRIBUTION ANALYSIS
## Method: Logistic regression with full interaction model
## QIs: gender, age, education (demographic identifiers only)
################################################################################
library(devtools)
load_all()
library(synthpop)
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)
source("dev/simulate_microdata.R")
set.seed(2025)

################################################################################
## REGRESSION ANALYSIS FUNCTION - FULL MODEL
################################################################################
run_regression_analysis <- function(kappa = 10, n = 1000, tau = 0.3,
                                    synthesizer = "cart", attacker = "rf",
                                    sim = 1, seed_base = 2025) {
  # Generate data
  orig_data <- simulate_microdata(n = n, dep = kappa, seed = seed_base + sim)
  syn_obj <- syn(orig_data, method = synthesizer, print.flag = FALSE,
                 seed = seed_base + sim + 1000)
  syn_data <- syn_obj$syn

  # Run RAPID with only demographic QIs
  rapid_result <- rapid(
    original_data = orig_data,
    synthetic_data = syn_data,
    quasi_identifiers = c("gender", "age", "education"),
    sensitive_attribute = "disease_status",
    model_type = attacker,
    cat_eval_method = "RCS_marginal",
    cat_tau = tau,
    return_all_records = TRUE,
    seed = seed_base + sim + 2000,
    trace = FALSE
  )

  # Extract records
  records <- rapid_result$risk$rows_risk_df

  # Set reference levels
  records <- records %>%
    mutate(
      education = relevel(education, ref = "low"),
      gender = relevel(gender, ref = "female")
    )

  # FIT FULL LOGISTIC REGRESSION with all interactions
  logit_model <- glm(
    at_risk ~ gender * age * education,
    data = records,
    family = binomial()
  )

  # Extract coefficients
  coef_results <- tidy(logit_model) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      # Classify effect type
      variable = case_when(
        # 3-way interaction
        grepl("gender.*:.*age.*:.*education|age.*:.*gender.*:.*education", term) ~
          "Gender × Age × Education",
        # 2-way interactions
        grepl("^gendermale:age$", term) ~ "Gender × Age",
        grepl("^gendermale:education", term) ~ "Gender × Education",
        grepl("^age:education", term) ~ "Age × Education",
        # Main effects
        grepl("^gendermale$", term) ~ "Gender",
        term == "age" ~ "Age",
        grepl("^education", term) ~ "Education",
        TRUE ~ "Other"
      ),
      category = case_when(
        # Main effects
        term == "gendermale" ~ "Male",
        term == "age" ~ "Per year",
        term == "educationmedium" ~ "Medium",
        term == "educationhigh" ~ "High",
        # 2-way: Gender × Age
        term == "gendermale:age" ~ "Male × Age",
        # 2-way: Gender × Education
        term == "gendermale:educationmedium" ~ "Male × Medium",
        term == "gendermale:educationhigh" ~ "Male × High",
        # 2-way: Age × Education
        term == "age:educationmedium" ~ "Age × Medium",
        term == "age:educationhigh" ~ "Age × High",
        # 3-way interactions
        term == "gendermale:age:educationmedium" ~ "Male × Age × Medium",
        term == "gendermale:age:educationhigh" ~ "Male × Age × High",
        TRUE ~ term
      ),
      is_reference = FALSE,
      sim = sim
    ) %>%
    select(sim, variable, category, estimate, std.error, p.value, is_reference)

  # ADD REFERENCE CATEGORIES (β = 0)
  reference_cats <- data.frame(
    sim = sim,
    variable = c("Gender", "Education",
                 "Gender × Age", "Gender × Education", "Age × Education",
                 "Gender × Age × Education"),
    category = c("Female", "Low",
                 "Female × Age", "Female × Low", "Age × Low",
                 "Female × Age × Low"),
    estimate = 0,
    std.error = 0,
    p.value = NA,
    is_reference = TRUE
  )

  # Combine
  all_results <- bind_rows(coef_results, reference_cats)

  return(all_results)
}

################################################################################
## RUN SIMULATION 3
################################################################################
kappa_fixed <- 10
tau_fixed <- 0.3
n_sim <- 50

sim3_results <- do.call(rbind, lapply(1:n_sim, function(s) {
  cat("Running simulation", s, "of", n_sim, "\n")
  run_regression_analysis(
    kappa = kappa_fixed,
    tau = tau_fixed,
    sim = s
  )
}))

save(sim3_results, file = "simulation3_qi_attribution.RData")

################################################################################
## ANALYZE
################################################################################
coef_summary <- sim3_results %>%
  group_by(variable, category, is_reference) %>%
  summarise(
    beta_mean = mean(estimate),
    beta_sd = sd(estimate),
    beta_lower = quantile(estimate, 0.025),
    beta_upper = quantile(estimate, 0.975),
    sig_pct = ifelse(is_reference[1], NA, mean(p.value < 0.05) * 100),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  arrange(variable, is_reference, desc(abs(beta_mean)))

# Variable importance by coefficient magnitude
var_importance <- coef_summary %>%
  filter(!is_reference) %>%
  group_by(variable) %>%
  summarise(
    max_abs_beta = max(abs(beta_mean)),
    mean_abs_beta = mean(abs(beta_mean)),
    strongest_category = category[which.max(abs(beta_mean))],
    .groups = "drop"
  ) %>%
  arrange(desc(max_abs_beta))

print("Variable Importance:")
print(var_importance)

################################################################################
## PREPARE PLOT DATA
################################################################################

# PLOT 1 DATA: Main Effects
main_effects_data <- sim3_results %>%
  filter(variable %in% c("Gender", "Age", "Education")) %>%
  mutate(
    variable = factor(variable, levels = c("Gender", "Age", "Education")),
    category = case_when(
      variable == "Gender" ~ factor(category, levels = c("Female", "Male")),
      variable == "Age" ~ factor(category, levels = c("Per year")),
      variable == "Education" ~ factor(category, levels = c("Low", "Medium", "High")),
      TRUE ~ factor(category)
    )
  )

main_effects_summary <- coef_summary %>%
  filter(variable %in% c("Gender", "Age", "Education")) %>%
  mutate(
    sig_stars = case_when(
      is_reference ~ "",
      sig_pct >= 95 ~ "***",
      sig_pct >= 90 ~ "**",
      sig_pct >= 80 ~ "*",
      TRUE ~ ""
    ),
    variable = factor(variable, levels = c("Gender", "Age", "Education")),
    category = case_when(
      variable == "Gender" ~ factor(category, levels = c("Female", "Male")),
      variable == "Age" ~ factor(category, levels = c("Per year")),
      variable == "Education" ~ factor(category, levels = c("Low", "Medium", "High")),
      TRUE ~ factor(category)
    )
  )

# PLOT 2 DATA: All Interactions
interaction_data <- sim3_results %>%
  filter(variable %in% c("Gender × Age", "Gender × Education",
                         "Age × Education", "Gender × Age × Education")) %>%
  mutate(
    variable = factor(variable,
                      levels = c("Gender × Age", "Gender × Education",
                                 "Age × Education", "Gender × Age × Education")),
    # Better category ordering
    category = case_when(
      variable == "Gender × Age" ~
        factor(category, levels = c("Female × Age", "Male × Age")),
      variable == "Gender × Education" ~
        factor(category, levels = c("Female × Low", "Male × Medium", "Male × High")),
      variable == "Age × Education" ~
        factor(category, levels = c("Age × Low", "Age × Medium", "Age × High")),
      variable == "Gender × Age × Education" ~
        factor(category, levels = c("Female × Age × Low", "Male × Age × Medium",
                                    "Male × Age × High")),
      TRUE ~ factor(category)
    )
  )

interaction_summary <- coef_summary %>%
  filter(variable %in% c("Gender × Age", "Gender × Education",
                         "Age × Education", "Gender × Age × Education")) %>%
  mutate(
    sig_stars = case_when(
      is_reference ~ "",
      sig_pct >= 95 ~ "***",
      sig_pct >= 90 ~ "**",
      sig_pct >= 80 ~ "*",
      TRUE ~ ""
    ),
    variable = factor(variable,
                      levels = c("Gender × Age", "Gender × Education",
                                 "Age × Education", "Gender × Age × Education")),
    category = case_when(
      variable == "Gender × Age" ~
        factor(category, levels = c("Female × Age", "Male × Age")),
      variable == "Gender × Education" ~
        factor(category, levels = c("Female × Low", "Male × Medium", "Male × High")),
      variable == "Age × Education" ~
        factor(category, levels = c("Age × Low", "Age × Medium", "Age × High")),
      variable == "Gender × Age × Education" ~
        factor(category, levels = c("Female × Age × Low", "Male × Age × Medium",
                                    "Male × Age × High")),
      TRUE ~ factor(category)
    )
  )


################################################################################
## CALCULATE FACET-SPECIFIC Y-LIMITS
################################################################################

# Function to calculate reasonable limits per facet
calculate_facet_limits <- function(data, lower_q = 0.01, upper_q = 0.99) {
  data %>%
    filter(!is_reference) %>%
    group_by(variable) %>%
    summarise(
      y_min = quantile(estimate, lower_q, na.rm = TRUE),
      y_max = quantile(estimate, upper_q, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # Add some padding
      y_min = y_min - 0.1 * abs(y_min),
      y_max = y_max + 0.1 * abs(y_max)
    )
}

# Calculate limits for main effects
main_limits <- calculate_facet_limits(main_effects_data)
print("Main Effects Y-Limits:")
print(main_limits)

# Calculate limits for interactions
interaction_limits <- calculate_facet_limits(interaction_data)
print("\nInteraction Effects Y-Limits:")
print(interaction_limits)

################################################################################
## CALCULATE FACET-SPECIFIC Y-LIMITS
################################################################################

# Function to calculate reasonable limits per facet
calculate_facet_limits <- function(data, lower_q = 0.01, upper_q = 0.99) {
  data %>%
    filter(!is_reference) %>%
    group_by(variable) %>%
    summarise(
      y_min = quantile(estimate, lower_q, na.rm = TRUE),
      y_max = quantile(estimate, upper_q, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # Add some padding
      y_min = pmin(y_min, 0) - 0.05 * abs(pmax(y_max, abs(y_min))),
      y_max = pmax(y_max, 0) + 0.05 * abs(pmax(y_max, abs(y_min)))
    )
}

# Calculate limits for main effects
main_limits <- calculate_facet_limits(main_effects_data,  lower_q = 0.05, upper_q = 0.95)

# Calculate limits for interactions
interaction_limits <- calculate_facet_limits(interaction_data,  lower_q = 0.05, upper_q = 0.95)

################################################################################
## PLOT 1: MAIN EFFECTS - with facet-specific limits
################################################################################

# Add limits to data
main_effects_data_limited <- main_effects_data %>%
  left_join(main_limits, by = "variable") %>%
  mutate(
    estimate_display = case_when(
      estimate < y_min ~ y_min,
      estimate > y_max ~ y_max,
      TRUE ~ estimate
    ),
    is_clipped = estimate < y_min | estimate > y_max
  )

# Calculate median positions for stars
star_positions <- main_effects_data_limited %>%
  filter(!is_reference) %>%
  group_by(variable, category) %>%
  summarise(
    median_estimate = median(estimate_display, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    filter(main_effects_summary, !is_reference, sig_stars != ""),
    by = c("variable", "category")
  ) %>%
  filter(!is.na(sig_stars))

p_main <- ggplot(main_effects_data_limited, aes(x = category, y = estimate_display)) +
  # Reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.6) +

  # Non-reference categories: boxplot
  geom_boxplot(data = filter(main_effects_data_limited, !is_reference),
               width = 0.6,
               outlier.size = 1, outlier.alpha = 0.4, alpha = 0.7) +

  # Reference categories as points
  geom_point(data = filter(main_effects_data_limited, is_reference),
             color = "red3", size = 4, shape = 18) +

  # Reference label (only once)
  geom_text(data = filter(main_effects_data_limited, is_reference, variable == "Gender") %>%
              distinct(variable, category),
            aes(label = "Reference"),
            y = 0, color = "red3", size = 3, fontface = "italic",
            hjust = 0.5, vjust = -1.2) +

  # Mean points (white diamonds)
  stat_summary(data = filter(main_effects_data_limited, !is_reference),
               fun = mean, geom = "point", shape = 23,
               size = 3, fill = "white", color = "black", stroke = 0.8) +

  # Add significance stars
  geom_text(data = star_positions,
            aes(x = category, y = median_estimate, label = sig_stars),
            size = 5, fontface = "bold", vjust = -0.5, nudge_y = 0) +

  # Faceting with FREE scales
  facet_wrap(~variable, scales = "free", nrow = 1) +

  # Color scheme
  scale_fill_manual() +

  # Use blank to allow different limits per facet
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max)) +

  # Labels
  labs(
    x = NULL,
    y = NULL, #expression(paste("Regression Coefficient ", beta)),
    #title = "Main Effects: Which QIs Drive Inference Risk?",
    #caption = paste0("n = ", n_sim, " simulations; κ = ", kappa_fixed,
    #                 "; y-axes trimmed at 5th & 95th percentiles per facet\n",
    #                 "Significance: *** p<0.05 in ≥95% of sims, ** ≥90%, * ≥80%")
  )+

  # Theme
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(1.5, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(15, 15, 15, 15),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.caption = element_text(size = 9, hjust = 1, color = "gray40",
                                margin = margin(t = 10))
  )

#ggsave("sim3_main_effects.pdf", p_main, width = 11, height = 6, device = cairo_pdf)

################################################################################
## PLOT 2: INTERACTION EFFECTS - with facet-specific limits
################################################################################

# Add limits to data
interaction_data_limited <- interaction_data %>%
  left_join(interaction_limits, by = "variable") %>%
  mutate(
    estimate_display = case_when(
      estimate < y_min ~ y_min,
      estimate > y_max ~ y_max,
      TRUE ~ estimate
    ),
    is_clipped = estimate < y_min | estimate > y_max
  )

# Calculate median positions for stars
star_positions_interact <- interaction_data_limited %>%
  filter(!is_reference) %>%
  group_by(variable, category) %>%
  summarise(
    median_estimate = median(estimate_display, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    filter(interaction_summary, !is_reference, sig_stars != ""),
    by = c("variable", "category")
  ) %>%
  filter(!is.na(sig_stars))

p_interact <- ggplot(interaction_data_limited, aes(x = category, y = estimate_display)) +
  # Reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.6) +

  # Non-reference categories: boxplot
  geom_boxplot(data = filter(interaction_data_limited, !is_reference),
               width = 0.6,
               outlier.size = 1, outlier.alpha = 0.4, alpha = 0.7) +

  # Reference categories as points
  geom_point(data = filter(interaction_data_limited, is_reference),
             color = "red3", size = 4, shape = 18) +

  # Reference label (only once)
  geom_text(data = filter(interaction_data_limited, is_reference, variable == "Gender × Age") %>%
              distinct(variable, category),
            aes(label = ""),
            y = 0, color = "red3", size = 3, fontface = "italic",
            hjust = 0.5, vjust = -1.2) +

  # Mean points (white diamonds)
  stat_summary(data = filter(interaction_data_limited, !is_reference),
               fun = mean, geom = "point", shape = 23,
               size = 3, fill = "white", color = "black", stroke = 0.8) +

  # Add significance stars
  geom_text(data = star_positions_interact,
            aes(x = category, y = median_estimate, label = sig_stars),
            size = 5, fontface = "bold", vjust = -0.5, nudge_y = 0) +

  # Faceting with FREE scales
  facet_wrap(~variable, scales = "free", nrow = 1) +

  # Use blank to allow different limits per facet
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max)) +

  # Labels
  labs(
    x = NULL,
    y = NULL,#expression(paste("Regression Coefficient ", beta)),
    #title = "Interaction Effects: How QIs Work Together",
    #caption = paste0("n = ", n_sim, " simulations; κ = ", kappa_fixed,
    #                 "; y-axes trimmed at 5th & 95th percentiles per facet\n",
    #                 "Significance: *** p<0.05 in ≥95% of sims, ** ≥90%, * ≥80%")
  ) +

  # Theme
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    strip.text = element_text(face = "bold", size = 11),
    panel.spacing = unit(1.2, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(15, 15, 15, 15),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.caption = element_text(size = 9, hjust = 1, color = "gray40",
                                margin = margin(t = 10))
  )

#ggsave("sim3_interaction_effects.pdf", p_interact, width = 14, height = 6.5, device = cairo_pdf)


library(patchwork)
p_combined <- p_main / p_interact +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  ) &
  labs(y = expression(paste("Regression Coefficient ", beta)))
p_combined


library(cowplot)

# Beide Plots ohne Y-Titel (Caption bleibt im 2. Plot)
p_main_clean <- p_main + labs(y = NULL, title = NULL)
p_interact_clean <- p_interact + labs(y = NULL, title = NULL)

# Gemeinsame Y-Achsen-Beschriftung
y_label <- ggdraw() +
  draw_label(
    expression(paste("Regression Coefficient ", hat(beta))),
    angle = 90,
    size = 12,
    fontface = "bold"
  )

# Plots vertikal stacken
plots_stacked <- plot_grid(
  p_main_clean,
  p_interact_clean,
  ncol = 1,
  labels = c("Main", "Interaction"),
  label_size = 10
)

# Y-Achse + Plots kombinieren
p_combined <- plot_grid(
  y_label,
  plots_stacked,
  ncol = 2,
  rel_widths = c(0.05, 1)
)
p_combined
ggsave("sim3_qi_attribution_combined.pdf", p_combined, width = 14, height = 10, device = cairo_pdf)
