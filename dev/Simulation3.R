################################################################################
## SIMULATION 3: QI ATTRIBUTION ANALYSIS
## Method: Logistic regression coefficients with visible reference categories
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
## REGRESSION ANALYSIS FUNCTION
################################################################################

run_regression_analysis <- function(kappa = 10, n = 1000, tau = 0.3,
                                    synthesizer = "cart", attacker = "rf",
                                    sim = 1, seed_base = 5025) {

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

  # Create age categories
  records <- records %>%
    mutate(
      age_cat = cut(age, breaks = c(18, 35, 55, 90),
                    labels = c("Young", "Middle", "Old"))
    ) %>%
    mutate(
      # Set reference levels
      education = relevel(education, ref = "low"),
      gender = relevel(gender, ref = "female"),
      age_cat = relevel(age_cat, ref = "Young")
    )

  # FIT LOGISTIC REGRESSION
  logit_model <- glm(
    at_risk ~ gender + age_cat + education,
    data = records,
    family = binomial()
  )

  # Extract coefficients
  coef_results <- tidy(logit_model) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      variable = case_when(
        grepl("^gender", term) ~ "Gender",
        grepl("^age_cat", term) ~ "Age",
        grepl("^education", term) ~ "Education"
      ),
      category = case_when(
        term == "gendermale" ~ "male",
        term == "age_catMiddle" ~ "Middle",
        term == "age_catOld" ~ "Old",
        term == "educationmedium" ~ "medium",
        term == "educationhigh" ~ "high"
      ),
      is_reference = FALSE,
      sim = sim
    ) %>%
    select(sim, variable, category, estimate, std.error, p.value, is_reference)

  # ADD REFERENCE CATEGORIES (β = 0)
  reference_cats <- data.frame(
    sim = sim,
    variable = c("Gender", "Age", "Education"),
    category = c("female", "Young", "low"),
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

# Variable importance by coefficient range
var_importance <- coef_summary %>%
  group_by(variable) %>%
  summarise(
    beta_range = max(beta_mean) - min(beta_mean),
    max_category = category[which.max(beta_mean)],
    .groups = "drop"
  ) %>%
  arrange(desc(beta_range))

################################################################################
## VISUALIZE - IMPROVED VERSION
################################################################################

plot_data <- sim3_results %>%
  mutate(
    variable = factor(variable, levels = c("Education", "Age", "Gender")),
    category = case_when(
      variable == "Education" ~ factor(category, levels = c("low", "medium", "high")),
      variable == "Age" ~ factor(category, levels = c("Young", "Middle", "Old")),
      variable == "Gender" ~ factor(category, levels = c("female", "male")),
      TRUE ~ factor(category)
    )
  )

# Add significance stars to summary data
plot_summary <- coef_summary %>%
  mutate(
    sig_stars = case_when(
      is_reference ~ "",
      sig_pct >= 95 ~ "***",
      sig_pct >= 90 ~ "**",
      sig_pct >= 80 ~ "*",
      TRUE ~ ""
    ),
    variable = factor(variable, levels = c("Education", "Age", "Gender")),
    category = case_when(
      variable == "Education" ~ factor(category, levels = c("low", "medium", "high")),
      variable == "Age" ~ factor(category, levels = c("Young", "Middle", "Old")),
      variable == "Gender" ~ factor(category, levels = c("female", "male")),
      TRUE ~ factor(category)
    )
  )

# Create improved plot
p <- ggplot(plot_data, aes(x = category, y = estimate)) +
  # Reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.6) +

  # Non-reference categories: violin + boxplot
  geom_violin(data = filter(plot_data, !is_reference),
              aes(fill = variable), alpha = 0.4, trim = TRUE) +
  geom_boxplot(data = filter(plot_data, !is_reference),
               aes(fill = variable), width = 0.2, outlier.size = 0.8,
               outlier.alpha = 0.3, notch = TRUE, alpha = 0.8) +

  # Reference categories as points (all panels)
  geom_point(data = filter(plot_data, is_reference),
             color = "red3", size = 3.5, shape = 18) +

  # Reference label ONLY in Education panel (leftmost)
  geom_text(data = filter(plot_data, is_reference, variable == "Education") %>%
              distinct(variable, category),
            aes(label = "Reference"),
            y = 0, color = "red3", size = 3, fontface = "italic",
            hjust = 5.5, vjust = -.5) +

  # Add significance stars
  geom_text(data = filter(plot_summary, !is_reference, sig_stars != ""),
            aes(label = sig_stars, y = beta_mean + 0.05),
            size = 4, fontface = "bold", vjust = 0) +

  # Faceting
  facet_wrap(~variable, scales = "free_x", nrow = 1) +

  # Color scheme (colorblind-friendly)
  scale_fill_manual(
    values = c("Education" = "#0072B2",  # blue
               "Age" = "#E69F00",        # orange
               "Gender" = "#D55E00")     # red-orange
  ) +

  # Labels - minimal for paper figure
  labs(
    x = NULL,
    y = expression(paste("Regression Coefficient ", beta))
  ) +

  # Theme improvements
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
    strip.text = element_text(face = "bold", size = 11),
    panel.spacing = unit(1.5, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave("sim3_coefficients_with_references.pdf", p, width = 12, height = 5, device = cairo_pdf)
################################################################################
## SIMPLE BOXPLOTS
################################################################################

p_minimal <- ggplot(plot_data, aes(x = category, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.6) +

  geom_boxplot(data = filter(plot_data, !is_reference),
               fill = "gray70", width = 0.5, outlier.size = 1,
               outlier.alpha = 0.4, alpha = 0.5) +

  geom_point(data = filter(plot_data, is_reference),
             color = "red3", size = 4, shape = 18) +

  geom_text(data = filter(plot_data, is_reference, variable == "Education") %>%
              distinct(variable, category),
            aes(label = "Reference"),
            y = 0, color = "red3", size = 3, fontface = "italic",
            hjust = 5, vjust = -0.5) +

  # Mean points (white diamonds)
  stat_summary(data = filter(plot_data, !is_reference),
               fun = mean, geom = "point", shape = 23,
               size = 3, fill = "white", color = "black", stroke = 0.8) +

  facet_wrap(~variable, scales = "free_x", nrow = 1) +

  # Y-axis scale with sufficient range
  scale_y_continuous(
    breaks = seq(-20, 5, by = 2.5),
    limits = c(-20, 5)
  ) +

  labs(x = NULL, y = expression(paste("Regression Coefficient ", beta))) +

  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
    strip.text = element_text(face = "bold", size = 11),
    panel.spacing = unit(1.5, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave("plot_sim3_coefficients_minimal.pdf", p_minimal, width = 12, height = 5, device = cairo_pdf)



################################################################################
## TABLES WITH IMPROVEMENTS
################################################################################

table_sim3 <- coef_summary %>%
  mutate(
    `β (SD)` = ifelse(is_reference,
                      "0.000 (Ref)",
                      sprintf("%.3f (%.3f)", beta_mean, beta_sd)),
    `95% CI` = ifelse(is_reference,
                      "—",
                      sprintf("[%.3f, %.3f]", beta_lower, beta_upper)),
    `Sig %` = ifelse(is_reference, "—", sprintf("%.0f%%", sig_pct)),
    `N` = n_obs
  ) %>%
  select(Variable = variable, Category = category, `β (SD)`, `95% CI`, `Sig %`, N)

write.csv(table_sim3, "sim3_coefficients_table.csv", row.names = FALSE)

table_importance <- var_importance %>%
  mutate(
    `β Range` = sprintf("%.3f", beta_range),
    `Highest Risk` = max_category
  ) %>%
  select(Variable = variable, `β Range`, `Highest Risk`)

write.csv(table_importance, "sim3_variable_importance.csv", row.names = FALSE)
