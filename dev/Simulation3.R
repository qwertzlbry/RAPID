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

  # Set reference levels (keep age continuous)
  records <- records %>%
    mutate(
      education = relevel(education, ref = "low"),
      gender = relevel(gender, ref = "female")
    )

  # FIT LOGISTIC REGRESSION with age × education interaction
  logit_model <- glm(
    at_risk ~ gender + age + education + age:education,
    data = records,
    family = binomial()
  )

  # Extract coefficients
  coef_results <- tidy(logit_model) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      variable = case_when(
        grepl("^age:education", term) ~ "Age × Education",
        grepl("^gender", term) ~ "Gender",
        term == "age" ~ "Age",
        grepl("^education", term) ~ "Education",
        TRUE ~ "Other"  # catch-all to avoid NAs
      ),
      category = case_when(
        term == "gendermale" ~ "Male",
        term == "age" ~ "Per year",
        term == "educationmedium" ~ "Medium",
        term == "educationhigh" ~ "High",
        term == "age:educationmedium" ~ "Age:Medium",
        term == "age:educationhigh" ~ "Age:High",
        TRUE ~ term  # fallback
      ),
      is_reference = FALSE,
      sim = sim
    ) %>%
    select(sim, variable, category, estimate, std.error, p.value, is_reference)

  # ADD REFERENCE CATEGORIES (β = 0)
  reference_cats <- data.frame(
    sim = sim,
    variable = c("Gender", "Education", "Age × Education"),
    category = c("Female", "Low", "Age:Low"),
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

# Variable importance by coefficient magnitude
var_importance <- coef_summary %>%
  filter(!is_reference) %>%
  group_by(variable) %>%
  summarise(
    max_abs_beta = max(abs(beta_mean)),
    strongest_category = category[which.max(abs(beta_mean))],
    .groups = "drop"
  ) %>%
  arrange(desc(max_abs_beta))

################################################################################
## PREPARE PLOT DATA
################################################################################
# Separate main effects and interactions
main_effects_data <- sim3_results %>%
  filter(variable %in% c("Education", "Age", "Gender")) %>%
  mutate(
    variable = factor(variable, levels = c("Education", "Age", "Gender")),
    category = case_when(
      variable == "Education" ~ factor(category, levels = c("Low", "Medium", "High")),
      variable == "Age" ~ factor(category, levels = c("Per year")),
      variable == "Gender" ~ factor(category, levels = c("Female", "Male")),
      TRUE ~ factor(category)
    )
  )

interaction_data <- sim3_results %>%
  filter(variable == "Age × Education") %>%
  mutate(
    category = factor(category, levels = c("Age:Low", "Age:Medium", "Age:High"))
  )

# Summary data for main effects
main_summary <- coef_summary %>%
  filter(variable %in% c("Education", "Age", "Gender")) %>%
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
      variable == "Education" ~ factor(category, levels = c("Low", "Medium", "High")),
      variable == "Age" ~ factor(category, levels = c("Per year")),
      variable == "Gender" ~ factor(category, levels = c("Female", "Male")),
      TRUE ~ factor(category)
    )
  )

# Summary data for interactions
interaction_summary <- coef_summary %>%
  filter(variable == "Age × Education") %>%
  mutate(
    sig_stars = case_when(
      is_reference ~ "",
      sig_pct >= 95 ~ "***",
      sig_pct >= 90 ~ "**",
      sig_pct >= 80 ~ "*",
      TRUE ~ ""
    ),
    category = factor(category, levels = c("Age:Low", "Age:Medium", "Age:High"))
  )

################################################################################
## PLOT 1: MAIN EFFECTS ONLY
################################################################################
p_main <- ggplot(main_effects_data, aes(x = category, y = estimate)) +
  # Reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.6) +

  # Non-reference categories: boxplot
  geom_boxplot(data = filter(main_effects_data, !is_reference),
               aes(fill = variable), width = 0.5, outlier.size = 1,
               outlier.alpha = 0.4, alpha = 0.6) +

  # Reference categories as points
  geom_point(data = filter(main_effects_data, is_reference),
             color = "red3", size = 4, shape = 18) +

  # Reference label
  geom_text(data = filter(main_effects_data, is_reference, variable == "Education") %>%
              distinct(variable, category),
            aes(label = "Reference"),
            y = 0, color = "red3", size = 3, fontface = "italic",
            hjust = 5, vjust = -0.5) +

  # Mean points (white diamonds)
  stat_summary(data = filter(main_effects_data, !is_reference),
               fun = mean, geom = "point", shape = 23,
               size = 3, fill = "white", color = "black", stroke = 0.8) +

  # Add significance stars
  geom_text(data = filter(main_summary, !is_reference, sig_stars != ""),
            aes(label = sig_stars, y = beta_mean),
            size = 4, fontface = "bold", vjust = -0.5) +

  # Faceting
  facet_wrap(~variable, scales = "free", nrow = 1) +

  # Color scheme (colorblind-friendly)
  scale_fill_manual(
    values = c("Education" = "#0072B2",  # blue
               "Age" = "#E69F00",        # orange
               "Gender" = "#D55E00")     # red-orange
  ) +

  # Labels
  labs(
    x = NULL,
    y = expression(paste("Regression Coefficient ", beta)),
    title = "Main Effects of QIs on Inference Risk"
  ) +

  # Theme
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
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
  )

ggsave("sim3_main_effects.pdf", p_main, width = 10, height = 5, device = cairo_pdf)

################################################################################
## PLOT 2: INTERACTION EFFECTS
################################################################################
p_interaction <- ggplot(interaction_data, aes(x = category, y = estimate)) +
  # Reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.6) +

  # Non-reference categories: boxplot
  geom_boxplot(data = filter(interaction_data, !is_reference),
               fill = "#009E73", width = 0.5, outlier.size = 1,
               outlier.alpha = 0.4, alpha = 0.6) +

  # Reference category as point
  geom_point(data = filter(interaction_data, is_reference),
             color = "red3", size = 4, shape = 18) +

  # Reference label
  geom_text(data = filter(interaction_data, is_reference) %>% distinct(category),
            aes(label = "Reference"),
            y = 0, color = "red3", size = 3, fontface = "italic",
            hjust = 0.5, vjust = -0.5) +

  # Mean points
  stat_summary(data = filter(interaction_data, !is_reference),
               fun = mean, geom = "point", shape = 23,
               size = 3, fill = "white", color = "black", stroke = 0.8) +

  # Add significance stars
  geom_text(data = filter(interaction_summary, !is_reference, sig_stars != ""),
            aes(label = sig_stars, y = beta_mean),
            size = 4, fontface = "bold", vjust = -0.5) +

  # Labels
  labs(
    x = "Age × Education Interaction",
    y = expression(paste("Interaction Coefficient ", beta)),
    title = "Age × Education Interaction Effects"
  ) +

  # Theme
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
  )

ggsave("sim3_interaction_effects.pdf", p_interaction, width = 6, height = 5, device = cairo_pdf)

################################################################################
## TABLES WITH IMPROVEMENTS
################################################################################
# Main effects table
table_main <- coef_summary %>%
  filter(variable %in% c("Education", "Age", "Gender")) %>%
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

write.csv(table_main, "sim3_main_effects_table.csv", row.names = FALSE)

# Interaction table
table_interaction <- coef_summary %>%
  filter(variable == "Age × Education") %>%
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
  select(Category = category, `β (SD)`, `95% CI`, `Sig %`, N)

write.csv(table_interaction, "sim3_interaction_table.csv", row.names = FALSE)

# Variable importance table
table_importance <- var_importance %>%
  mutate(
    `|β| Max` = sprintf("%.3f", max_abs_beta),
    `Strongest Category` = strongest_category
  ) %>%
  select(Variable = variable, `|β| Max`, `Strongest Category`)

write.csv(table_importance, "sim3_variable_importance.csv", row.names = FALSE)


