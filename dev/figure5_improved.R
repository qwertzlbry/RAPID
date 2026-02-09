# Matthias Templ, 05.02.2026
# Rapid Paper, first submission to CSDA

# ============================================================
# Figure 5: Improved Visualization of QI Attribution
# Including Age and all interactions (2-way and 3-way)
# ============================================================

library(ggplot2)
library(data.table)
library(patchwork)

set.seed(42)

# ============================================================
# 1. SIMULATION PARAMETERS
# ============================================================

n_sims <- 50      # Number of simulations
n_obs <- 1000     # Observations per simulation
kappa <- 10       # Dependency strength

# Age levels for visualization
age_levels <- c(30, 50, 70)

# ============================================================
# 2. DATA GENERATION FUNCTION
# ============================================================

generate_data <- function(n, kappa) {

  dat <- data.table(
    gender = factor(sample(c("Female", "Male"), n, replace = TRUE),
                    levels = c("Female", "Male")),
    education = factor(sample(c("Low", "Medium", "High"), n, replace = TRUE),
                       levels = c("Low", "Medium", "High")),
    age = round(rnorm(n, mean = 50, sd = 15))
  )

  # Standardize age for model (centered at 50)
  dat[, age_c := (age - 50) / 10]

  # True log-odds with all effects (scaled by kappa)
  # Higher kappa = stronger dependencies
  scale_factor <- sqrt(kappa / (1 + kappa))

  dat[, log_odds := scale_factor * (
    # Intercept (Female, Low, Age=50)
    -0.5 +

    # Main effects
    0.4 * (gender == "Male") +
    0.3 * (education == "Medium") +
    0.7 * (education == "High") +
    0.15 * age_c +

    # 2-way interactions: Gender × Education
    -0.2 * (gender == "Male") * (education == "Medium") +
    -0.1 * (gender == "Male") * (education == "High") +

    # 2-way interactions: Gender × Age
    0.1 * (gender == "Male") * age_c +

    # 2-way interactions: Education × Age
    -0.15 * (education == "Medium") * age_c +
    -0.25 * (education == "High") * age_c +

    # 3-way interaction: Gender × Education × Age
    0.05 * (gender == "Male") * (education == "Medium") * age_c +
    0.08 * (gender == "Male") * (education == "High") * age_c
  )]

  # Add noise and generate outcome
  dat[, prob := plogis(log_odds + rnorm(n, 0, 0.05))]
  dat[, at_risk := rbinom(n, 1, prob)]

  return(dat)
}

# ============================================================
# 3. RUN SIMULATIONS
# ============================================================
results_list <- list()

for (sim in 1:n_sims) {

  if (sim %% 10 == 0) cat("  Simulation", sim, "/", n_sims, "\n")

  # Generate data
  dat <- generate_data(n_obs, kappa)

  # Fit full model with all interactions
  mod <- glm(at_risk ~ gender * education * age_c,
             data = dat, family = binomial)

  # ----------------------------------------------------------
  # A) Extract coefficients
  # ----------------------------------------------------------
  coefs <- coef(mod)
  se <- sqrt(diag(vcov(mod)))

  coef_dt <- data.table(
    term = names(coefs),
    estimate = coefs,
    se = se,
    sim = sim
  )

  # ----------------------------------------------------------
  # B) Compute predicted log-odds for all group combinations
  # ----------------------------------------------------------

  # Grid: all combinations of gender, education, and age levels
  pred_grid <- as.data.table(expand.grid(
    gender = c("Female", "Male"),
    education = c("Low", "Medium", "High"),
    age = age_levels
  ))
  pred_grid[, age_c := (age - 50) / 10]

  # Predictions
  pred_grid[, log_odds := predict(mod, newdata = pred_grid, type = "link")]
  pred_grid[, prob := predict(mod, newdata = pred_grid, type = "response")]
  pred_grid[, sim := sim]

  # ----------------------------------------------------------
  # C) Compute marginal effects (averaged over other variables)
  # ----------------------------------------------------------

  # For marginal effects, we need to average over the other variables
  # using the observed distribution in the data

  # Gender effect (averaged over education and age)
  marg_gender <- dat[, .(
    Female = mean(predict(mod, newdata = copy(.SD)[, gender := "Female"], type = "link")),
    Male = mean(predict(mod, newdata = copy(.SD)[, gender := "Male"], type = "link"))
  )]
  marg_gender[, sim := sim]

  # Education effect (averaged over gender and age)
  marg_education <- dat[, .(
    Low = mean(predict(mod, newdata = copy(.SD)[, education := "Low"], type = "link")),
    Medium = mean(predict(mod, newdata = copy(.SD)[, education := "Medium"], type = "link")),
    High = mean(predict(mod, newdata = copy(.SD)[, education := "High"], type = "link"))
  )]
  marg_education[, sim := sim]

  # Age effect (slope, averaged over gender and education)
  # Compute as difference in prediction for age+10 vs age
  marg_age <- dat[, {
    pred_base <- predict(mod, newdata = .SD, type = "link")
    newdata_plus <- copy(.SD)
    newdata_plus[, age_c := age_c + 1]  # +10 years
    pred_plus <- predict(mod, newdata = newdata_plus, type = "link")
    .(age_effect = mean(pred_plus - pred_base))
  }]
  marg_age[, sim := sim]

  # Store results
  results_list[[sim]] <- list(
    coefficients = coef_dt,
    predictions = pred_grid,
    marginal_gender = marg_gender,
    marginal_education = marg_education,
    marginal_age = marg_age
  )
}

cat("Done!\n\n")

# ============================================================
# 4. COMBINE RESULTS
# ============================================================

all_coefs <- rbindlist(lapply(results_list, `[[`, "coefficients"))
all_preds <- rbindlist(lapply(results_list, `[[`, "predictions"))
all_marg_gender <- rbindlist(lapply(results_list, `[[`, "marginal_gender"))
all_marg_education <- rbindlist(lapply(results_list, `[[`, "marginal_education"))
all_marg_age <- rbindlist(lapply(results_list, `[[`, "marginal_age"))

# ============================================================
# 5. PLOT 1: ORIGINAL STYLE (Coefficients, Reference = 0)
# ============================================================

# Prepare coefficient data
coef_plot <- copy(all_coefs)
coef_plot[, term_clean := gsub("gender|education|age_c|:", " × ", term)]
coef_plot[, term_clean := gsub("^ × | × $", "", term_clean)]
coef_plot[term == "(Intercept)", term_clean := "Intercept"]

# Classify terms
coef_plot[, term_type := fcase(
  term == "(Intercept)", "Intercept",
  !grepl(":", term) & term != "(Intercept)", "Main",
  grepl(":", term) & !grepl(".*:.*:", term), "2-way",
  grepl(".*:.*:", term), "3-way"
)]

# Add reference categories (as 0)
ref_terms <- c("genderFemale", "educationLow")
ref_data <- rbindlist(lapply(1:n_sims, function(s) {
  data.table(
    term = ref_terms,
    estimate = 0,
    se = 0,
    sim = s,
    term_clean = c("Female", "Low"),
    term_type = "Main"
  )
}))

coef_plot_full <- rbind(coef_plot, ref_data, fill = TRUE)
coef_plot_full[, is_reference := term %in% ref_terms]

# Separate by type for faceting
coef_main <- coef_plot_full[term_type == "Main"]
coef_main[, variable := fcase(
  grepl("Female|Male", term_clean), "Gender",
  grepl("Low|Medium|High", term_clean), "Education",
  grepl("age", term, ignore.case = TRUE), "Age"
)]

coef_2way <- coef_plot_full[term_type == "2-way"]
coef_2way[, variable := fcase(
  grepl("gender.*education|education.*gender", term, ignore.case = TRUE), "Gender × Education",
  grepl("gender.*age|age.*gender", term, ignore.case = TRUE), "Gender × Age",
  grepl("education.*age|age.*education", term, ignore.case = TRUE), "Age × Education"
)]

coef_3way <- coef_plot_full[term_type == "3-way"]

# Plot main effects
p_main_old <- ggplot(coef_main[!is.na(variable)],
                     aes(x = term_clean, y = estimate, fill = is_reference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "gray70"), guide = "none") +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +
  labs(subtitle = "Main Effects", y = expression(hat(beta)), x = NULL) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 2-way interactions
p_2way_old <- ggplot(coef_2way[!is.na(variable)],
                     aes(x = term_clean, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(fill = "steelblue", alpha = 0.8, outlier.size = 0.5) +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +
  labs(subtitle = "2-Way Interactions", y = expression(hat(beta)), x = NULL) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 3-way interactions
p_3way_old <- ggplot(coef_3way, aes(x = term_clean, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(fill = "steelblue", alpha = 0.8, outlier.size = 0.5) +
  labs(subtitle = "3-Way Interactions", y = expression(hat(beta)), x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine old-style plot
p_old_combined <- p_main_old / p_2way_old / p_3way_old +
  plot_annotation(
    title = "ORIGINAL: Regression Coefficients (Reference Categories = 0)",
    subtitle = paste0("n = ", n_obs, ", ", n_sims, " simulations, kappa = ", kappa)
  )

# ============================================================
# 6. PLOT 2: NEW STYLE (Predicted Log-Odds, all groups have variance)
# ============================================================

# A) Marginal effects for main effects
marg_gender_long <- melt(all_marg_gender, id.vars = "sim",
                         variable.name = "gender", value.name = "log_odds")
marg_edu_long <- melt(all_marg_education, id.vars = "sim",
                      variable.name = "education", value.name = "log_odds")

p_gender_new <- ggplot(marg_gender_long, aes(x = gender, y = log_odds)) +
  geom_boxplot(fill = "forestgreen", alpha = 0.7, outlier.size = 0.5) +
  labs(subtitle = "Gender", y = "Marginal Log-Odds", x = NULL) +
  theme_minimal()

p_edu_new <- ggplot(marg_edu_long, aes(x = education, y = log_odds)) +
  geom_boxplot(fill = "forestgreen", alpha = 0.7, outlier.size = 0.5) +
  labs(subtitle = "Education", y = "Marginal Log-Odds", x = NULL) +
  theme_minimal()

p_age_new <- ggplot(all_marg_age, aes(x = "Age (per 10y)", y = age_effect)) +
  geom_boxplot(fill = "forestgreen", alpha = 0.7, outlier.size = 0.5) +
  labs(subtitle = "Age Effect", y = "Change in Log-Odds per 10 years", x = NULL) +
  theme_minimal()

# B) 2-way interactions: Gender × Education (at mean age)
pred_age50 <- all_preds[age == 50]
pred_age50[, group := paste(gender, education, sep = "\n")]

p_genderedu_new <- ggplot(pred_age50, aes(x = group, y = log_odds, fill = gender)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_manual(values = c("Female" = "#E69F00", "Male" = "#56B4E9")) +
  labs(subtitle = "Gender × Education (Age = 50)",
       y = "Predicted Log-Odds", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8))

# C) 2-way interactions: Gender × Age
pred_gender_age <- all_preds[education == "Medium"]  # Fix education at Medium
pred_gender_age[, age_f := factor(age)]

p_genderage_new <- ggplot(pred_gender_age, aes(x = age_f, y = log_odds, fill = gender)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("Female" = "#E69F00", "Male" = "#56B4E9")) +
  labs(subtitle = "Gender × Age (Education = Medium)",
       y = "Predicted Log-Odds", x = "Age") +
  theme_minimal() +
  theme(legend.position = "bottom")

# D) 2-way interactions: Education × Age
pred_edu_age <- all_preds[gender == "Female"]  # Fix gender at Female

p_eduage_new <- ggplot(pred_edu_age, aes(x = factor(age), y = log_odds, fill = education)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, position = position_dodge(0.8)) +
  scale_fill_brewer(palette = "Set2") +
  labs(subtitle = "Education × Age (Gender = Female)",
       y = "Predicted Log-Odds", x = "Age") +
  theme_minimal() +
  theme(legend.position = "bottom")

# E) 3-way interaction: Full grid
all_preds[, group_full := paste(gender, education, paste0("Age ", age), sep = "\n")]
all_preds[, age_f := factor(age)]

p_3way_new <- ggplot(all_preds, aes(x = education, y = log_odds, fill = gender)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.3, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("Female" = "#E69F00", "Male" = "#56B4E9")) +
  facet_wrap(~ paste("Age =", age), nrow = 1) +
  labs(subtitle = "Gender × Education × Age (Full Interaction)",
       y = "Predicted Log-Odds", x = "Education") +
  theme_minimal() +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

# Combine new-style plots
p_new_main <- (p_gender_new | p_edu_new | p_age_new) +
  plot_layout(widths = c(1, 1.5, 0.8))

p_new_2way <- (p_genderedu_new | p_genderage_new) /
  (p_eduage_new | plot_spacer())

p_new_combined <- p_new_main / p_3way_new +
  plot_annotation(
    title = "NEW: Predicted Log-Odds (All Groups Have Variance)",
    subtitle = paste0("n = ", n_obs, ", ", n_sims, " simulations, kappa = ", kappa)
  ) +
  plot_layout(heights = c(1, 1.5))

# ============================================================
# 7. PLOT 3: COMPARISON - Reference category variance
# ============================================================

# Show that reference category (Female, Low, Age=50) has variance
# when showing predicted values, but not when showing coefficients

# Predicted log-odds for reference group
ref_group_pred <- all_preds[gender == "Female" & education == "Low" & age == 50]

# Intercept from coefficients (= predicted log-odds for reference)
intercept_coef <- all_coefs[term == "(Intercept)"]

comparison_data <- rbind(
  data.table(type = "Intercept\n(from coefficients)",
             value = intercept_coef$estimate),
  data.table(type = "Predicted Log-Odds\n(Female, Low, Age=50)",
             value = ref_group_pred$log_odds)
)

p_comparison <- ggplot(comparison_data, aes(x = type, y = value)) +
  geom_boxplot(fill = c("steelblue", "forestgreen"), alpha = 0.7) +
  labs(title = "Reference Category: Same Values, Different Perspectives",
       subtitle = "Both show the log-odds for (Female, Low, Age=50)",
       y = "Log-Odds", x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10))

# ============================================================
# 8. SAVE PLOTS
# ============================================================

# Save individual plots
ggsave("dev/figure5_original_coefficients.pdf",
       p_old_combined, width = 12, height = 10)

ggsave("dev/igure5_new_predicted.pdf",
       p_new_combined, width = 12, height = 8)

ggsave("dev/figure5_comparison_reference.pdf",
       p_comparison, width = 6, height = 5)

# Combined comparison
p_final <- (p_old_combined | p_new_combined) +
  plot_annotation(
    title = "Figure 5: Comparison of Visualization Approaches",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

ggsave("dev/figure5_full_comparison.pdf",
       p_final, width = 20, height = 10)

cat("\nPlots saved to ~/Downloads/\n")
cat("  - figure5_original_coefficients.pdf\n")
cat("  - figure5_new_predicted.pdf\n")
cat("  - figure5_comparison_reference.pdf\n")
cat("  - figure5_full_comparison.pdf\n")

# ============================================================
# 9. PRINT SUMMARY STATISTICS
# ============================================================

cat("\n============================================================\n")
cat("SUMMARY: Reference Category Variance\n")
cat("============================================================\n\n")

cat("Intercept (= Reference Group Log-Odds):\n")
cat("  Mean:", round(mean(intercept_coef$estimate), 3), "\n")
cat("  SD:  ", round(sd(intercept_coef$estimate), 3), "\n")
cat("  Range:", round(min(intercept_coef$estimate), 3), "to",
    round(max(intercept_coef$estimate), 3), "\n\n")

cat("In the ORIGINAL plot, 'Female' and 'Low' appear as points at 0.\n")
cat("In the NEW plot, they show the same variance as the intercept.\n")

# ============================================================
# 10. DISPLAY PLOTS
# ============================================================

print(p_old_combined)
print(p_new_combined)
print(p_comparison)

# ============================================================
# 11. EFFECT CODING (Sum-to-Zero Contrasts)
# ============================================================

cat("\n============================================================\n")
cat("Running simulations with EFFECT CODING (Sum-to-Zero)...\n")
cat("============================================================\n\n")

results_effect <- list()

for (sim in 1:n_sims) {

  if (sim %% 10 == 0) cat("  Simulation", sim, "/", n_sims, "\n")

  # Generate fresh data
  dat <- generate_data(n_obs, kappa)

  # Set SUM-TO-ZERO contrasts (effect coding)
  contrasts(dat$gender) <- contr.sum(2)
  contrasts(dat$education) <- contr.sum(3)

  # Fit model with effect coding
  mod_effect <- glm(at_risk ~ gender * education * age_c,
                    data = dat, family = binomial)

  coefs <- coef(mod_effect)

  # With effect coding:
  # - Intercept = grand mean (across all groups)
  # - gender1 = effect of Female (Male = -Female)
  # - education1 = effect of Low, education2 = effect of Medium
  #   (High = -Low - Medium)

  # Extract and transform coefficients
  # Main effects
  eff_female <- coefs["gender1"]
  eff_male <- -coefs["gender1"]

  eff_low <- -coefs["education1"] - coefs["education2"]  # Sum to zero

  eff_medium <- coefs["education1"]
  eff_high <- coefs["education2"]

  eff_age <- coefs["age_c"]

  # 2-way: Gender x Education
  # gender1:education1 = Female x Low effect deviation
  # gender1:education2 = Female x Medium effect deviation
  # Male x Low = -Female x Low, etc.
  int_female_low <- -coefs["gender1:education1"] - coefs["gender1:education2"]
  int_female_medium <- coefs["gender1:education1"]
  int_female_high <- coefs["gender1:education2"]
  int_male_low <- -int_female_low
  int_male_medium <- -int_female_medium
  int_male_high <- -int_female_high

  # 2-way: Gender x Age
  int_female_age <- coefs["gender1:age_c"]
  int_male_age <- -coefs["gender1:age_c"]

  # 2-way: Education x Age
  int_low_age <- -coefs["education1:age_c"] - coefs["education2:age_c"]
  int_medium_age <- coefs["education1:age_c"]
  int_high_age <- coefs["education2:age_c"]

  # 3-way: Gender x Education x Age
  int3_female_low_age <- -coefs["gender1:education1:age_c"] - coefs["gender1:education2:age_c"]
  int3_female_medium_age <- coefs["gender1:education1:age_c"]
  int3_female_high_age <- coefs["gender1:education2:age_c"]
  int3_male_low_age <- -int3_female_low_age
  int3_male_medium_age <- -int3_female_medium_age
  int3_male_high_age <- -int3_female_high_age

  results_effect[[sim]] <- data.table(
    sim = sim,
    intercept = coefs["(Intercept)"],
    # Main effects
    Female = eff_female, Male = eff_male,
    Low = eff_low, Medium = eff_medium, High = eff_high,
    Age = eff_age,
    # 2-way Gender x Education
    `Female x Low` = int_female_low, `Female x Medium` = int_female_medium,
    `Female x High` = int_female_high,
    `Male x Low` = int_male_low, `Male x Medium` = int_male_medium,
    `Male x High` = int_male_high,
    # 2-way Gender x Age
    `Female x Age` = int_female_age, `Male x Age` = int_male_age,
    # 2-way Education x Age
    `Low x Age` = int_low_age, `Medium x Age` = int_medium_age,
    `High x Age` = int_high_age,
    # 3-way
    `Female x Low x Age` = int3_female_low_age,
    `Female x Medium x Age` = int3_female_medium_age,
    `Female x High x Age` = int3_female_high_age,
    `Male x Low x Age` = int3_male_low_age,
    `Male x Medium x Age` = int3_male_medium_age,
    `Male x High x Age` = int3_male_high_age
  )
}

effect_results <- rbindlist(results_effect)

cat("Done!\n\n")

# Reshape for plotting
effect_long <- melt(effect_results, id.vars = "sim",
                    variable.name = "term", value.name = "estimate")

# Classify terms
effect_long[, term_type := fcase(
  term == "intercept", "Intercept",
  term %in% c("Female", "Male", "Low", "Medium", "High", "Age"), "Main",
  grepl("^(Female|Male) x (Low|Medium|High)$", term), "2-way (Gender x Edu)",
  grepl("^(Female|Male) x Age$", term), "2-way (with Age)",
  grepl("^(Low|Medium|High) x Age$", term), "2-way (with Age)",
  grepl(" x .+ x ", term), "3-way"
)]

# ============================================================
# 12. PLOT: Effect Coding Results
# ============================================================

# Main effects
effect_main <- effect_long[term_type == "Main"]
effect_main[, variable := fcase(
  term %in% c("Female", "Male"), "Gender",
  term %in% c("Low", "Medium", "High"), "Education",
  term == "Age", "Age"
)]

p_effect_main <- ggplot(effect_main, aes(x = term, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(fill = "darkorange", alpha = 0.7, outlier.size = 0.5) +
  facet_wrap(~ variable, scales = "free_x", nrow = 1) +
  labs(subtitle = "Main Effects (all estimated, sum to zero within factor)",
       y = "Effect (deviation from grand mean)", x = NULL) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

# 2-way: Gender x Education
effect_2way_ge <- effect_long[term_type == "2-way (Gender x Edu)"]

p_effect_2way_ge <- ggplot(effect_2way_ge, aes(x = term, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(fill = "darkorange", alpha = 0.7, outlier.size = 0.5) +
  labs(subtitle = "2-Way: Gender x Education",
       y = "Interaction Effect", x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2-way with Age
effect_2way_age <- effect_long[term_type == "2-way (with Age)"]

p_effect_2way_age <- ggplot(effect_2way_age, aes(x = term, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(fill = "darkorange", alpha = 0.7, outlier.size = 0.5) +
  labs(subtitle = "2-Way Interactions with Age",
       y = "Interaction Effect", x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3-way
effect_3way <- effect_long[term_type == "3-way"]

p_effect_3way <- ggplot(effect_3way, aes(x = term, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(fill = "darkorange", alpha = 0.7, outlier.size = 0.5) +
  labs(subtitle = "3-Way: Gender x Education x Age",
       y = "Interaction Effect", x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine effect coding plots
p_effect_combined <- (p_effect_main) /
  (p_effect_2way_ge | p_effect_2way_age) /
  p_effect_3way +
  plot_annotation(
    title = "EFFECT CODING (Sum-to-Zero): All Categories Estimated",
    subtitle = paste0("n = ", n_obs, ", ", n_sims,
                      " simulations, kappa = ", kappa,
                      "\nNo reference category - all effects relative to grand mean")
  ) +
  plot_layout(heights = c(1, 1.2, 1))

# Save
ggsave("dev/figure5_effect_coding.pdf",
       p_effect_combined, width = 14, height = 12)

cat("Effect coding plot saved to:\n")
cat("  - figure5_effect_coding.pdf\n")

# ============================================================
# 13. COMPARISON: Three Approaches Side by Side
# ============================================================

# Summary statistics for main effects across approaches
cat("\n============================================================\n")
cat("COMPARISON: Main Effect Estimates Across Coding Schemes\n")
cat("============================================================\n\n")

cat("GENDER EFFECT:\n")
cat("  Dummy coding (Male coefficient):\n")
male_coef <- all_coefs[term == "genderMale", estimate]
cat("    Mean:", round(mean(male_coef), 3),
    " SD:", round(sd(male_coef), 3), "\n")

cat("  Effect coding:\n")
cat("    Female - Mean:", round(mean(effect_results$Female), 3),
    " SD:", round(sd(effect_results$Female), 3), "\n")
cat("    Male   - Mean:", round(mean(effect_results$Male), 3),
    " SD:", round(sd(effect_results$Male), 3), "\n")
cat("    (Note: Male = -Female, both estimated with variance)\n\n")

cat("EDUCATION EFFECT:\n")
cat("  Dummy coding:\n")
med_coef <- all_coefs[term == "educationMedium", estimate]
high_coef <- all_coefs[term == "educationHigh", estimate]
cat("    Low (reference): fixed at 0\n")
cat("    Medium - Mean:", round(mean(med_coef), 3),
    " SD:", round(sd(med_coef), 3), "\n")
cat("    High   - Mean:", round(mean(high_coef), 3),
    " SD:", round(sd(high_coef), 3), "\n")

cat("  Effect coding:\n")
cat("    Low    - Mean:", round(mean(effect_results$Low), 3),
    " SD:", round(sd(effect_results$Low), 3), "\n")
cat("    Medium - Mean:", round(mean(effect_results$Medium), 3),
    " SD:", round(sd(effect_results$Medium), 3), "\n")
cat("    High   - Mean:", round(mean(effect_results$High), 3),
    " SD:", round(sd(effect_results$High), 3), "\n")
cat("    (Note: Low + Medium + High = 0, all estimated with variance)\n")

# ============================================================
# 14. FINAL COMBINED PLOT: All Three Approaches
# ============================================================

p_all_three <- p_old_combined + p_new_combined + p_effect_combined +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Figure 5: Three Visualization Approaches Compared",
    subtitle = "Left: Dummy Coding | Center: Predicted Values | Right: Effect Coding"
  )

ggsave("dev/figure5_all_three_approaches.pdf",
       p_all_three, width = 30, height = 12)

cat("\n\nFinal comparison saved to:\n")
cat("  - figure5_all_three_approaches.pdf\n")

print(p_effect_combined)
