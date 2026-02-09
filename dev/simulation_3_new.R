################################################################################
## SIMULATION 3: QI ATTRIBUTION ANALYSIS - FIGURE 5 (FINAL VERSION)
## Using predicted log-odds approach (like Matthias NEW style)
################################################################################

library(devtools)
load_all()
library(synthpop)
library(data.table)
library(ggplot2)
library(patchwork)
source("dev/simulate_microdata.R")

set.seed(2025)

################################################################################
## PARAMETERS
################################################################################
kappa_fixed <- 10
tau_fixed <- 0.3
n_sim <- 50
n_obs <- 1000

age_levels <- c(30, 50, 70)
age_center <- 50     # Fixed centering point

################################################################################
## RUN SIMULATIONS AND COLLECT PREDICTIONS
################################################################################
cat("Running", n_sim, "simulations...\n")
results_list <- vector("list", n_sim)

for (sim in 1:n_sim) {

  if (sim %% 10 == 0) cat("  Simulation", sim, "/", n_sim, "\n")

  # Generate data
  orig_data <- simulate_microdata(n = n_obs, dep = kappa_fixed, seed = 2025 + sim)

  syn_obj <- syn(orig_data, method = "cart", print.flag = FALSE,
                 seed = 2025 + sim + 1000)
  syn_data <- syn_obj$syn

  # Run RAPID
  rapid_result <- rapid(
    original_data = orig_data,
    synthetic_data = syn_data,
    quasi_identifiers = c("gender", "age", "education"),
    sensitive_attribute = "disease_status",
    model_type = "rf",
    cat_eval_method = "RCS_marginal",
    cat_tau = tau_fixed,
    return_all_records = TRUE,
    seed = 2025 + sim + 2000,
    trace = FALSE
  )

  # Extract records
  dat <- as.data.table(rapid_result$risk$rows_risk_df)

  # Set factors with nice labels
  dat[, gender := factor(gender, levels = c("female", "male"),
                         labels = c("Female", "Male"))]
  dat[, education := factor(education, levels = c("low", "medium", "high"),
                            labels = c("Low", "Medium", "High"))]

  # Standardize age: centered at 50, scaled per 10 years
  dat[, age_c := (age - age_center) / 10]

  # Fit model
  mod <- glm(at_risk ~ gender * education * age_c, data = dat, family = binomial)

  # ----------------------------------------------------------
  # Predicted grid for ALL combinations × age levels
  # ----------------------------------------------------------
  pred_grid <- as.data.table(expand.grid(
    gender = levels(dat$gender),
    education = levels(dat$education),
    age = age_levels
  ))
  pred_grid[, age_c := (age - age_center) / 10]
  pred_grid[, log_odds := predict(mod, newdata = pred_grid, type = "link")]
  pred_grid[, sim := sim]

  # ----------------------------------------------------------
  # Marginal effects (averaged over observed distribution)
  # ----------------------------------------------------------

  # Gender (marginalized over education and age)
  marg_gender <- dat[, .(
    Female = mean(predict(mod, newdata = copy(.SD)[, gender := "Female"], type = "link")),
    Male   = mean(predict(mod, newdata = copy(.SD)[, gender := "Male"],   type = "link"))
  )]
  marg_gender[, sim := sim]

  # Education (marginalized over gender and age)
  marg_education <- dat[, .(
    Low    = mean(predict(mod, newdata = copy(.SD)[, education := "Low"],    type = "link")),
    Medium = mean(predict(mod, newdata = copy(.SD)[, education := "Medium"], type = "link")),
    High   = mean(predict(mod, newdata = copy(.SD)[, education := "High"],   type = "link"))
  )]
  marg_education[, sim := sim]

  # Age effect: change in log-odds per 10 years
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
    predictions = pred_grid,
    marginal_gender = marg_gender,
    marginal_education = marg_education,
    marginal_age = marg_age
  )
}

cat("Done!\n\n")

################################################################################
## COMBINE RESULTS
################################################################################
all_preds <- rbindlist(lapply(results_list, `[[`, "predictions"))
all_marg_gender <- rbindlist(lapply(results_list, `[[`, "marginal_gender"))
all_marg_education <- rbindlist(lapply(results_list, `[[`, "marginal_education"))
all_marg_age <- rbindlist(lapply(results_list, `[[`, "marginal_age"))

################################################################################
## PLOTS - ROW 1: MAIN EFFECTS
################################################################################

# Gender
marg_gender_long <- melt(all_marg_gender, id.vars = "sim",
                         variable.name = "gender", value.name = "log_odds")
marg_gender_long[, gender := factor(gender, levels = c("Female", "Male"))]

p_gender <- ggplot(marg_gender_long, aes(x = gender, y = log_odds)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  labs(subtitle = "Gender", y = "Marginal Log-Odds", x = NULL) +
  theme_minimal(base_size = 11)

# Education
marg_edu_long <- melt(all_marg_education, id.vars = "sim",
                      variable.name = "education", value.name = "log_odds")
marg_edu_long[, education := factor(education, levels = c("Low", "Medium", "High"))]

p_edu <- ggplot(marg_edu_long, aes(x = education, y = log_odds)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  labs(subtitle = "Education", y = NULL, x = NULL) +
  theme_minimal(base_size = 11)

# Age effect
p_age <- ggplot(all_marg_age, aes(x = "Age\n(per 10y)", y = age_effect)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  labs(subtitle = "Age Effect", y = NULL, x = NULL) +
  theme_minimal(base_size = 11)

################################################################################
## PLOTS - ROW 2: 2-WAY INTERACTIONS
################################################################################

# A) Gender × Education (at Age = 50)
pred_age50 <- all_preds[age == 50]
pred_age50[, group := factor(
  paste(gender, education, sep = "\n"),
  levels = c("Female\nLow","Female\nMedium","Female\nHigh",
             "Male\nLow","Male\nMedium","Male\nHigh")
)]

p_gender_edu <- ggplot(pred_age50, aes(x = group, y = log_odds)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  coord_cartesian(ylim = c(-8, 5))+
  labs(subtitle = "Gender × Education\n(Age = 50)",
       y = "Predicted Log-Odds", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(size = 8, hjust = 0.5))

# B) Gender × Age (fix Education = Medium)
pred_gender_age <- all_preds[education == "Medium"]
pred_gender_age[, age_f := factor(age, levels = age_levels)]

p_gender_age <- ggplot(pred_gender_age,
                       aes(x = age_f, y = log_odds, fill = gender
                           )) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("Female" = "#E69F00", "Male" = "#56B4E9")) +
  coord_cartesian(ylim = c(-3, 3))+
  labs(subtitle = "Gender × Age\n(Education = Medium)",
       y = NULL, x = "Age") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank())

# C) Education × Age (fix Gender = Female)
pred_edu_age <- all_preds[gender == "Female"]
pred_edu_age[, age_f := factor(age, levels = age_levels)]

p_edu_age <- ggplot(pred_edu_age,
                    aes(x = age_f, y = log_odds ,fill = education
                        )) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, position = position_dodge(0.8)) +
  coord_cartesian(ylim = c(-10, 10))+
  scale_fill_brewer(palette = "YlOrRd") +
  #scale_fill_manual(values = c("Female" = "#E69F00", "Male" = "#56B4E9")) +
  labs(subtitle = "Education × Age\n(Gender = Female)",
       y = NULL, x = "Age") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", legend.title = element_blank())

################################################################################
## PLOTS - ROW 3: 3-WAY INTERACTION
################################################################################

p_3way <- ggplot(all_preds, aes(x = education, y = log_odds, fill = gender
                                )) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.3, position = position_dodge(0.8)) +
  coord_cartesian(ylim = c(-10, 10))+
  scale_fill_manual(values = c("Female" = "#E69F00", "Male" = "#56B4E9")) +
  facet_wrap(~ paste("Age =", age), nrow = 1) +
  labs(subtitle = "Gender × Education × Age",
       y = "Predicted Log-Odds", x = "Education") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        strip.text = element_text(face = "bold"))

################################################################################
## COMBINE ALL PLOTS - FIGURE 5
################################################################################

# Row 1: Main effects
p_main <- (p_gender | p_edu | p_age) +
  plot_layout(widths = c(1, 1.5, 0.9))

# Row 2: 2-way interactions
p_2way <- (p_gender_edu | p_gender_age | p_edu_age) +
  plot_layout(widths = c(1.5, 1, 1))

# Combined figure
p_combined <- p_main / p_2way / p_3way +
  plot_annotation(
    #title = "Quasi-Identifier Attribution Analysis",
    #subtitle = paste0("Predicted log-odds from logistic regression (n=", n_sim,
      #                " simulations, κ=", kappa_fixed, ", τ=", tau_fixed, ")"),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  ) +
  plot_layout(heights = c(1, 1.5, 1.3))

# Save
ggsave("assets/plot_sim3_qi_attribution_combined_new.pdf",
       p_combined, width = 14, height = 12, device = cairo_pdf)




# ################################################################################
# ## SIMULATION 3: QI ATTRIBUTION ANALYSIS - FIGURE 5
# ################################################################################
# library(devtools)
# load_all()
# library(synthpop)
# library(dplyr)
# library(ggplot2)
# library(tidyr)
# library(broom)
# library(data.table)
# library(patchwork)
# source("dev/simulate_microdata.R")
# set.seed(2025)
#
# ################################################################################
# ## PARAMETERS
# ################################################################################
# kappa_fixed <- 10
# tau_fixed <- 0.3
# n_sim <- 50
# n_obs <- 1000
#
# # Age levels for predictions (for gender×edu interaction plot)
# age_levels <- c(30, 50, 70)
#
# ################################################################################
# ## RUN SIMULATIONS AND COLLECT PREDICTIONS
# ################################################################################
#
# results_list <- list()
#
# for (sim in 1:n_sim) {
#
#   if (sim %% 10 == 0) cat("  Simulation", sim, "/", n_sim, "\n")
#
#   # Generate data
#   orig_data <- simulate_microdata(n = n_obs, dep = kappa_fixed, seed = 2025 + sim)
#
#   syn_obj <- syn(orig_data, method = "cart", print.flag = FALSE,
#                  seed = 2025 + sim + 1000)
#   syn_data <- syn_obj$syn
#
#   # Run RAPID
#   rapid_result <- rapid(
#     original_data = orig_data,
#     synthetic_data = syn_data,
#     quasi_identifiers = c("gender", "age", "education"),
#     sensitive_attribute = "disease_status",
#     model_type = "rf",
#     cat_eval_method = "RCS_marginal",
#     cat_tau = tau_fixed,
#     return_all_records = TRUE,
#     seed = 2025 + sim + 2000,
#     trace = FALSE
#   )
#
#   # Extract records
#   records <- rapid_result$risk$rows_risk_df
#
#   # Convert to data.table
#   dat <- as.data.table(records)
#
#   # Set reference levels
#   dat[, gender := factor(gender, levels = c("female", "male"))]
#   dat[, education := factor(education, levels = c("low", "medium", "high"))]
#
#   # Standardize age (center at mean)
#   age_mean <- mean(dat$age)
#   #age_sd <- sd(dat$age)
#   dat[, age_c := (age - age_mean) / 10]
#
#   # FIT MODEL
#   mod <- glm(at_risk ~ gender * education * age_c,
#              data = dat, family = binomial)
#
#   # ----------------------------------------------------------
#   # A) Extract coefficients
#   # ----------------------------------------------------------
#   coefs <- coef(mod)
#   se <- sqrt(diag(vcov(mod)))
#
#   coef_dt <- data.table(
#     term = names(coefs),
#     estimate = coefs,
#     se = se,
#     sim = sim
#   )
#
#   # ----------------------------------------------------------
#   # B) Compute predicted log-odds for grid
#   # ----------------------------------------------------------
#
#   pred_grid <- as.data.table(expand.grid(
#     gender = factor(c("female", "male"), levels = c("female", "male")),
#     education = factor(c("low", "medium", "high"), levels = c("low", "medium", "high")),
#     age = age_levels
#   ))
#
#   pred_grid[, age_c := (age - age_mean) / 10]
#   pred_grid[, log_odds := predict(mod, newdata = pred_grid, type = "link")]
#   pred_grid[, sim := sim]
#
#   # ----------------------------------------------------------
#   # C) Compute marginal effects
#   # ----------------------------------------------------------
#
#   # Gender (averaged over education and age)
#   marg_gender <- dat[, .(
#     female = mean(predict(mod, newdata = copy(.SD)[, gender := "female"], type = "link")),
#     male = mean(predict(mod, newdata = copy(.SD)[, gender := "male"], type = "link"))
#   )]
#   marg_gender[, sim := sim]
#
#   # Education (averaged over gender and age)
#   marg_education <- dat[, .(
#     low = mean(predict(mod, newdata = copy(.SD)[, education := "low"], type = "link")),
#     medium = mean(predict(mod, newdata = copy(.SD)[, education := "medium"], type = "link")),
#     high = mean(predict(mod, newdata = copy(.SD)[, education := "high"], type = "link"))
#   )]
#   marg_education[, sim := sim]
#
#   # Age effect (slope per SD)
#   marg_age <- dat[, {
#     pred_base <- predict(mod, newdata = .SD, type = "link")
#     newdata_plus <- copy(.SD)
#     newdata_plus[, age_c := age_c + 1]  # +1 SD
#     pred_plus <- predict(mod, newdata = newdata_plus, type = "link")
#     .(age_effect = mean(pred_plus - pred_base))
#   }]
#   marg_age[, sim := sim]
#
#   # Store results
#   results_list[[sim]] <- list(
#     coefficients = coef_dt,
#     predictions = pred_grid,
#     marginal_gender = marg_gender,
#     marginal_education = marg_education,
#     marginal_age = marg_age
#   )
# }
#
#
# # ============================================================
# # COMBINE RESULTS
# # ============================================================
#
# all_coefs <- rbindlist(lapply(results_list, `[[`, "coefficients"))
# all_preds <- rbindlist(lapply(results_list, `[[`, "predictions"))
# all_marg_gender <- rbindlist(lapply(results_list, `[[`, "marginal_gender"))
# all_marg_education <- rbindlist(lapply(results_list, `[[`, "marginal_education"))
# all_marg_age <- rbindlist(lapply(results_list, `[[`, "marginal_age"))
#
# # ============================================================
# # PLOT: MAIN EFFECTS
# # ============================================================
#
# # Gender
# marg_gender_long <- melt(all_marg_gender, id.vars = "sim",
#                          variable.name = "gender", value.name = "log_odds")
# marg_gender_long[, gender := factor(gender, levels = c("female", "male"),
#                                     labels = c("Female", "Male"))]
#
# p_gender <- ggplot(marg_gender_long, aes(x = gender, y = log_odds)) +
#   geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
#   labs(subtitle = "Gender", y = "Predicted Log-Odds", x = NULL) +
#   theme_minimal(base_size = 11)
#
# # Education
# marg_edu_long <- melt(all_marg_education, id.vars = "sim",
#                       variable.name = "education", value.name = "log_odds")
# marg_edu_long[, education := factor(education, levels = c("low", "medium", "high"),
#                                     labels = c("Low", "Medium", "High"))]
#
# p_edu <- ggplot(marg_edu_long, aes(x = education, y = log_odds)) +
#   geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
#   labs(subtitle = "Education", y = NULL, x = NULL) +
#   theme_minimal(base_size = 11)
#
# # Age slope
# p_age <- ggplot(all_marg_age, aes(x = "Age", y = age_effect)) +
#   geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
#   labs(subtitle = "Age Effect",
#        y = "Log-Odds Change\n(per 10 years)",  # ← Eigenes Label!
#        x = NULL) +
#   theme_minimal(base_size = 11)
#
# # ============================================================
# # PLOT: 2-WAY INTERACTIONS
# # ============================================================
#
# # A) Gender × Education (at age = 50)
# age_mean_closest <- age_levels[which.min(abs(age_levels - 50))]
# pred_gender_edu <- all_preds[age == age_mean_closest]
#
# pred_gender_edu[, gender_label := factor(gender, levels = c("female", "male"),
#                                          labels = c("Female", "Male"))]
# pred_gender_edu[, education_label := factor(education, levels = c("low", "medium", "high"),
#                                             labels = c("Low", "Medium", "High"))]
# pred_gender_edu[, group := paste(gender_label, education_label, sep = "\n")]
#
# pred_gender_edu[, group := factor(group, levels = c(
#   "Female\nLow", "Female\nMedium", "Female\nHigh",
#   "Male\nLow",   "Male\nMedium",   "Male\nHigh"
# ))]
#
# p_gender_edu <- ggplot(pred_gender_edu, aes(x = group, y = log_odds,# fill = gender_label
# )) +
#   geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
#   #scale_fill_manual(values = c("Female" = "#E69F00", "Male" = "#56B4E9")) +
#   labs(subtitle = paste0("Gender × Education\n(Age = ", age_mean_closest, ")"),
#        y = "Predicted Log-Odds", x = NULL) +
#   theme_minimal(base_size = 11) +
#   theme(legend.position = "none",
#         axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5))
#
# # B) Gender × Age (slopes)
# gender_age_slopes <- rbindlist(lapply(1:n_sim, function(s) {
#   coefs <- results_list[[s]]$coefficients
#
#   slope_female <- coefs[term == "age_c", estimate]
#
#   # Find male:age interaction
#   male_age_term <- coefs[grepl("gendermale.*age_c|age_c.*gendermale", term) &
#                            !grepl("education", term), ]
#
#   slope_male <- if(nrow(male_age_term) > 0) {
#     slope_female + male_age_term$estimate
#   } else {
#     slope_female  # No interaction found
#   }
#
#   data.table(sim = s, Female = slope_female, Male = slope_male)
# }))
#
# gender_age_long <- melt(gender_age_slopes, id.vars = "sim",
#                         variable.name = "gender", value.name = "age_slope")
#
# p_gender_age <- ggplot(gender_age_long, aes(x = gender, y = age_slope#, fill = gender
# )) +
#   geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
#   coord_cartesian(ylim = c(-10, 10)) +
#   #scale_fill_manual(values = c("Female" = "#E69F00", "Male" = "#56B4E9")) +
#   labs(subtitle = "Gender × Age",
#        y = NULL, x = NULL) +
#   theme_minimal(base_size = 11) +
#   theme(legend.position = "none")
#
# # C) Education × Age (slopes)
# edu_age_slopes <- rbindlist(lapply(1:n_sim, function(s) {
#   coefs <- results_list[[s]]$coefficients
#
#   slope_low <- coefs[term == "age_c", estimate]
#
#   # Find education:age interactions
#   #medium_age_term <- coefs[grepl("educationmedium.*age_c|age_c.*educationmedium", term) &
#   #                           !grepl("gender", term), ]
#
#   slope_medium_adj <- coefs[term == "educationmedium:age_c", estimate]
#   slope_medium <- slope_low + slope_medium_adj
#
#   #high_age_term <- coefs[grepl("educationhigh.*age_c|age_c.*educationhigh", term) &
#   #                         !grepl("gender", term), ]
#
#   slope_high_adj <- coefs[term == "educationhigh:age_c", estimate]
#   slope_high <- slope_low + slope_high_adj
#
#
#   data.table(sim = s, Low = slope_low, Medium = slope_medium, High = slope_high)
# }))
#
# edu_age_long <- melt(edu_age_slopes, id.vars = "sim",
#                      variable.name = "education", value.name = "age_slope")
# edu_age_long[, education := factor(education, levels = c("Low", "Medium", "High"))]
#
# p_edu_age <- ggplot(edu_age_long, aes(x = education, y = age_slope)) +
#   geom_boxplot(outlier.size = 0.5) +
#   coord_cartesian(ylim = c(-10, 10)) +  # ← CLIP!
#   labs(subtitle = "Education × Age",
#        y = "Age Effect (per 10 years)", x = NULL) +
#   theme_minimal(base_size = 11)
#
# # ============================================================
# # PLOT: 3-WAY INTERACTION
# # ============================================================
#
# gender_edu_age_slopes <- rbindlist(lapply(1:n_sim, function(s) {
#   coefs <- results_list[[s]]$coefficients
#
#   # Base slope (Female, Low)
#   base <- coefs[term == "age_c", estimate]
#
#   # 2-way adjustments
#   gender_adj <- coefs[grepl("gendermale.*age_c|age_c.*gendermale", term) &
#                         !grepl("education", term), estimate]
#   gender_adj <- ifelse(length(gender_adj) == 0, 0, gender_adj)
#
#   edu_med_adj <- coefs[grepl("educationmedium.*age_c|age_c.*educationmedium", term) &
#                          !grepl("gender", term), estimate]
#   edu_med_adj <- ifelse(length(edu_med_adj) == 0, 0, edu_med_adj)
#
#   edu_high_adj <- coefs[grepl("educationhigh.*age_c|age_c.*educationhigh", term) &
#                           !grepl("gender", term), estimate]
#   edu_high_adj <- ifelse(length(edu_high_adj) == 0, 0, edu_high_adj)
#
#   # 3-way adjustments
#   male_med_adj <- coefs[grepl("gendermale.*educationmedium.*age_c", term) |
#                           grepl("educationmedium.*gendermale.*age_c", term), estimate]
#   male_med_adj <- ifelse(length(male_med_adj) == 0, 0, male_med_adj)
#
#   male_high_adj <- coefs[grepl("gendermale.*educationhigh.*age_c", term) |
#                            grepl("educationhigh.*gendermale.*age_c", term), estimate]
#   male_high_adj <- ifelse(length(male_high_adj) == 0, 0, male_high_adj)
#
#   data.table(
#     sim = s,
#     `Female × Low` = base,
#     `Female × Medium` = base + edu_med_adj,
#     `Female × High` = base + edu_high_adj,
#     `Male × Low` = base + gender_adj,
#     `Male × Medium` = base + gender_adj + edu_med_adj + male_med_adj,
#     `Male × High` = base + gender_adj + edu_high_adj + male_high_adj
#   )
# }))
#
# gender_edu_age_long <- melt(gender_edu_age_slopes, id.vars = "sim",
#                             variable.name = "group", value.name = "age_slope")
#
# gender_edu_age_long[, gender := ifelse(grepl("Female", group), "Female", "Male")]
# gender_edu_age_long[, group := factor(group,
#                                       levels = c("Female × Low", "Female × Medium", "Female × High",
#                                                  "Male × Low", "Male × Medium", "Male × High"))]
#
# p_3way <- ggplot(gender_edu_age_long, aes(x = group, y = age_slope,# fill = gender
# )) +
#   geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
#   coord_cartesian(ylim = c(-10, 10)) +
#   #scale_fill_manual(values = c("Female" = "#E69F00", "Male" = "#56B4E9")) +
#   labs(subtitle = "Gender × Education × Age",
#        y = "Age Effect (per 10 years)", x = NULL) +
#   theme_minimal(base_size = 11) +
#   theme(legend.position = "bottom",
#         axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9))
#
# # ============================================================
# # COMBINE ALL PLOTS - FIGURE 5
# # ============================================================
#
# # Main effects row
# p_main <- (p_gender | p_edu | p_age) +
#   plot_layout(widths = c(1, 1.5, 0.8))
#
# # 2-way interactions row
# p_2way <- (p_edu_age  | p_gender_age | p_gender_edu) +#/  (p_gender_edu | plot_spacer()) +
#   plot_layout(heights = c(1, 1))
#
# # Full combined plot
# p_combined <- p_main / p_2way / p_3way +
#   plot_annotation(
#     title = "Quasi-Identifier Attribution Analysis",
#     subtitle = paste0("Predicted log-odds from logistic regression (n=", n_sim,
#                       " simulations, κ=", kappa_fixed, ", τ=", tau_fixed, ")"),
#     theme = theme(plot.title = element_text(size = 14, face = "bold"))
#   ) +
#   plot_layout(heights = c(1, 1.9, 1.3))
# p_combined
# # Save
# ggsave("dev/sim3_qi_attribution_combined.pdf",
#        p_combined, width = 14, height = 12, device = cairo_pdf)
#
# cat("\nFigure 5 saved to: assets/plot_sim3_qi_attribution_combined.pdf\n")
# print(p_combined)
