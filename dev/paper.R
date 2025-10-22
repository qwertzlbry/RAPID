#  source("~/workspace25/SDG_Comparison2/code/auxiliary_functions/risk/_risk_inferential.R", echo = TRUE)
# source("~/workspace25/SDG_Comparison2/code/auxiliary_functions/risk/handle_known_na.R", echo = TRUE)
# load function 
source("code/auxiliary_functions/risk/_risk_inferential.R", echo = TRUE)
# source("code/auxiliary_functions/risk/handle_known_na.R", echo = TRUE)


library(ranger)

simulate_microdata <- function(
    n = 1000,
    seed = NULL,
    # ----- dependency scaling (0 = independent; 1 = default; >1 = stronger) -----
    dep = 1,
    # ----- AGE -----
    age_mean = 45, age_sd = 15, age_min = 18, age_max = 90,
    # ----- EDUCATION (ordered: low < medium < high); latent ~ N(mu, sd^2) -----
    edu_age_beta = -0.15,      # younger cohorts more likely high ed (negative wrt age_z)
    edu_latent_sd = 1.0,
    edu_cutpoints = c(-0.3, 0.7),  # cutpoints for low|med|high on the latent scale
    # ----- INCOME (log-income model) -----
    log_income_mean = log(40000),  # baseline mean income around 40k
    income_beta_age = 0.30,        # effect of age_z
    income_beta_age2 = -0.15,      # life-cycle concavity
    income_beta_edu = 0.35,        # per education step (low=0, med=1, high=2)
    income_noise_sd = 0.35,        # residual sd on log-income
    # ----- HEALTH SCORE (0..100; higher = better) -----
    health_mean = 70,
    health_beta_age = -8.0,        # per SD of age (z-scale)
    health_beta_log_income = 6.0,  # per SD of log-income (z-scale)
    health_beta_edu = 4.0,         # per education step
    health_noise_sd = 10.0,
    # ----- DISEASE STATUS (multinomial: healthy, diabetic, hypertensive) -----
    # Baseline class = "healthy"; intercepts target ~70/15/15 when predictors at mean
    disease_intercepts = c(diabetic = -1.54, hypertensive = -1.54),
    disease_beta_age = c(diabetic = 0.80, hypertensive = 1.00),       # age_z effects
    disease_beta_log_income = c(diabetic = -0.30, hypertensive = -0.20),# log_inc_z
    disease_beta_health = c(diabetic = -1.20, hypertensive = -1.00),    # (health - mean)/15
    disease_beta_edu = c(diabetic = -0.20, hypertensive = -0.10),       # per edu step
    # ----- miscellaneous -----
    return_numeric_edu = FALSE
) {
  if (!is.null(seed)) set.seed(seed)
  
  # --- helpers ---
  clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))
  rtruncnorm_base <- function(n, mean, sd, lo, hi) {
    x <- rnorm(n, mean, sd)
    oob <- which(x < lo | x > hi)
    while (length(oob) > 0) {
      x[oob] <- rnorm(length(oob), mean, sd)
      oob <- which(x < lo | x > hi)
    }
    x
  }
  draw_multinomial <- function(p_mat, labels) {
    # p_mat: n x K matrix of probabilities
    u <- runif(nrow(p_mat))
    cs <- t(apply(p_mat, 1, cumsum))
    idx <- rowSums(u > cs) + 1L
    factor(labels[idx], levels = labels)
  }
  
  # 1) AGE ----------------------------------------------------------------------
  age <- rtruncnorm_base(n, age_mean, age_sd, age_min, age_max)
  age_z <- (age - age_mean) / age_sd
  
  # 2) EDUCATION (ordinal via latent variable + cutpoints) ----------------------
  # latent mean depends on age_z; sd = edu_latent_sd
  edu_latent_mu <- dep * (edu_age_beta * age_z)
  edu_latent <- rnorm(n, mean = edu_latent_mu, sd = edu_latent_sd)
  
  # cut into low/medium/high
  cut1 <- edu_cutpoints[1]; cut2 <- edu_cutpoints[2]
  edu_num <- ifelse(edu_latent < cut1, 0L,
                    ifelse(edu_latent < cut2, 1L, 2L))
  education <- factor(c("low","medium","high")[edu_num + 1L],
                      levels = c("low","medium","high"))
  
  # 3) INCOME (log-linear; depends on age_z, age_z^2, education) ---------------
  log_inc <- log_income_mean +
    dep * (income_beta_age  * age_z +
             income_beta_age2 * (age_z^2) +
             income_beta_edu  * edu_num) +
    rnorm(n, sd = income_noise_sd)
  
  income <- as.numeric(exp(log_inc))  # positively skewed, realistic
  # standardize log-income for downstream models
  log_inc_z <- (log_inc - log_income_mean) / income_noise_sd
  
  # 4) HEALTH SCORE (0..100; higher = better) ----------------------------------
  # Use z-scales for age and log-income; education as steps (0,1,2)
  health_raw <- health_mean +
    dep * (health_beta_age        * age_z +
             health_beta_log_income * log_inc_z +
             health_beta_edu        * edu_num) +
    rnorm(n, sd = health_noise_sd)
  
  health_score <- clamp(health_raw, 0, 100)
  # modest scale for health in disease logits
  health_z <- (health_score - health_mean) / 15
  
  # 5) DISEASE STATUS (multinomial logistic vs "healthy") ----------------------
  # linear predictors for diabetic & hypertensive relative to healthy
  eta_diab <- disease_intercepts["diabetic"] +
    dep * (disease_beta_age["diabetic"]         * age_z +
             disease_beta_log_income["diabetic"]  * log_inc_z +
             disease_beta_health["diabetic"]      * health_z +
             disease_beta_edu["diabetic"]         * edu_num)
  
  eta_hyp  <- disease_intercepts["hypertensive"] +
    dep * (disease_beta_age["hypertensive"]        * age_z +
             disease_beta_log_income["hypertensive"] * log_inc_z +
             disease_beta_health["hypertensive"]     * health_z +
             disease_beta_edu["hypertensive"]        * edu_num)
  
  # softmax with healthy as baseline (eta_healthy = 0)
  num_diab <- exp(eta_diab)
  num_hyp  <- exp(eta_hyp)
  denom    <- 1 + num_diab + num_hyp
  p_healthy <- 1 / denom
  p_diab    <- num_diab / denom
  p_hyp     <- num_hyp  / denom
  
  probs <- cbind(healthy = p_healthy,
                 diabetic = p_diab,
                 hypertensive = p_hyp)
  
  disease_status <- draw_multinomial(probs, labels = c("healthy","diabetic","hypertensive"))
  
  # gender
  
  # Gender: conditional Bernoulli with logistic model
  # Adjusted coefficients based on variable scales and desired balance
  eta_gender <- -3 + 
    0.03 * (age - 45) + 
    0.6  * edu_num + 
    0.00004 * income - 
    0.02 * (health_score - 70)
  
  prob_male <- 1 / (1 + exp(-eta_gender))
  gender <- factor(rbinom(n, 1, prob = prob_male), labels = c("female", "male"))
  
  # OUTPUT ----------------------------------------------------------------------
  out <- data.frame(
    gender = as.factor(gender),
    age = age,
    education = as.factor(education),
    income = income,
    health_score = health_score,
    disease_status = as.factor(disease_status)
  )
  if (return_numeric_edu) out$education_num <- edu_num
  out
}

# new version (test)
simulate_microdata <- function(
    n = 1000,
    seed = NULL,
    dep = 1,  # signal-to-noise scaling; dep = 1 gives ~50% explained variance
    age_mean = 45,
    age_sd = 15,
    age_min = 18,
    age_max = 90,
    return_latent = FALSE
) {
  if (!is.null(seed)) set.seed(seed)
  
  # ----- Dependency Scaling -----
  # dep = 0 → no signal (noise only); dep = Inf → fully deterministic
  signal_weight <- sqrt(dep / (1 + dep))
  noise_weight  <- sqrt(1 / (1 + dep))
  
  # ----- AGE (independent) -----
  rtruncnorm <- function(n, mean, sd, lo, hi) {
    x <- rnorm(n, mean, sd)
    while (any(out <- x < lo | x > hi)) {
      x[out] <- rnorm(sum(out), mean, sd)
    }
    x
  }
  age <- rtruncnorm(n, mean = age_mean, sd = age_sd, lo = age_min, hi = age_max)
  age_z <- scale(age)[, 1]
  
  # ----- LATENT SES (shared dependency source) -----
  latent_SES <- rnorm(n)
  
  # ----- EDUCATION (ordinal: low/med/high via latent + thresholds) -----
  edu_latent <- signal_weight * (0.8 * latent_SES - 0.4 * age_z) + noise_weight * rnorm(n)
  education_num <- cut(edu_latent, breaks = c(-Inf, -0.3, 0.7, Inf), labels = c(0,1,2)) |> as.integer() - 1
  education <- factor(c("low", "medium", "high")[education_num + 1], levels = c("low", "medium", "high"))
  
  # ----- INCOME (log-linear) -----
  lin_income <- signal_weight * (0.5 * latent_SES + 0.3 * age_z + 0.25 * education_num) +
    noise_weight * rnorm(n)
  income <- exp(10 + lin_income)  # ~ realistic log-income (centered around ~20-40k)
  
  # ----- HEALTH SCORE (bounded: 0–100, sigmoid link) -----
  lin_health <- signal_weight * (0.6 * latent_SES - 0.5 * age_z + 0.2 * education_num + 0.2 * scale(log(income))[,1]) +
    noise_weight * rnorm(n)
  health_score <- 100 / (1 + exp(-lin_health))  # sigmoid scale
  
  # ----- DISEASE STATUS (softmax over logits) -----
  health_z <- scale(health_score)[, 1]
  log_income_z <- scale(log(income))[, 1]
  
  # NEW disease logits: scaled linearly with dep (not sqrt)
  eta_diab <- -1.5 + dep * (
    0.8 * age_z - 0.3 * log_income_z - 0.2 * education_num
  )
  eta_hyp <- -1.3 + dep * (
    1.0 * age_z - 0.2 * log_income_z - 0.1 * education_num
  )
  eta_healthy <- 0  # baseline
  
  # softmax
  exp_eta <- cbind(
    healthy = exp(eta_healthy),
    diabetic = exp(eta_diab),
    hypertensive = exp(eta_hyp)
  )
  prob_mat <- exp_eta / rowSums(exp_eta)
  draw_multinom <- function(p, levels) {
    idx <- apply(p, 1, function(prob) sample(seq_along(prob), size = 1, prob = prob))
    factor(levels[idx], levels = levels)
  }
  disease_status <- draw_multinom(prob_mat, c("healthy", "diabetic", "hypertensive"))
  
  # ----- GENDER (optional, mild dependency) -----
  eta_gender <- signal_weight * (0.3 * latent_SES - 0.2 * age_z + 0.2 * education_num)
  prob_male <- 1 / (1 + exp(-eta_gender))
  gender <- factor(rbinom(n, 1, prob_male), labels = c("female", "male"))
  
  # ----- OUTPUT -----
  out <- data.frame(
    gender = gender,
    age = age,
    education = education,
    income = income,
    health_score = health_score,
    disease_status = disease_status
  )
  if (return_latent) {
    out$latent_SES <- latent_SES
    out$edu_latent <- edu_latent
  }
  return(out)
}


truth <- simulate_microdata(dep = 1, n = 1000)
synth <- synthpop::syn(truth, m = 1, seed = 123)$syn

# check visually
x <- simulate_microdata(dep = 100); pairs(x, col = as.integer(x$gender))


sapply(list.files("code/auxiliary_functions/risk_paper", pattern = "\\.R$", full.names = TRUE), source)

rapid_result <- rapid(
  original_data = truth,
  synthetic_data = synth,
  quasi_identifiers = c("age", "income", "education"),
  sensitive_attribute = "disease_status",
  model_type = "rf",
  tau = 1.25,
  trace = FALSE
)

rapid_result$risk$confidence_rate

# vary dep from 0 (independent) to 2 (stronger dep)
library(dplyr)
library(progressr)
library(purrr)
library(tibble)

dep_values <- seq(0, 100, by = 0.5)
n_rep <- 3

handlers(global = TRUE)  # enable progress bars (or choose a specific handler below)
# handlers("cli")       # for modern RStudio progress bar
# handlers("txtprogressbar") # fallback terminal progress bar
# handlers("beepr")     # tik-tok sounds, if you want audio

# Total number of iterations
total_steps <- length(dep_values) * n_rep

set.seed(42)
sim_seeds <- sample.int(1e6, size = total_steps)
# Use progressor
results <- with_progress({
  p <- progressor(steps = total_steps)
  
  map_dfr(dep_values, function(dep) {
    map_dfr(1:n_rep, function(rep_id) {
      p(message = sprintf("dep=%.1f, rep=%d", dep, rep_id))  # update progress bar
      
      rep_index <- (which(dep_values == dep) - 1) * n_rep + rep_id
      myseed <- sim_seeds[rep_index]
      truth <- simulate_microdata(dep = dep, n = 1000, seed = myseed)
      synth <- synthpop::syn(truth, m = 1, print.flag = FALSE, seed = myseed + 1)$syn
      
      rapid_result <- rapid(
        original_data = truth,
        synthetic_data = synth,
        quasi_identifiers = c("age", "income", "education"),
        sensitive_attribute = "disease_status",
        model_type = "rf",
        tau = 1.1,
        trace = FALSE
      )
      
        # # health_score R2
        # lm_health <- lm(health_score ~ age + income + education, data = x)
        # R2_health <- summary(lm_health)$r.squared
        
        # disease_status pseudo R2
        mod_disease <- nnet::multinom(disease_status ~ age + income + education + health_score, data = synth, trace = FALSE)
        ll_null <- logLik(update(mod_disease, . ~ 1))
        ll_full <- logLik(mod_disease)
        R2_disease <- 1 - as.numeric(ll_full) / as.numeric(ll_null)
      
      
      tibble(
        rep = rep_id,
        dep = dep,
        confidence_rate = rapid_result$risk$confidence_rate,
        accuracy = rapid_result$metrics$accuracy,
        R2 = R2_disease
      )
    })
  })
})
results

summary_results <- results %>%
  group_by(dep) %>%
  summarise(
    mean_conf = mean(confidence_rate),
    sd_conf = sd(confidence_rate),
    mean_acc = mean(accuracy),
    sd_acc = sd(accuracy),
    .groups = "drop"
  )

library(ggplot2)
ggplot(summary_results, aes(x = dep)) +
  # Confidence Rate
  geom_line(aes(y = mean_conf), color = "steelblue", size = 1.2) +
  geom_ribbon(aes(ymin = mean_conf - sd_conf,
                  ymax = mean_conf + sd_conf),
              fill = "steelblue", alpha = 0.2) +
  
  # Accuracy
  geom_line(aes(y = mean_acc), color = "firebrick", size = 1.2) +
  geom_ribbon(aes(ymin = mean_acc - sd_acc,
                  ymax = mean_acc + sd_acc),
              fill = "firebrick", alpha = 0.2) +
  
  labs(
    title = "Attribute Disclosure Risk and Accuracy vs Dependency Strength",
    subtitle = "Mean ± 1 SD over repeated simulations (n = 100 per dep)",
    x = expression(kappa~"(dependency strength)"),
    y = "Rate",
    caption = "Metrics: Confidence Rate (steelblue), Accuracy (firebrick)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = unique(results$dep), minor_breaks = NULL) +
  coord_cartesian(ylim = c(0, 1))


library(ggplot2)
library(tidyr)
library(dplyr)

# Convert to long format for faceting
results_long <- results %>%
  pivot_longer(cols = c(confidence_rate, accuracy), names_to = "metric", values_to = "value")

# Plot
ggplot(results_long, aes(x = factor(dep), y = value)) +
  geom_boxplot(aes(fill = metric), alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  labs(
    title = "Distribution of Disclosure Risk and Accuracy Across Dependency Strengths",
    x = expression(kappa~"(dependency strength)"),
    y = NULL
  ) +
  scale_fill_manual(values = c("confidence_rate" = "steelblue", "accuracy" = "firebrick")) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )


# diagnostics:

estimate_R2s <- function(dep) {
  x <- simulate_microdata(dep = dep, n = 1000)
  
  # health_score R2
  lm_health <- lm(health_score ~ age + income + education, data = x)
  R2_health <- summary(lm_health)$r.squared
  
  # disease_status pseudo R2
  mod_disease <- nnet::multinom(disease_status ~ age + income + education + health_score, data = x, trace = FALSE)
  ll_null <- logLik(update(mod_disease, . ~ 1))
  ll_full <- logLik(mod_disease)
  R2_disease <- 1 - as.numeric(ll_full) / as.numeric(ll_null)
  
  tibble(dep = dep, R2_health = R2_health, R2_disease = R2_disease)
}

dep_values <- seq(0, 5, by = 0.2)
r2_results <- purrr::map_dfr(dep_values, estimate_R2s)

ggplot(r2_results, aes(x = dep)) +
  geom_line(aes(y = R2_health, color = "health_score")) +
  geom_line(aes(y = R2_disease, color = "disease_status")) +
  labs(x = expression(kappa), y = expression(R^2),
       color = "Response variable",
       title = "Explained Variance (R²) Increases with Dependency Strength") +
  theme_minimal(base_size = 13)


# # Generate toy original dataset (truth)
# set.seed(123)
# n <- 100
# truth <- data.frame(
#    age = sample(20:70, n, replace = TRUE),
#    income = round(rnorm(n, mean = 50000, sd = 10000)),
#    education = factor(sample(c("low", "medium", "high"), n, replace = TRUE)),
#    health_score = round(rnorm(n, mean = 75, sd = 10), 1)
# )
# 
# # Generate a synthetic dataset with some variability
# synth <- data.frame(
#   age = truth$age + sample(-2:2, n, replace = TRUE),
#   income = round(truth$income * runif(n, 0.9, 1.1)),
#   education = factor(sample(c("low", "medium", "high"), n, replace = TRUE)),
#   health_score = round(truth$health_score + rnorm(n, 0, 5), 1)
# )


# new: more realistic:

set.seed(123)

n <- 1000

## 1) Age: truncate a normal to 20–70, then round to years
age <- round(pmin(pmax(rnorm(n, mean = 44, sd = 12), 20), 70))

age_s  <- as.numeric(scale(age))           # standardized
age_c  <- (age - 40) / 10                  # centered/scaled for smoother effects

## 2) Education depends on age (younger cohorts more likely to be highly educated)
#    Use a softmax over three logits to get valid class probabilities for each row
l_low   <- -0.2 +  0.8 * age_s             # older -> more "low"
l_med   <-  0.2 +  0.0 * age_s             # baseline
l_high  <- -0.1 -  0.8 * age_s             # younger -> more "high"

logits  <- cbind(l_low, l_med, l_high)
probs   <- exp(logits)
probs   <- probs / rowSums(probs)

draw_one <- function(p) sample(c("low","medium","high"), size = 1, prob = p)
education <- factor(apply(probs, 1, draw_one),
                    levels = c("low","medium","high"), ordered = TRUE)

edu_idx <- as.integer(education) - 1L       # low=0, medium=1, high=2

## 3) Income: log-normal; log-income depends on education and age with a hump (quadratic)
#   - Each step in education lifts log-income
#   - Age shows increasing earnings early, then flatten/decline (− age^2 term)
mu_log   <- 10.7 + 0.22 * edu_idx + 0.12 * age_c - 0.04 * (age_c^2)
sd_log   <- 0.30
income   <- round(rlnorm(n, meanlog = mu_log, sdlog = sd_log))

## 4) Health score: 0–100, depends on age (worse with age),
##    education (better with higher edu), and income (better, diminishing via log)
log_income_std <- as.numeric(scale(log(pmax(income, 1))))  # guard against log(0)
z <- 0.7 - 0.6 * ((age - 50) / 15) + 0.25 * edu_idx + 0.25 * log_income_std + rnorm(n, 0, 0.7)
health_score <- round(100 * plogis(z), 1)  # squash to [0,100]

## Final data frame
truth <- data.frame(
  age = age,
  income = income,
  education = education,
  health_score = health_score
)


library(synthpop)

synth <- syn(truth, m = 1, seed = 123)$syn 

## Make sure factor levels are consistent
# synth$education <- factor(synth$education, levels = levels(truth$education))

# Run the inferential disclosure risk function
result <- risk_inferential(
  truth = truth,
  synth = synth,
  known_vars = c("age", "income", "education"),
  sensitive_var = "health_score",
  model = "rf",
  method = "rmse",                      # root mean squared error
  numeric_threshold_type = "Percentage",
  numeric_threshold = 10,             # 10% relative error threshold
  trace = TRUE
) # TODO: check the warning

# Output
print(result$risk$rows_risk_n)
print(result$risk$rows_risk_p)
head(result$risk$rows_risk_df)


## Categorical sensitive variable:
# # Create toy truth data (original)
# set.seed(2025)
# n <- 100
# truth_cat <- data.frame(
#   age = sample(20:70, n, replace = TRUE),
#   income = round(rnorm(n, 50000, 10000)),
#   education = sample(c("low", "medium", "high"), n, replace = TRUE),
#   disease_status = factor(sample(c("healthy", "diabetic", "hypertensive"), 
#   n, replace = TRUE, prob = c(0.6, 0.2, 0.2)))
# )
# 
# # Create toy synth data (synthetic)
# synth_cat <- truth_cat
# synth_cat$age <- synth_cat$age + sample(-2:2, n, replace = TRUE)
# synth_cat$income <- round(synth_cat$income * runif(n, 0.9, 1.1))
# synth_cat$education <- sample(c("low", "medium", "high"), n, replace = TRUE)
# synth_cat$disease_status <- factor(sample(c("healthy", "diabetic", "hypertensive"), 
#   n, replace = TRUE, prob = c(0.55, 0.25, 0.2)))

# new: more realistic:

add_disease_status <- function(df,
                               health_higher_is_better = TRUE,
                               seed = NULL) {
  stopifnot(all(c("age", "education", "income", "health_score") %in% names(df)))
  if (!is.null(seed)) set.seed(seed)
  
  # helper: stable z-score (handles zero variance)
  z <- function(x) {
    mu <- mean(x, na.rm = TRUE)
    s  <- sd(x,  na.rm = TRUE)
    if (is.na(s) || s == 0) x - mu else (x - mu) / s
  }
  
  # build numeric proxies
  z_age  <- z(df$age)
  z_inc  <- z(log1p(df$income))  # income is typically skewed
  edu_num <- if (is.numeric(df$education)) df$education else as.numeric(as.ordered(df$education))
  z_edu  <- z(edu_num)
  
  # make "higher = worse" health index for easier sign control
  health_bad <- if (health_higher_is_better) -df$health_score else df$health_score
  z_hbad <- z(health_bad)
  
  # Multinomial logits: baseline "healthy" has eta = 0
  # Coefficients chosen for plausible marginal rates and directions:
  #   age (+), worse health (+), income (-), education (-)
  eta_diab <- -2.0 + 0.9 * z_age - 0.3 * z_inc - 0.4 * z_edu + 1.0 * z_hbad
  eta_hyp  <- -1.0 + 1.1 * z_age - 0.2 * z_inc - 0.2 * z_edu + 0.6 * z_hbad
  
  # Softmax to probabilities
  e_d <- exp(eta_diab)
  e_h <- exp(eta_hyp)
  denom <- 1 + e_d + e_h
  
  p_healthy      <- 1 / denom
  p_diabetic     <- e_d / denom
  p_hypertensive <- e_h / denom
  
  probs <- cbind(p_healthy, p_diabetic, p_hypertensive)
  levs  <- c("healthy", "diabetic", "hypertensive")
  
  # Draw one category per row
  draw_one <- function(p) sample(levs, size = 1, prob = p, replace = TRUE)
  status <- apply(probs, 1, draw_one)
  
  df$disease_status <- factor(status, levels = levs)
  
  # (optional) attach quick prevalence check
  attr(df, "disease_status_prevalence") <- setNames(colMeans(probs), levs)
  
  df
}


truth <- add_disease_status(truth, health_higher_is_better = FALSE, seed = 123)
synth <- syn(truth, m = 1, seed = 123)$syn 


# Apply the inferential_risk function with categorical sensitive variable
result_cat <- risk_inferential(
  truth = truth,
  synth = synth,
  known_vars = c("age", "income", "education"),
  sensitive_var = "disease_status",
  model = "rf",
  method = "rmse",
  numeric_threshold_type = "Percentage",
  numeric_threshold = 5,
  categorical_threshold = 1.2
)

# Print results
str(result_cat)

