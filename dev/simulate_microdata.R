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


# simulate_microdata <- function(
#     n = 1000,
#     seed = NULL,
#     # ----- dependency scaling (0 = independent; 1 = default; >1 = stronger) -----
#     dep = 1,
#     # ----- AGE -----
#     age_mean = 45, age_sd = 15, age_min = 18, age_max = 90,
#     # ----- EDUCATION (ordered: low < medium < high); latent ~ N(mu, sd^2) -----
#     edu_age_beta = -0.15,      # younger cohorts more likely high ed (negative wrt age_z)
#     edu_latent_sd = 1.0,
#     edu_cutpoints = c(-0.3, 0.7),  # cutpoints for low|med|high on the latent scale
#     # ----- INCOME (log-income model) -----
#     log_income_mean = log(40000),  # baseline mean income around 40k
#     income_beta_age = 0.30,        # effect of age_z
#     income_beta_age2 = -0.15,      # life-cycle concavity
#     income_beta_edu = 0.35,        # per education step (low=0, med=1, high=2)
#     income_noise_sd = 0.35,        # residual sd on log-income
#     # ----- HEALTH SCORE (0..100; higher = better) -----
#     health_mean = 70,
#     health_beta_age = -8.0,        # per SD of age (z-scale)
#     health_beta_log_income = 6.0,  # per SD of log-income (z-scale)
#     health_beta_edu = 4.0,         # per education step
#     health_noise_sd = 10.0,
#     # ----- DISEASE STATUS (multinomial: healthy, diabetic, hypertensive) -----
#     # Baseline class = "healthy"; intercepts target ~70/15/15 when predictors at mean
#     disease_intercepts = c(diabetic = -1.54, hypertensive = -1.54),
#     disease_beta_age = c(diabetic = 0.80, hypertensive = 1.00),       # age_z effects
#     disease_beta_log_income = c(diabetic = -0.30, hypertensive = -0.20),# log_inc_z
#     disease_beta_health = c(diabetic = -1.20, hypertensive = -1.00),    # (health - mean)/15
#     disease_beta_edu = c(diabetic = -0.20, hypertensive = -0.10),       # per edu step
#     # ----- miscellaneous -----
#     return_numeric_edu = FALSE
# ) {
#   if (!is.null(seed)) set.seed(seed)
#
#   # --- helpers ---
#   clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))
#   rtruncnorm_base <- function(n, mean, sd, lo, hi) {
#     x <- rnorm(n, mean, sd)
#     oob <- which(x < lo | x > hi)
#     while (length(oob) > 0) {
#       x[oob] <- rnorm(length(oob), mean, sd)
#       oob <- which(x < lo | x > hi)
#     }
#     x
#   }
#   draw_multinomial <- function(p_mat, labels) {
#     # p_mat: n x K matrix of probabilities
#     u <- runif(nrow(p_mat))
#     cs <- t(apply(p_mat, 1, cumsum))
#     idx <- rowSums(u > cs) + 1L
#     factor(labels[idx], levels = labels)
#   }
#
#   # 1) AGE ----------------------------------------------------------------------
#   age <- rtruncnorm_base(n, age_mean, age_sd, age_min, age_max)
#   age_z <- (age - age_mean) / age_sd
#
#   # 2) EDUCATION (ordinal via latent variable + cutpoints) ----------------------
#   # latent mean depends on age_z; sd = edu_latent_sd
#   edu_latent_mu <- dep * (edu_age_beta * age_z)
#   edu_latent <- rnorm(n, mean = edu_latent_mu, sd = edu_latent_sd)
#
#   # cut into low/medium/high
#   cut1 <- edu_cutpoints[1]; cut2 <- edu_cutpoints[2]
#   edu_num <- ifelse(edu_latent < cut1, 0L,
#                     ifelse(edu_latent < cut2, 1L, 2L))
#   education <- factor(c("low","medium","high")[edu_num + 1L],
#                       levels = c("low","medium","high"))
#
#   # 3) INCOME (log-linear; depends on age_z, age_z^2, education) ---------------
#   log_inc <- log_income_mean +
#     dep * (income_beta_age  * age_z +
#              income_beta_age2 * (age_z^2) +
#              income_beta_edu  * edu_num) +
#     rnorm(n, sd = income_noise_sd)
#
#   income <- as.numeric(exp(log_inc))  # positively skewed, realistic
#   # standardize log-income for downstream models
#   log_inc_z <- (log_inc - log_income_mean) / income_noise_sd
#
#   # 4) HEALTH SCORE (0..100; higher = better) ----------------------------------
#   # Use z-scales for age and log-income; education as steps (0,1,2)
#   health_raw <- health_mean +
#     dep * (health_beta_age        * age_z +
#              health_beta_log_income * log_inc_z +
#              health_beta_edu        * edu_num) +
#     rnorm(n, sd = health_noise_sd)
#
#   health_score <- clamp(health_raw, 0, 100)
#   # modest scale for health in disease logits
#   health_z <- (health_score - health_mean) / 15
#
#   # 5) DISEASE STATUS (multinomial logistic vs "healthy") ----------------------
#   # linear predictors for diabetic & hypertensive relative to healthy
#   eta_diab <- disease_intercepts["diabetic"] +
#     dep * (disease_beta_age["diabetic"]         * age_z +
#              disease_beta_log_income["diabetic"]  * log_inc_z +
#              disease_beta_health["diabetic"]      * health_z +
#              disease_beta_edu["diabetic"]         * edu_num)
#
#   eta_hyp  <- disease_intercepts["hypertensive"] +
#     dep * (disease_beta_age["hypertensive"]        * age_z +
#              disease_beta_log_income["hypertensive"] * log_inc_z +
#              disease_beta_health["hypertensive"]     * health_z +
#              disease_beta_edu["hypertensive"]        * edu_num)
#
#   # softmax with healthy as baseline (eta_healthy = 0)
#   num_diab <- exp(eta_diab)
#   num_hyp  <- exp(eta_hyp)
#   denom    <- 1 + num_diab + num_hyp
#   p_healthy <- 1 / denom
#   p_diab    <- num_diab / denom
#   p_hyp     <- num_hyp  / denom
#
#   probs <- cbind(healthy = p_healthy,
#                  diabetic = p_diab,
#                  hypertensive = p_hyp)
#
#   disease_status <- draw_multinomial(probs, labels = c("healthy","diabetic","hypertensive"))
#
#   # gender
#
#   # Gender: conditional Bernoulli with logistic model
#   # Adjusted coefficients based on variable scales and desired balance
#   eta_gender <- -3 +
#     0.03 * (age - 45) +
#     0.6  * edu_num +
#     0.00004 * income -
#     0.02 * (health_score - 70)
#
#   prob_male <- 1 / (1 + exp(-eta_gender))
#   gender <- factor(rbinom(n, 1, prob = prob_male), labels = c("female", "male"))
#
#   # OUTPUT ----------------------------------------------------------------------
#   out <- data.frame(
#     gender = as.factor(gender),
#     age = age,
#     education = as.factor(education),
#     income = income,
#     health_score = health_score,
#     disease_status = as.factor(disease_status)
#   )
#   if (return_numeric_edu) out$education_num <- edu_num
#   out
# }
